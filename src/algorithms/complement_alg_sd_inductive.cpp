// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_sd_inductive.hpp"

#include <cassert>
#include <map>
#include <stdexcept>
#include <type_traits>

#include "../util/helpers.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

namespace kofola {
namespace sd_inductive {

namespace {

using base_tree = kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf>;

check_macrostate assign_leaf_ids(const check_macrostate& tree, unsigned& next_id) {
  if (tree.is_leaf()) {
    return std::visit(
      [&](const auto& leaf) -> check_macrostate {
        using leaf_t = std::decay_t<decltype(leaf)>;
        const unsigned id = next_id++;
        if constexpr (std::is_same_v<leaf_t, fin_leaf>) {
          return check_macrostate::fin(tree.get_options_ptr(), leaf.safe, leaf.color, id);
        } else {
          static_assert(std::is_same_v<leaf_t, inf_leaf>, "Unexpected leaf type");
          return check_macrostate::inf(tree.get_options_ptr(), leaf.track, leaf.breakpoint, leaf.color, id);
        }
      },
      tree.leaf_value());
  }

  const check_macrostate left_ms(tree.get_options_ptr(), base_tree(tree.left()));
  const check_macrostate right_ms(tree.get_options_ptr(), base_tree(tree.right()));
  const auto node_type = tree.type();
  return check_macrostate::make(tree.get_options_ptr(), node_type, assign_leaf_ids(left_ms, next_id), assign_leaf_ids(right_ms, next_id), tree.node_value().context);
}

/**
 * @brief Initialize (or re-initialize) `NodeContext` payloads inside a check tree.
 *
 * When shared-breakpoint mode is enabled (`opts->use_shared_breakpoint == true`),
 * each internal node's `NodeContext` is recomputed from the subtree so that it
 * references the current set of `Inf` leaf IDs.
 *
 * @param tree Check tree to update in-place.
 * @param opts Options controlling whether shared-breakpoint contexts are used.
 *             Must be non-null.
 *
 * @warning This function traverses children by casting `base_tree::left()/right()`
 *          to `check_macrostate&`. The underlying tree stores children as
 *          `base_tree`, not `check_macrostate`, so these casts are undefined
 *          behavior and can manifest as crashes (e.g., seemingly-null
 *          `get_options_ptr()`). Prefer a functional rebuild that wraps children
 *          as `check_macrostate(opts, base_tree(child))` if you need a safe
 *          traversal.
 */
void init_contexts_in_tree(check_macrostate& tree, options_ptr opts) {
  if (tree.is_leaf()) {
    return;
  }

  auto& base = static_cast<base_tree&>(tree);
  init_contexts_in_tree(static_cast<check_macrostate&>(base.left()), opts);
  init_contexts_in_tree(static_cast<check_macrostate&>(base.right()), opts);

  if (opts->use_shared_breakpoint) {
    NodeContext ctx = NodeContext::create_subtree_sh_context(tree.type(), tree);
    tree.node_value().set_context(ctx);
  }
  
}

/**
 * Return a copy of `tree` with all states from `forbidden` removed
 * from any leaf state-sets.
 */
static sd_inductive::check_macrostate restrict_states_in_tree(
  const sd_inductive::check_macrostate& tree,
  const std::set<unsigned>& forbidden) {

  using check_macrostate = sd_inductive::check_macrostate;
  using base_tree = kofola::types::binary_tree<sd_inductive::TreeType, sd_inductive::AndOrNode, sd_inductive::fin_leaf, sd_inductive::inf_leaf>;

  if (tree.is_leaf()) {
    return std::visit(
      [&](const auto& leaf) -> check_macrostate {
        using leaf_t = std::decay_t<decltype(leaf)>;
        if constexpr (std::is_same_v<leaf_t, sd_inductive::fin_leaf>) {
          return check_macrostate::fin(tree.get_options_ptr(), get_set_difference(leaf.safe, forbidden), leaf.color, leaf.id);
        } else {
          static_assert(std::is_same_v<leaf_t, sd_inductive::inf_leaf>, "Unexpected leaf type");
          return check_macrostate::inf(
            tree.get_options_ptr(),
            get_set_difference(leaf.track, forbidden),
            get_set_difference(leaf.breakpoint, forbidden),
            leaf.color,
            leaf.id);
        }
      },
      tree.leaf_value());
  }

  const auto node_type = tree.type();
  check_macrostate left(tree.get_options_ptr(), base_tree(tree.left()));
  check_macrostate right(tree.get_options_ptr(), base_tree(tree.right()));
  NodeContext context = tree.node_value().context;
  context.restrict_states(forbidden);
  return check_macrostate::make(
    tree.get_options_ptr(),
    node_type,
    restrict_states_in_tree(left, forbidden),
    restrict_states_in_tree(right, forbidden),
    context);
}

/**
 * Collect all top-level OR disjuncts by flattening consecutive OR nodes
 * into @p parts.  Non-OR nodes (including leaves and AND) are appended as-is.
 */
static void collect_or_disjuncts(const check_macrostate& tree, std::vector<check_macrostate>& parts) {
  if (tree.is_leaf() || tree.type() != TreeType::Or) {
    parts.push_back(tree);
    return;
  }
  collect_or_disjuncts(check_macrostate(tree.get_options_ptr(), base_tree(tree.left())), parts);
  collect_or_disjuncts(check_macrostate(tree.get_options_ptr(), base_tree(tree.right())), parts);
}

/**
 * Build a left-associative OR chain from the non-empty subrange
 * @p parts[@p begin .. @p end).
 */
static check_macrostate build_or_chain(options_ptr opts,std::vector<check_macrostate>& parts, size_t begin, size_t end) {
  check_macrostate result = parts[begin];
  for (size_t i = begin + 1; i < end; ++i) {
    result = check_macrostate::make(opts, TreeType::Or, std::move(result), parts[i]);
  }
  return result;
}

/**
 * Reorganize a check tree so that in every OR node all subtrees that consist
 * only of Fin leaves (satisfying `has_only_fin_leaves_no_inner_or`) are
 * placed in the left subtree, and all remaining (Inf-containing) subtrees
 * are placed in the right subtree.  AND nodes are left structurally
 * unchanged; their children are reorganized recursively.
 *
 * This is a prerequisite for the OR-FIN optimization in
 * `check_macrostate::get_succ`, which assumes that the left subtree of
 * every Or node satisfies `has_only_fin_leaves_no_inner_or`.
 *
 * @param tree Tree to reorganize.
 * @return A semantically equivalent tree with FIN-only disjuncts on the
 *         left of every OR node.
 */
static check_macrostate reorganize_fins_left(const check_macrostate& tree) {
  if (tree.is_leaf()) return tree;

  const auto node_type = tree.type();
  const check_macrostate left_ms(tree.get_options_ptr(), base_tree(tree.left()));
  const check_macrostate right_ms(tree.get_options_ptr(), base_tree(tree.right()));

  if (node_type == TreeType::And) {
    return check_macrostate::make(
        tree.get_options_ptr(), TreeType::And,
        reorganize_fins_left(left_ms), reorganize_fins_left(right_ms),
        tree.node_value().context);
  }

  // OR node: flatten all consecutive OR disjuncts, reorganize each
  // recursively, bucket into FIN-only vs other, then rebuild.
  std::vector<check_macrostate> parts;
  collect_or_disjuncts(tree, parts);

  std::vector<check_macrostate> fin_parts, other_parts;
  for (const auto& p : parts) {
    auto reorg = reorganize_fins_left(p);
    if (has_only_fin_leaves_no_inner_or(reorg)) {
      fin_parts.push_back(std::move(reorg));
    } else {
      other_parts.push_back(std::move(reorg));
    }
  }

  if (fin_parts.empty()) {
    return build_or_chain(tree.get_options_ptr(), other_parts, 0, other_parts.size());
  }

  auto fin_chain = build_or_chain(tree.get_options_ptr(), fin_parts, 0, fin_parts.size());

  if (other_parts.empty()) {
    return fin_chain;
  }

  auto other_chain = build_or_chain(tree.get_options_ptr(), other_parts, 0, other_parts.size());
  return check_macrostate::make(
      tree.get_options_ptr(), TreeType::Or,
      std::move(fin_chain), std::move(other_chain));
}

} // namespace

/**
 * Parse a Spot acceptance condition and build an equivalent
 * `check_macrostate` tree.
 *
 * The function supports the following acceptance grammar (no negation
 * expected):
 * - conjunction: `&`  -> `TreeType::And`
 * - disjunction: `|`  -> `TreeType::Or`
 * - leaves: `Fin(m)`  -> `TreeType::Fin` (mark `m` stored in
 *   `fin_leaf::color`)
 * - leaves: `Inf(m)`  -> `TreeType::Inf` (mark `m` stored in
 *   `inf_leaf::color`)
 *
 * This method is a thin public wrapper that delegates the actual
 * parsing work to `from_acc_code_impl`.
 *
 * @param code The Spot acceptance formula to convert.
 * @return A `check_macrostate` representing the parsed acceptance
 *         formula.
 * @throws std::invalid_argument if `code` is empty or contains an
 *         unsupported/top-level operator that cannot be represented.
 */
check_macrostate check_macrostate::from_acc_code(options_ptr opts, const spot::acc_cond::acc_code& code) {
  auto parsed = from_acc_code_impl(opts, code);
  if (opts && opts->use_or_fin_opt) {
    parsed = reorganize_fins_left(parsed);
  }
  unsigned next_id = 0;
  // First assign stable leaf IDs, then (re)initialize NodeContext using the
  // fully ID-annotated subtrees.
  return assign_leaf_ids(parsed, next_id).init_contexts();
}

/**
 * @brief Return a copy of this check tree with refreshed `NodeContext` values.
 *
 * Intended usage is after parsing acceptance and assigning stable leaf IDs.
 * In shared-breakpoint mode, `NodeContext` depends on `inf_leaf::id`, so any
 * transformation that changes leaf IDs should be followed by `init_contexts()`.
 *
 * @note The current implementation delegates to `init_contexts_in_tree()`,
 *       which performs an in-place traversal using casts that assume children
 *       are `check_macrostate`. If you observe intermittent crashes here,
 *       refactor to a traversal that rebuilds the tree while wrapping children
 *       with `check_macrostate(opts, base_tree(child))`. 
 */
check_macrostate check_macrostate::init_contexts() const {
  check_macrostate out = *this;
  init_contexts_in_tree(out, this->opts_);
  return out;
}

/**
 * Compute successor macrostate(s) for this `check_macrostate` node.
 *
 * Behavior:
 * - Leaf: delegate to the leaf's `get_succ()` implementation.
 * - And: compute successors of both children using the same
 *   `check_states` and combine results via the Cartesian product.
 * - Or: nondeterministically partition `check_states` into two parts,
 *   compute successors for each child using the corresponding partition,
 *   then combine and deduplicate the results.
 *
 * @param aut Pointer to the Spot automaton graph.
 * @param scc_info SCC information for the automaton (used to restrict
 *                 transitions to the current SCC).
 * @param check_states Set of automaton states currently being checked/tracked.
 * @param bdd BDD representing the current input/condition used to test
 *            transition guards.
 * @return A vector of successor `check_macrostate` instances produced by
 *         advancing this check node under the given automaton/transitions.
 * @throws std::logic_error If an unexpected internal node type is
 *         encountered or if `nondet_split_set` yields a non-binary
 *         partition.
 */
std::vector<std::pair<check_macrostate, NodeContext>> check_macrostate::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  const bdd&                        bdd,
  bool                              resample,
  NodeContext                       parent_context) const {

  if (this->is_leaf()) {
    return std::visit(
      [&](const auto& leaf) {
        return leaf.get_succ(aut, scc_info, check_states, this->opts_, bdd, resample, parent_context);
      },
      this->leaf_value());
  }

  NodeContext context_sent = parent_context;
  NodeContext actual_node_context = this->node_value().get_context();
  bool is_scope_root = actual_node_context.is_scope_root(parent_context);
  if (is_scope_root) {
    context_sent = this->node_value().get_succ_context(resample);
    actual_node_context = context_sent;
  }

  const auto node_type = this->type();
  const check_macrostate left_ms(this->opts_, base_tree(this->left()));
  const check_macrostate right_ms(this->opts_, base_tree(this->right()));

  if (node_type == TreeType::And) {
    const auto left_succ = left_ms.get_succ(aut, scc_info, check_states, bdd, resample, context_sent);
    if (left_succ.empty()) return {};
    const auto right_succ = right_ms.get_succ(aut, scc_info, check_states, bdd, resample, context_sent);

    return cartesian_product<std::pair<check_macrostate, NodeContext>, std::pair<check_macrostate, NodeContext>>(
      left_succ,
      right_succ,
      [opts = this->opts_, is_scope_root, &actual_node_context](const std::pair<check_macrostate, NodeContext>& l, const std::pair<check_macrostate, NodeContext>& r) -> std::pair<check_macrostate, NodeContext> {
        NodeContext local = actual_node_context;
        NodeContext merge = l.second.union_contexts(r.second);
        if (is_scope_root) {
          auto viol_states = merge.violating_states;
          auto viol_preds  = merge.violating_predecessors;
          if (!merge.is_none()) local = merge;
          merge = NodeContext{};
          merge.violating_states       = std::move(viol_states);
          merge.violating_predecessors = std::move(viol_preds);
        }
        auto tmp = check_macrostate::make(opts, TreeType::And, l.first, r.first, local);
        return {tmp, merge};
      });
  }

  if (node_type == TreeType::Or) {
    // OR-FIN optimization: when the left subtree contains only FIN leaves (no inner Or
    // nodes), send ALL check_states to the left subtree with collect_violating=true.
    // States that fire a Fin-colored transition are collected as "violating" and are
    // subsequently removed from the left subtree result and forwarded to the right
    // subtree as check_states. If the right subtree itself cannot accommodate those
    // states (returns empty), this path produces no successor. Any violating states
    // that propagate back from the right subtree are passed up to our caller.
    if (this->opts_ && this->opts_->use_or_fin_opt && has_only_fin_leaves_no_inner_or(left_ms)) {
      NodeContext left_ctx = context_sent;
      left_ctx.collect_violating = true;

      const auto left_succ = left_ms.get_succ(aut, scc_info, check_states, bdd, resample, left_ctx);

      std::set<std::pair<check_macrostate, NodeContext>> out;
      for (const auto& [left_tree, left_node_ctx] : left_succ) {
        const std::set<unsigned>& viol_succs = left_node_ctx.violating_states;
        const std::set<unsigned>& viol_preds = left_node_ctx.violating_predecessors;

        // Remove successors of violating states from the left subtree (left_tree has successor states).
        const check_macrostate left_restricted =
          viol_succs.empty() ? left_tree : restrict_states_in_tree(left_tree, viol_succs);

        // Route violating predecessor states to the right subtree as fresh check_states.
        const auto right_succ = right_ms.get_succ(aut, scc_info, viol_preds, bdd, resample, context_sent);
        if (right_succ.empty()) {
          // Right cannot accommodate the violating states: no valid successor for this case.
          continue;
        }

        for (const auto& [right_tree, right_node_ctx] : right_succ) {
          // Propagate any further violations from right up to the parent.
          NodeContext out_ctx{};
          out_ctx.violating_states = right_node_ctx.violating_states;
          out_ctx.violating_predecessors = right_node_ctx.violating_predecessors;
          auto combined = check_macrostate::make(this->opts_, TreeType::Or,
            left_restricted, right_tree, actual_node_context);
          out.insert({combined.reduce(), out_ctx});
        }
      }
      return std::vector<std::pair<check_macrostate, NodeContext>>(out.begin(), out.end());
    }

    // Optimization: reduce check_states to behavior-equivalent representatives.
    // In a deterministic SCC, each state has at most one successor per symbol.
    // States with the same (successor_state, transition_marks) pair produce
    // identical effects on all leaves, so only one representative per group
    // is needed. States with no SCC-internal successor are dropped entirely.
    using sig_t = std::pair<unsigned, spot::acc_cond::mark_t>;
    std::map<sig_t, unsigned> sig_to_repr;
    std::set<unsigned> repr_states;

    for (unsigned s : check_states) {
      for (const auto& t : aut->out(s)) {
        if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
          sig_t sig{t.dst, t.acc};
          auto [it, inserted] = sig_to_repr.try_emplace(sig, s);
          if (inserted) {
            repr_states.insert(s);
          }
          break; // deterministic SCC: at most one successor per symbol
        }
      }
    }

    std::set<std::pair<check_macrostate, NodeContext>> out;
    const auto partitions = nondet_split_set(repr_states, 2);
    for (const auto& part : partitions) {
      if (part.size() != 2) {
        throw std::logic_error("check_macrostate::get_succ: nondet_split_set returned non-binary partition");
      }
      const auto& left_states = part[0];
      const auto& right_states = part[1];
      const auto left_succ = left_ms.get_succ(aut, scc_info, left_states, bdd, resample, context_sent);
      if (left_succ.empty()) continue;
      const auto right_succ = right_ms.get_succ(aut, scc_info, right_states, bdd, resample, context_sent);
      if (right_succ.empty()) continue;

      // TODO: it is not efficient to call reduce here
      const auto combined = cartesian_product<std::pair<check_macrostate, NodeContext>, std::pair<check_macrostate, NodeContext>>(
        left_succ,
        right_succ,
        [opts = this->opts_, is_scope_root, &actual_node_context](const std::pair<check_macrostate, NodeContext>& l, const std::pair<check_macrostate, NodeContext>& r) -> std::pair<check_macrostate, NodeContext> {
          NodeContext local = actual_node_context;
          NodeContext merge = l.second.union_contexts(r.second);
          if (is_scope_root) {
            if (!merge.is_none()) local = merge;
            merge = NodeContext{};
          }
          auto tmp = check_macrostate::make(opts, TreeType::Or, l.first, r.first, local);
          return { tmp.reduce(), merge };
        });
      out.insert(combined.begin(), combined.end());
    }
    return std::vector<std::pair<check_macrostate, NodeContext>>(out.begin(), out.end());
  }

  throw std::logic_error("check_macrostate::get_succ: unexpected internal node type");
}

/**
 * Determine whether this `check_macrostate` tree is satisfied.
 *
 * For leaves this delegates to the leaf's `is_satisfied()` method. For
 * internal nodes both children are required to be satisfied (this mirrors
 * the structural check used by the complement construction).
 *
 * @return `true` if the entire check tree is satisfied, `false` otherwise.
 * @throws std::logic_error If the node contains an unexpected internal
 *         type.
 */
bool check_macrostate::is_satisfied() const {
  if (this->is_leaf()) {
    return std::visit(
      [&](const auto& leaf) {
        return leaf.is_satisfied();
      },
      this->leaf_value());
  }

  const auto node_type = this->type();
  const check_macrostate left_ms(this->opts_, base_tree(this->left()));
  const check_macrostate right_ms(this->opts_, base_tree(this->right()));

  if (node_type == TreeType::And || node_type == TreeType::Or) {
    return this->node_value().is_satisfied() && left_ms.is_satisfied() && right_ms.is_satisfied();
  }

  throw std::logic_error("check_macrostate::is_satisfied: unexpected internal node type");
}

/**
 * Gather all automaton states referenced by this `check_macrostate`.
 *
 * For a leaf node:
 * - `Fin` leaf: returns the `safe` set.
 * - `Inf` leaf: returns the `track` set.
 *
 * For an internal `And` or `Or` node: returns the union of the
 * states gathered from both children.
 *
 * @return A `std::set<unsigned>` containing the state indices referenced
 *         in this check tree.
 * @throws std::logic_error if the node contains an unexpected internal
 *         type (neither `And` nor `Or`).
 */
std::set<unsigned> check_macrostate::gather_states() const {
  if (this->is_leaf()) {
    return std::visit(
      [](const auto& leaf) -> std::set<unsigned> {
        using leaf_t = std::decay_t<decltype(leaf)>;
        if constexpr (std::is_same_v<leaf_t, fin_leaf>) {
          return leaf.safe;
        } else {
          static_assert(std::is_same_v<leaf_t, inf_leaf>, "Unexpected leaf type");
          return leaf.track;
        }
      },
      this->leaf_value());
  }

  const auto node_type = this->type();
  if (node_type != TreeType::And && node_type != TreeType::Or) {
    throw std::logic_error("check_macrostate::gather_states: unexpected internal node type");
  }

  const check_macrostate left_ms(this->opts_, base_tree(this->left()));
  const check_macrostate right_ms(this->opts_, base_tree(this->right()));
  return get_set_union(left_ms.gather_states(), right_ms.gather_states());
}

/**
 * Reduce/simplify this `check_macrostate` tree.
 *
 * Reduction rules applied:
 * - Leaf nodes are returned unchanged.
 * - `And` nodes: both children are reduced and the `And` node is
 *   reconstructed from the reduced children.
 * - `Or` nodes: the left child is reduced, then the set of states
 *   referenced by the reduced left child is gathered and removed from
 *   the right child (via `restrict_states_in_tree`). The right child
 *   (after restriction) is reduced and the `Or` node is rebuilt from
 *   the reduced left and reduced/right-restricted right child.
 *
 * This transformation ensures that states mentioned on the left of an
 * `Or` do not persist in the right subtree, avoiding redundant
 * checking and improving determinism of the check tree.
 *
 * @return A reduced `check_macrostate` equivalent to this tree but with
 *         redundant state references removed where applicable.
 * @throws std::logic_error if the node contains an unexpected internal
 *         type.
 */
check_macrostate check_macrostate::reduce() const {
  if (this->is_leaf()) {
    return *this;
  }

  const auto node_type = this->type();
  const check_macrostate left_ms(this->opts_, base_tree(this->left()));
  const check_macrostate right_ms(this->opts_, base_tree(this->right()));

  if (node_type == TreeType::And) {
    return check_macrostate::make(this->opts_, TreeType::And, left_ms.reduce(), right_ms.reduce(), this->node_value().context);
  }

  if (node_type == TreeType::Or) {
    // OR-FIN optimization: the right subtree holds states that were moved there
    // because they fired a Fin-colored transition (i.e., they violated the left/Fin
    // condition). These states must remain in the right subtree so that the Inf check
    // can track them. Priority is therefore reversed: remove from the *left* any
    // states already present in the *right*, not the other way around.
    if (this->opts_ && this->opts_->use_or_fin_opt) {
      const check_macrostate right_red = right_ms.reduce();
      const auto right_states = right_red.gather_states();
      const check_macrostate left_restricted = restrict_states_in_tree(left_ms, right_states);
      return check_macrostate::make(this->opts_, TreeType::Or, left_restricted.reduce(), right_red, this->node_value().context);
    }

    const check_macrostate left_red = left_ms.reduce();
    const auto left_states = left_red.gather_states();
    const check_macrostate right_restricted = restrict_states_in_tree(right_ms, left_states);
    return check_macrostate::make(this->opts_, TreeType::Or, left_red, right_restricted.reduce(), this->node_value().context);
  }

  throw std::logic_error("check_macrostate::reduce: unexpected internal node type");
}

namespace {

/**
 * Recursively collect all `inf_leaf` payloads from the tree in in-order traversal.
 *
 * @param tree The check tree (or subtree) to traverse.
 * @param out  Output vector to append leaf copies into.
 */
void collect_inf_leaves_inorder(const check_macrostate& tree, std::vector<inf_leaf>& out) {
  if (tree.is_leaf()) {
    if (tree.type() == TreeType::Inf) {
      out.push_back(std::get<inf_leaf>(tree.leaf_value()));
    }
    return;
  }

  const check_macrostate left_ms(tree.get_options_ptr(), base_tree(tree.left()));
  const check_macrostate right_ms(tree.get_options_ptr(), base_tree(tree.right()));

  collect_inf_leaves_inorder(left_ms, out);
  collect_inf_leaves_inorder(right_ms, out);
}

} // namespace

std::optional<unsigned> check_macrostate::find_first_inf_leaf() const {
  std::vector<inf_leaf> inf_leaves;
  collect_inf_leaves_inorder(*this, inf_leaves);

  if (inf_leaves.empty()) {
    return std::nullopt;
  }

  return inf_leaves.front().id;
}

std::optional<unsigned> check_macrostate::find_next_inf_leaf(unsigned current_id) const {  
  std::vector<inf_leaf> inf_leaves;
  collect_inf_leaves_inorder(*this, inf_leaves);

  if (inf_leaves.empty()) {
    return std::nullopt;
  }

  // Find the current leaf in the collected list
  int current_idx = -1;
  for (size_t i = 0; i < inf_leaves.size(); ++i) {
    if (inf_leaves[i].id == current_id) {
      current_idx = i;
      break;
    }
  }

  // If current not found, return first leaf
  if (current_idx == -1) {
    for (const auto& leaf : inf_leaves) {
      if (!leaf.track.empty()) {
        return leaf.id;
      }
    }
    // No leaf with nonempty track found, return first leaf anyway
    return inf_leaves.front().id;
  }

  // Look for next leaf after current with nonempty track set
  for (size_t j = current_idx + 1; j < inf_leaves.size(); ++j) {
    if (!inf_leaves[j].track.empty()) {
      return inf_leaves[j].id;
    }
  }

  // No next leaf with nonempty track set found; wrap around and return first leaf
  for (const auto& leaf : inf_leaves) {
    if (!leaf.track.empty()) {
      return leaf.id;
    }
  }

  // No leaf with nonempty track found, return first leaf
  return inf_leaves.front().id;
}

std::optional<inf_leaf> check_macrostate::find_inf_leaf_by_id(unsigned id) const {
  std::vector<inf_leaf> inf_leaves;
  collect_inf_leaves_inorder(*this, inf_leaves);
  
  for (const auto& leaf : inf_leaves) {
    if (leaf.id == id) {
      return leaf;
    }
  }
  return std::nullopt;
}

/**
 * Fold a sequence of acceptance code parts into a single
 * `check_macrostate` tree using the specified binary operator.
 *
 * The function constructs a left-associated tree by converting the
 * first element with `from_acc_code_impl` and successively combining
 * remaining parts with `check_macrostate::make(op, left, right)`.
 *
 * @param op The binary tree operator (`TreeType::And` or `TreeType::Or`)
 *           used to combine the parts.
 * @param parts Vector of Spot acceptance `acc_code` parts to fold.
 * @return A `check_macrostate` representing the folded acceptance
 *         formula.
 * @throws std::invalid_argument if `parts` is empty.
 */
check_macrostate check_macrostate::fold(options_ptr opts, TreeType op, const std::vector<spot::acc_cond::acc_code>& parts) {
  if (parts.empty()) {
    throw std::invalid_argument("check_macrostate: empty And/Or in acceptance formula");
  }
  check_macrostate acc = from_acc_code_impl(opts, parts.front());
  for (size_t i = 1; i < parts.size(); ++i) {
    acc = check_macrostate::make(opts, op, std::move(acc), from_acc_code_impl(opts, parts[i]));
  }
  return acc;
}

/**
 * Internal parser that builds a `check_macrostate` from a Spot
 * `acc_code` acceptance formula. Handles leaf forms (`Fin`/`Inf`),
 * and flattens top-level conjunctions/disjunctions via
 * `top_conjuncts()` / `top_disjuncts()`.
 *
 * The function attempts a single-level unwrap of parenthesised
 * singletons and delegates to `fold()` for multi-operand operators.
 *
 * @param code The Spot acceptance code to parse.
 * @return A `check_macrostate` representing the parsed acceptance.
 * @throws std::invalid_argument if `code` is empty or contains an
 *         unsupported acceptance operator.
 */
check_macrostate check_macrostate::from_acc_code_impl(options_ptr opts, const spot::acc_cond::acc_code& code) {
  if (code.empty()) {
    throw std::invalid_argument("check_macrostate: empty acceptance formula");
  }

  // Prefer top-level flattening (Spot returns a singleton vector when the operator is not present at top-level).
  const auto conjuncts = code.top_conjuncts();
  if (conjuncts.size() > 1) {
    return fold(opts, TreeType::And, conjuncts);
  }
  const auto disjuncts = code.top_disjuncts();
  if (disjuncts.size() > 1) {
    return fold(opts, TreeType::Or, disjuncts);
  }

  if (code.size() == 2) {
    const auto op = code[1].sub.op;
    const auto mark = code[0].mark;
    if (op == spot::acc_cond::acc_op::Fin) {
      return check_macrostate(opts, base_tree::leaf(TreeType::Fin, fin_leaf{{}, mark, 0}));
    }
    if (op == spot::acc_cond::acc_op::Inf) {
      return check_macrostate(opts, base_tree::leaf(TreeType::Inf, inf_leaf{{}, {}, mark, 0}));
    }
  }

  throw std::invalid_argument("check_macrostate: unsupported acceptance formula operator");
}

std::string check_macrostate::tree_type_to_string(TreeType t) {
  switch (t) {
    case TreeType::Fin:
      return "Fin";
    case TreeType::Inf:
      return "Inf";
    case TreeType::And:
      return "And";
    case TreeType::Or:
      return "Or";
  }
  return "?";
}

std::string check_macrostate::to_string_impl(const base_tree& tree, bool show_shared_breakpoint) {
  const std::string head = tree_type_to_string(tree.type());
  if (tree.is_leaf()) {
    const std::string payload = std::visit(
      [](const auto& leaf) { return leaf.to_string(); },
      tree.leaf_value());
    return head + "(" + payload + ")";
  }

  std::string extra;
  if (show_shared_breakpoint && (tree.type() == TreeType::And || tree.type() == TreeType::Or)) {
    extra = "[sb=" + tree.node_value().context.to_string() + "]";
  }

  return head + extra + "(" + to_string_impl(tree.left(), show_shared_breakpoint) + ", " + to_string_impl(tree.right(), show_shared_breakpoint) + ")";
}

/**
 * Compute successor macrostate(s) for a `Fin` leaf.
 *
 * The method builds the union of the leaf's `safe` set and the
 * incoming `check_states`, then explores outgoing transitions from
 * each state in that set. Only transitions that remain within the
 * same SCC and whose guard is implied by `bdd` are considered.
 * If any considered transition carries an accepting mark that
 * intersects this leaf's `color`, the function returns an empty
 * vector (no valid successors). Otherwise the set of destination
 * states is collected and returned as a single `check_macrostate::fin`.
 *
 * @param aut Pointer to the Spot automaton graph.
 * @param scc_info SCC information for the automaton (used to restrict
 *                 transitions to the current SCC).
 * @param check_states Set of automaton states currently being checked/tracked.
 * @param bdd BDD representing the current input/condition used to test
 *            transition guards.
 * @return Vector containing a single `check_macrostate::fin` with the
 *         successor states, or an empty vector if an accepting transition
 *         matching `color` is encountered (no successors).
 */
std::vector<std::pair<check_macrostate, NodeContext>> fin_leaf::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  options_ptr                        opts,
  const bdd&                        bdd,
  bool                              resample,
  NodeContext                       context) const {

  (void)resample; // unused
  std::set<unsigned> st = get_set_union(this->safe, check_states);
  std::set<unsigned> succs {};

  // OR-FIN optimization: when collecting violations, inspect both the leaf's
  // `safe` set and the incoming `check_states`. Any source state that fires a
  // Fin-colored transition is recorded: its successor targets are placed into
  // `violating_succs` and the source is recorded in `violating_predecessors`.
  // The leaf returns the non-violating successor set while reporting these
  // violating sets so the caller (an Or-node) can remove violating successors
  // from the left subtree and route the violating predecessor states to the
  // right subtree for further processing. States in `safe` are not treated
  // specially here; they are reported like any other source and handled by
  // the caller.
  if (opts && opts->use_or_fin_opt && context.collect_violating) {
    // Newly-arriving check_states: reroute violating ones to Fin_R.
    std::set<unsigned> violating_preds{};
    std::set<unsigned> violating_succs{};
    for (unsigned s : st) {
      bool is_viol = false;
      for (const auto& t : aut->out(s)) {
        if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
          if (t.acc & this->color) {
            is_viol = true;
            violating_succs.insert(t.dst);
          } else {
            succs.insert(t.dst);
          }
        }
      }
      if (is_viol) violating_preds.insert(s);
    }
    NodeContext out_ctx{};
    out_ctx.violating_predecessors = violating_preds;
    out_ctx.violating_states = violating_succs;
    return {{check_macrostate::fin(opts, std::move(succs), this->color, this->id), out_ctx}};
  }

  (void)opts;
  (void)context;  // unused in normal path
  for (unsigned s : st) {
    for (const auto& t : aut->out(s)) {
      if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
        if (t.acc & this->color) {
          return {};
        }
        succs.insert(t.dst);
      }
    }
  }
  return {{check_macrostate::fin(std::move(opts), std::move(succs), this->color, this->id), NodeContext{}}};
}

bool fin_leaf::is_satisfied() const {
  return true;
}

/**
 * Compute successor macrostate(s) for an `Inf` leaf.
 *
 * The method forms the union of this leaf's `track` set and the
 * incoming `check_states`, and collects all reachable destinations
 * (restricted to the same SCC and where the transition guard is
 * implied by `bdd`). If `check_states` is empty the function also
 * evaluates transitions from the leaf's `breakpoint` set and
 * collects `succ_break` while ignoring transitions that carry an
 * accepting mark intersecting `color`.
 *
 * - If `check_states` is empty: returns a single
 *   `check_macrostate::inf(succs, succ_break)` where `succ_break` is
 *   filtered by acceptance marks.
 * - Otherwise: returns `check_macrostate::inf(succs, succs)` (second
 *   component is a copy of `succs`).
 *
 * @param aut Pointer to the Spot automaton graph.
 * @param scc_info SCC information for the automaton (used to restrict
 *                 transitions to the current SCC).
 * @param check_states Set of automaton states currently being checked/tracked.
 * @param bdd BDD representing the current input/condition used to test
 *            transition guards.
 * @return Vector containing a single `check_macrostate::inf` with the
 *         successor sets `(succs, succ_break)` as described above.
 */
std::vector<std::pair<check_macrostate, NodeContext>> inf_leaf::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  options_ptr                       opts,
  const bdd&                        bdd,
  bool                              resample,
  NodeContext                      context
) const {

  std::set<unsigned> st = get_set_union(this->track, check_states);
  std::set<unsigned> succs = kofola::get_all_successors_in_scc(aut, scc_info, st, bdd);
  std::set<unsigned> succ_break {};

  return {{check_macrostate::inf(std::move(opts), std::move(succs), std::move(succ_break), this->color, this->id), context}};
}

std::set<unsigned> inf_leaf::get_succ_breakpoint(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const bdd&                        bdd,
  const std::set<unsigned>&         breakpoint
) const {

  std::set<unsigned> succ_break {};

  for (unsigned s : breakpoint)
    for (const auto& t : aut->out(s)) {
      if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
        if (t.acc & this->color) {
          continue;
        }
        succ_break.insert(t.dst);
      }
    }

  return succ_break;
}

bool inf_leaf::is_satisfied() const {
  return this->breakpoint.empty();
}

std::string mstate_sd_inductive::to_string() const {
  std::string res = "[SD-INDUCTIVE: ";
  res += "C=" + std::to_string(this->check_);
  res += ", Tree=" + this->check_tree_.to_string();
  res += " | ";
  res += std::to_string(this->breakpoint_);
  res += ", ";
  res += "inf_id:" + std::to_string(this->current_active_inf_id_);
  res += "]";
  return res;
}

bool mstate_sd_inductive::eq(const mstate& rhs) const {
  const auto* rhs_sd = dynamic_cast<const mstate_sd_inductive*>(&rhs);
  assert(rhs_sd);

  return (this->check_ == rhs_sd->check_) &&
         (this->check_tree_ == rhs_sd->check_tree_) &&
         (this->current_active_inf_id_ == rhs_sd->current_active_inf_id_) &&
         (this->breakpoint_ == rhs_sd->breakpoint_);
}

bool mstate_sd_inductive::lt(const mstate& rhs) const {
  const auto* rhs_sd = dynamic_cast<const mstate_sd_inductive*>(&rhs);
  assert(rhs_sd);

  if (this->check_ != rhs_sd->check_) {
    return this->check_ < rhs_sd->check_;
  }
  if (this->check_tree_ != rhs_sd->check_tree_) {
    return this->check_tree_ < rhs_sd->check_tree_;
  }
  if (this->current_active_inf_id_ != rhs_sd->current_active_inf_id_) {
    return this->current_active_inf_id_ < rhs_sd->current_active_inf_id_;
  }
  if (this->breakpoint_ != rhs_sd->breakpoint_) {
    return this->breakpoint_ < rhs_sd->breakpoint_;
  }

  return false;
}

} // namespace sd_inductive


complement_sd_inductive::complement_sd_inductive(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index) { 

  this->first_inf_leaf_id_ = 0;
  
  spot::acc_cond::acc_code acc = this->info_.part_to_acc_map_.at(part_index_).get_acceptance();
  this->acc_cond_ = acc.complement();
  
  this->is_inf_leaf_ = true;
  if(acc.fin_one() == -1) // fin is inf_leaf in complement, it is easier to obtain this info before complementation of the acc condition
    is_inf_leaf_ = false;

  this->opts_ = std::make_shared<sd_inductive::options>(sd_inductive::options{
    .use_shared_breakpoint = false,
    .use_or_fin_opt = (kofola::OPTIONS.params["sd_ind_or_opt"] == "yes")
  });
}

/**
 * Build the initial set of macrostates for this partition.
 *
 * The implementation creates a single `mstate_sd_inductive` in the
 * `GUESS` mode whose `check_` set contains the automaton's original
 * initial state only if that state belongs to this partition
 * (`part_index_`). The check-tree is constructed from the stored
 * acceptance condition `acc_cond_`.
 *
 * @return A `mstate_set` containing the initial macrostate.
 */
mstate_set complement_sd_inductive::get_init() { // {{
  DEBUG_PRINT_LN("init SD-INDUCTIVE for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state {};

  std::shared_ptr<sd_inductive::mstate_sd_inductive> derived_ms(
                                                                new sd_inductive::mstate_sd_inductive(
                                                                init_state,
                                                                sd_inductive::check_macrostate::from_acc_code(this->opts_, this->acc_cond_)));
  derived_ms->current_active_inf_id_ = 0;
  derived_ms->breakpoint_ = std::set<unsigned>{};
  if(this->is_inf_leaf_) {
    this->first_inf_leaf_id_ = derived_ms->check_tree_.find_first_inf_leaf().value();
    derived_ms->current_active_inf_id_ = this->first_inf_leaf_id_;
    auto first_leaf = derived_ms->check_tree_.find_inf_leaf_by_id(this->first_inf_leaf_id_);
    derived_ms->breakpoint_ = first_leaf.value().track;
  }

  std::shared_ptr<mstate> new_ms = derived_ms;
  mstate_set result = {new_ms};
  return result;
} // get_init() }}}


/**
 * Compute successor macrostates for an active source mstate.
 *
 * Behavior summary:
 * - `resample` is ignored (currently unused).
 * - Computes `succ_check` as all successors (within the same SCC)
 *   of the source's `check_` set under `symbol`.
 * - Computes successor check-trees from the source's `check_tree_`.
 * - If the source is a `GUESS` state: returns successors with
 *   `GUESS` for the trees computed from an empty incoming check set,
 *   and `CHECK` for trees computed with the source's `check_` set.
 * - If the source is a `CHECK` state and its check-tree is satisfied:
 *   converts the provided `glob_reached` (restricted to the partition)
 *   into a `GUESS` successor and attaches color set `{0}`.
 * - Otherwise returns `CHECK` successors built from `succ_check` and
 *   the successor trees.
 *
 * @param glob_reached set of globally reached states (used when
 *        promoting satisfied CHECK trees to GUESS with full SCC reach).
 * @param src pointer to the source `mstate` (expected to be
 *        `mstate_sd_inductive`).
 * @param symbol BDD describing the input/transition condition.
 * @param resample boolean flag (ignored in this implementation).
 * @return a collection of pairs `(mstate, colors)` describing successor
 *         macrostates and their associated color sets.
 */
mstate_col_set complement_sd_inductive::get_succ_active(
  const std::set<unsigned>& glob_reached,
  const mstate* src,
  const bdd& symbol,
  bool resample) {
  
  (void)resample; // it should be true as shared breakpoint is not used
  DEBUG_PRINT_LN("computing successor for glob_reached = " + std::to_string(glob_reached) +
    ", " + std::to_string(*src) + " over " + std::to_string(symbol));

  const sd_inductive::mstate_sd_inductive* src_mst = dynamic_cast<const sd_inductive::mstate_sd_inductive*>(src);
  assert(src_mst);

  mstate_col_set result {};
  std::set<unsigned> empty{};
  sd_inductive::NodeContext context{};

  std::set<unsigned> succ_check = kofola::get_all_successors_in_scc(this->info_.aut_, this->info_.scc_info_, src_mst->check_, symbol);
  auto succ_trees = src_mst->check_tree_.get_succ(this->info_.aut_, this->info_.scc_info_, empty, symbol, false, context);

  kofola::sd_inductive::inf_leaf inf_leaf;
  std::set<unsigned> succ_breakpoint = {};
  if(this->is_inf_leaf_) {
    inf_leaf = src_mst->check_tree_.find_inf_leaf_by_id(src_mst->current_active_inf_id_).value();
    succ_breakpoint = inf_leaf.get_succ_breakpoint(this->info_.aut_, this->info_.scc_info_, symbol, src_mst->breakpoint_);
  }
  bool forward_br = succ_breakpoint.empty();

  // For the FALSE acceptance condition we generate accepting mark only if we reachable set of states is empty. 
  // Because for non-complete automata the FALSE condition satisfies runs that are not in the automaton structure.
  bool completed_rr_cycle = !this->is_inf_leaf_  || (forward_br && src_mst->check_tree_.find_next_inf_leaf(src_mst->current_active_inf_id_).value() == this->first_inf_leaf_id_);
  
  if(src_mst->check_.empty() && completed_rr_cycle) {

    std::set<unsigned> colors = {0};
    std::set<unsigned> full_scc_reach = {};
    for (unsigned s : glob_reached) {
      if (this->info_.st_to_part_map_.at(s) == static_cast<int>(this->part_index_)) {
        full_scc_reach.insert(s);
      }
    }

    if(!this->acc_cond_.is_f() || full_scc_reach.empty()) {
      for(const auto& tree : succ_trees) {
        // OR-FIN opt: discard results with unhandled violating predecessor states
        if (!tree.second.violating_predecessors.empty()) continue;
        std::shared_ptr<sd_inductive::mstate_sd_inductive> derived_ms(
          new sd_inductive::mstate_sd_inductive(full_scc_reach, tree.first));

        if (this->is_inf_leaf_) {
          derived_ms->current_active_inf_id_ = this->first_inf_leaf_id_;
          auto first_leaf = derived_ms->check_tree_.find_inf_leaf_by_id(this->first_inf_leaf_id_);
          derived_ms->breakpoint_ = first_leaf.value().track;
        } else {
          derived_ms->current_active_inf_id_ = 0;
          derived_ms->breakpoint_ = std::set<unsigned>{};
        }

        std::shared_ptr<mstate> new_ms = derived_ms;
        result.push_back({new_ms, colors});
      }

      return result;
    }
  }

  if(!src_mst->check_.empty()) {
    // Pre-filter check_ to only states with SCC-internal successors for this symbol.
    // States with no SCC-internal successor contribute nothing to any leaf's
    // successor set and don't trigger any acceptance marks, so they can be
    // safely removed before the expensive Or-node partitioning.
    std::set<unsigned> relevant_check;
    for (unsigned s : src_mst->check_) {
      for (const auto& t : this->info_.aut_->out(s)) {
        if (this->info_.scc_info_.scc_of(s) == this->info_.scc_info_.scc_of(t.dst)
            && bdd_implies(symbol, t.cond)) {
          relevant_check.insert(s);
          break;
        }
      }
    }
    std::vector<std::pair<sd_inductive::check_macrostate, sd_inductive::NodeContext>> succ_check_trees = src_mst->check_tree_.get_succ(this->info_.aut_,
      this->info_.scc_info_, relevant_check, symbol, true, context);
    
    for(const auto& tree : succ_trees) {
      // OR-FIN opt: discard results with unhandled violating predecessor states
      if (!tree.second.violating_predecessors.empty()) continue;

      std::shared_ptr<sd_inductive::mstate_sd_inductive> derived_ms(
              new sd_inductive::mstate_sd_inductive(succ_check, tree.first));

      // check phase does not matter
      derived_ms->current_active_inf_id_ = this->first_inf_leaf_id_;
      derived_ms->breakpoint_ = std::set<unsigned>{};
        
      std::shared_ptr<mstate> new_ms = derived_ms;

      result.push_back({new_ms, {}});
    }

    for(const auto& tree : succ_check_trees) {
      // OR-FIN opt: discard results with unhandled violating predecessor states
      if (!tree.second.violating_predecessors.empty()) continue;
      
      std::shared_ptr<sd_inductive::mstate_sd_inductive> derived_ms(
              new sd_inductive::mstate_sd_inductive(empty, tree.first));

      std::set<unsigned> colors = {};
      if (this->is_inf_leaf_) {
        derived_ms->current_active_inf_id_ = this->first_inf_leaf_id_;
        auto first_leaf = derived_ms->check_tree_.find_inf_leaf_by_id(this->first_inf_leaf_id_);
        derived_ms->breakpoint_ = first_leaf.value().track;
      } else {
        derived_ms->current_active_inf_id_ = 0;
        derived_ms->breakpoint_ = std::set<unsigned>{};
        colors.insert(0);
      }
        
      std::shared_ptr<mstate> new_ms = derived_ms;

      result.push_back({new_ms, colors});
    }
    return result;
  }

  for(const auto& tree : succ_trees) {
    // OR-FIN opt: discard results with unhandled violating predecessor states
    if (!tree.second.violating_predecessors.empty()) continue;
    
    std::shared_ptr<sd_inductive::mstate_sd_inductive> derived_ms(
              new sd_inductive::mstate_sd_inductive(empty, tree.first));

    std::set<unsigned> colors = {};
    if(forward_br && this->is_inf_leaf_) { 
      auto next_id = derived_ms->check_tree_.find_next_inf_leaf(src_mst->current_active_inf_id_).value();
      auto next_leaf = derived_ms->check_tree_.find_inf_leaf_by_id(next_id);
      derived_ms->current_active_inf_id_ = next_id;
      derived_ms->breakpoint_ = next_leaf.value().track;
    } else {
      derived_ms->current_active_inf_id_ = src_mst->current_active_inf_id_;
      derived_ms->breakpoint_ = succ_breakpoint;
    }
    
    std::shared_ptr<mstate> new_ms = derived_ms;
    result.push_back({new_ms, {}});
  }
  return result;
}

} // namespace kofola
