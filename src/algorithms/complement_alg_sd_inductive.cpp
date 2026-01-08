// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_sd_inductive.hpp"

#include <cassert>
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
          return check_macrostate::fin(leaf.safe, leaf.color, id);
        } else {
          static_assert(std::is_same_v<leaf_t, inf_leaf>, "Unexpected leaf type");
          return check_macrostate::inf(leaf.track, leaf.breakpoint, leaf.color, id);
        }
      },
      tree.leaf_value());
  }

  const check_macrostate left_ms(base_tree(tree.left()));
  const check_macrostate right_ms(base_tree(tree.right()));
  const auto node_type = tree.type();
  return check_macrostate::make(node_type, assign_leaf_ids(left_ms, next_id), assign_leaf_ids(right_ms, next_id));
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
check_macrostate check_macrostate::from_acc_code(const spot::acc_cond::acc_code& code) {
  const auto parsed = from_acc_code_impl(code);
  unsigned next_id = 0;
  return assign_leaf_ids(parsed, next_id);
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
std::vector<check_macrostate> check_macrostate::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  const bdd&                        bdd,
  bool                              resample,
  const NodeContext&                parent_context) const {

  if (this->is_leaf()) {
    return std::visit(
      [&](const auto& leaf) {
        return leaf.get_succ(aut, scc_info, check_states, bdd, resample, parent_context);
      },
      this->leaf_value());
  }

  // Propagate/merge context from parent into this node.
  NodeContext context = parent_context;
  if (this->type() == TreeType::And || this->type() == TreeType::Or) {
    const NodeContext local = this->node_value().get_context();
    context = local.merge_contexts(parent_context);
  }

  const auto node_type = this->type();
  const check_macrostate left_ms(base_tree(this->left()));
  const check_macrostate right_ms(base_tree(this->right()));

  if (node_type == TreeType::And) {
    const auto left_succ = left_ms.get_succ(aut, scc_info, check_states, bdd, resample, context);
    const auto right_succ = right_ms.get_succ(aut, scc_info, check_states, bdd, resample, context);
    return cartesian_product<check_macrostate, check_macrostate>(
      left_succ,
      right_succ,
      [](const check_macrostate& l, const check_macrostate& r) {
        return check_macrostate::make(TreeType::And, l, r);
      });
  }

  if (node_type == TreeType::Or) {
    std::set<check_macrostate> out;
    const auto partitions = nondet_split_set(check_states, 2);
    for (const auto& part : partitions) {
      if (part.size() != 2) {
        throw std::logic_error("check_macrostate::get_succ: nondet_split_set returned non-binary partition");
      }
      const auto& left_states = part[0];
      const auto& right_states = part[1];
      const auto left_succ = left_ms.get_succ(aut, scc_info, left_states, bdd, resample, context);
      const auto right_succ = right_ms.get_succ(aut, scc_info, right_states, bdd, resample, context);

      // TODO: it is not efficient to call reduce here
      const auto combined = cartesian_product<check_macrostate, check_macrostate>(
        left_succ,
        right_succ,
        [](const check_macrostate& l, const check_macrostate& r) {
          return check_macrostate::make(TreeType::Or, l, r).reduce();
        });
      out.insert(combined.begin(), combined.end());
    }
    return std::vector<check_macrostate>(out.begin(), out.end());
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
  const check_macrostate left_ms(base_tree(this->left()));
  const check_macrostate right_ms(base_tree(this->right()));

  if (node_type == TreeType::And || node_type == TreeType::Or) {
    return left_ms.is_satisfied() && right_ms.is_satisfied();
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

  const check_macrostate left_ms(base_tree(this->left()));
  const check_macrostate right_ms(base_tree(this->right()));
  return get_set_union(left_ms.gather_states(), right_ms.gather_states());
}

namespace {

/**
 * Return a copy of `tree` with all states from `forbidden` removed
 * from any leaf state-sets.
 *
 * - For a `fin_leaf`, the `safe` set is replaced with
 *   `get_set_difference(leaf.safe, forbidden)`.
 * - For an `inf_leaf`, the `track` and `breakpoint` sets are
 *   replaced with `get_set_difference(..., forbidden)`.
 *
 * Internal nodes (`And` / `Or`) are processed recursively and rebuilt
 * with the restricted children.
 *
 * @param tree The input `check_macrostate` to restrict.
 * @param forbidden Set of automaton state indices to remove from
 *                  leaf sets.
 * @return A new `check_macrostate` instance equivalent to `tree` but
 *         with `forbidden` removed from all leaf state-sets.
 */
sd_inductive::check_macrostate restrict_states_in_tree(
  const sd_inductive::check_macrostate& tree,
  const std::set<unsigned>& forbidden) {

  using check_macrostate = sd_inductive::check_macrostate;
  using base_tree = kofola::types::binary_tree<sd_inductive::TreeType, sd_inductive::AndOrNode, sd_inductive::fin_leaf, sd_inductive::inf_leaf>;

  if (tree.is_leaf()) {
    return std::visit(
      [&](const auto& leaf) -> check_macrostate {
        using leaf_t = std::decay_t<decltype(leaf)>;
        if constexpr (std::is_same_v<leaf_t, sd_inductive::fin_leaf>) {
          return check_macrostate::fin(get_set_difference(leaf.safe, forbidden), leaf.color, leaf.id);
        } else {
          static_assert(std::is_same_v<leaf_t, sd_inductive::inf_leaf>, "Unexpected leaf type");
          return check_macrostate::inf(
            get_set_difference(leaf.track, forbidden),
            get_set_difference(leaf.breakpoint, forbidden),
            leaf.color,
            leaf.id);
        }
      },
      tree.leaf_value());
  }

  const auto node_type = tree.type();
  check_macrostate left(base_tree(tree.left()));
  check_macrostate right(base_tree(tree.right()));
  return check_macrostate::make(
    node_type,
    restrict_states_in_tree(left, forbidden),
    restrict_states_in_tree(right, forbidden));
}

} // namespace

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
  const check_macrostate left_ms(base_tree(this->left()));
  const check_macrostate right_ms(base_tree(this->right()));

  if (node_type == TreeType::And) {
    return check_macrostate::make(TreeType::And, left_ms.reduce(), right_ms.reduce());
  }

  if (node_type == TreeType::Or) {
    const check_macrostate left_red = left_ms.reduce();
    const auto left_states = left_red.gather_states();
    const check_macrostate right_restricted = restrict_states_in_tree(right_ms, left_states);
    return check_macrostate::make(TreeType::Or, left_red, right_restricted.reduce());
  }

  throw std::logic_error("check_macrostate::reduce: unexpected internal node type");
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
check_macrostate check_macrostate::fold(TreeType op, const std::vector<spot::acc_cond::acc_code>& parts) {
  if (parts.empty()) {
    throw std::invalid_argument("check_macrostate: empty And/Or in acceptance formula");
  }
  check_macrostate acc = from_acc_code_impl(parts.front());
  for (size_t i = 1; i < parts.size(); ++i) {
    acc = check_macrostate::make(op, std::move(acc), from_acc_code_impl(parts[i]));
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
check_macrostate check_macrostate::from_acc_code_impl(const spot::acc_cond::acc_code& code) {
  if (code.empty()) {
    throw std::invalid_argument("check_macrostate: empty acceptance formula");
  }

  // Leaf: [mark][op]
  if (code.size() == 2 && false) {
    const auto op = code[1].sub.op;
    const auto mark = code[0].mark;
    if (op == spot::acc_cond::acc_op::Fin) {
      return check_macrostate(base_tree::leaf(TreeType::Fin, fin_leaf{{}, mark, 0}));
    }
    if (op == spot::acc_cond::acc_op::Inf) {
      return check_macrostate(base_tree::leaf(TreeType::Inf, inf_leaf{{}, {}, mark, 0}));
    }
  }

  // Prefer top-level flattening (Spot returns a singleton vector when the operator is not present at top-level).
  const auto conjuncts = code.top_conjuncts();
  if (conjuncts.size() > 1) {
    return fold(TreeType::And, conjuncts);
  }
  const auto disjuncts = code.top_disjuncts();
  if (disjuncts.size() > 1) {
    return fold(TreeType::Or, disjuncts);
  }

  if (code.size() == 2) {
    const auto op = code[1].sub.op;
    const auto mark = code[0].mark;
    if (op == spot::acc_cond::acc_op::Fin) {
      return check_macrostate(base_tree::leaf(TreeType::Fin, fin_leaf{{}, mark, 0}));
    }
    if (op == spot::acc_cond::acc_op::Inf) {
      return check_macrostate(base_tree::leaf(TreeType::Inf, inf_leaf{{}, {}, mark, 0}));
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

std::string check_macrostate::to_string_impl(const base_tree& tree) {
  const std::string head = tree_type_to_string(tree.type());
  if (tree.is_leaf()) {
    const std::string payload = std::visit(
      [](const auto& leaf) { return leaf.to_string(); },
      tree.leaf_value());
    return head + "(" + payload + ")";
  }

  return head + "(" + to_string_impl(tree.left()) + ", " + to_string_impl(tree.right()) + ")";
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
std::vector<check_macrostate> fin_leaf::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  const bdd&                        bdd,
  bool                              resample,
  const NodeContext&                context) const {

  (void)resample; // unused
  (void)context;  // unused
  std::set<unsigned> st = get_set_union(this->safe, check_states);
  std::set<unsigned> succs {};
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
  return {check_macrostate::fin(std::move(succs), this->color, this->id)};
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
std::vector<check_macrostate> inf_leaf::get_succ(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         check_states,
  const bdd&                        bdd,
  bool                              resample,
  const NodeContext&                context) const {

  std::set<unsigned> st = get_set_union(this->track, check_states);
  std::set<unsigned> succs {};
  std::set<unsigned> succ_break {};
  for (unsigned s : st) {
    for (const auto& t : aut->out(s)) {
      if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
        succs.insert(t.dst);
      }
    }
  }

  if (!resample) {
    const std::set<unsigned>* breakpoint_src = &this->breakpoint;
    if (context.type == NodeContextType::SHARED_BREAKPOINT &&
        context.breakpoint.has_value() &&
        context.leaf_id == this->id) {
      breakpoint_src = &context.breakpoint->get();
    }

    for (unsigned s : *breakpoint_src) {
      for (const auto& t : aut->out(s)) {
        if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
          if (t.acc & this->color) {
            continue;
          }
          succ_break.insert(t.dst);
        }
      }
    }

    // If using a shared breakpoint context, update it in-place.
    if (context.type == NodeContextType::SHARED_BREAKPOINT &&
        context.breakpoint.has_value() &&
        context.leaf_id == this->id) {
      context.breakpoint->get() = succ_break;
    }
    return {check_macrostate::inf(std::move(succs), std::move(succ_break), this->color, this->id)};
  }

  auto succs_copy = succs;
  return {check_macrostate::inf(std::move(succs), std::move(succs_copy), this->color, this->id)};
}

bool inf_leaf::is_satisfied() const {
  return this->breakpoint.empty();
}

namespace {

const char* mstate_type_to_string(mstate_type t) {
  switch (t) {
    case mstate_type::GUESS:
      return "GUESS";
    case mstate_type::CHECK:
      return "CHECK";
  }
  return "?";
}

} // namespace

std::string mstate_sd_inductive::to_string() const {
  std::string res = "[SD-INDUCTIVE: ";
  res += std::string("Type=") + mstate_type_to_string(this->type_);
  res += ", C=" + std::to_string(this->check_);
  res += ", Tree=" + this->check_tree_.to_string();
  res += "]";
  return res;
}

bool mstate_sd_inductive::eq(const mstate& rhs) const {
  const auto* rhs_sd = dynamic_cast<const mstate_sd_inductive*>(&rhs);
  assert(rhs_sd);
  return (this->type_ == rhs_sd->type_) &&
         (this->check_ == rhs_sd->check_) &&
         (this->check_tree_ == rhs_sd->check_tree_);
}

bool mstate_sd_inductive::lt(const mstate& rhs) const {
  const auto* rhs_sd = dynamic_cast<const mstate_sd_inductive*>(&rhs);
  assert(rhs_sd);

  if (this->type_ != rhs_sd->type_) {
    return this->type_ < rhs_sd->type_;
  }
  if (this->check_ != rhs_sd->check_) {
    return this->check_ < rhs_sd->check_;
  }
  if (this->check_tree_ != rhs_sd->check_tree_) {
    return this->check_tree_ < rhs_sd->check_tree_;
  }

  return false;
}

} // namespace sd_inductive


complement_sd_inductive::complement_sd_inductive(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index) { 
  
  spot::acc_cond::acc_code acc = this->info_.part_to_acc_map_.at(part_index_).get_acceptance();
  this->acc_cond_ = acc.complement();
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
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new sd_inductive::mstate_sd_inductive(init_state, sd_inductive::check_macrostate::from_acc_code(this->acc_cond_), sd_inductive::mstate_type::GUESS));
  mstate_set result = {ms};
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

  std::set<unsigned> succ_check = kofola::get_all_successors_in_scc(
      this->info_.aut_, this->info_.scc_info_, src_mst->check_, symbol);
  std::vector<sd_inductive::check_macrostate> succ_trees = src_mst->check_tree_.get_succ(this->info_.aut_, 
      this->info_.scc_info_, empty, symbol, false);
  if(src_mst->type_ == sd_inductive::mstate_type::GUESS) {
    std::vector<sd_inductive::check_macrostate> succ_check_trees = src_mst->check_tree_.get_succ(this->info_.aut_, 
      this->info_.scc_info_, src_mst->check_, symbol, true);
    for(const auto& tree : succ_trees) {
      std::shared_ptr<mstate> new_ms(new sd_inductive::mstate_sd_inductive(
          succ_check, tree, sd_inductive::mstate_type::GUESS));
      result.push_back({new_ms, {}});
    }
    for(const auto& tree : succ_check_trees) {
      std::shared_ptr<mstate> new_ms(new sd_inductive::mstate_sd_inductive(
          succ_check, tree, sd_inductive::mstate_type::CHECK));
      result.push_back({new_ms, {}});
    }
  } else if(src_mst->check_tree_.is_satisfied()) {
    std::set<unsigned> colors = {0};
    std::set<unsigned> full_scc_reach = {};
    for (unsigned s : glob_reached) {
      if (this->info_.st_to_part_map_.at(s) == static_cast<int>(this->part_index_)) {
        full_scc_reach.insert(s);
      }
    }
    for(const auto& tree : succ_trees) {
      std::shared_ptr<mstate> new_ms(new sd_inductive::mstate_sd_inductive(
          full_scc_reach, tree, sd_inductive::mstate_type::GUESS));
      result.push_back({new_ms, colors});
    }
  } else {
    for(const auto& tree : succ_trees) {
      std::shared_ptr<mstate> new_ms(new sd_inductive::mstate_sd_inductive(
          succ_check, tree, sd_inductive::mstate_type::CHECK));
      result.push_back({new_ms, {}});
    }
  }

  return result;
}

} // namespace kofola
