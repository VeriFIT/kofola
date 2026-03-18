// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

#include <compare>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../types/binary_tree.hpp"
#include "../util/sets.hpp"

// SPOT
#include <spot/twa/acc.hh>

namespace kofola { // {{{

namespace sd_inductive {

  struct options {
    bool use_shared_breakpoint{false};
    bool use_or_fin_opt{false};
  };

  using options_ptr = std::shared_ptr<const options>;

  /**
   * @brief Type of a node in a `check_macrostate` tree.
   *
   * Leaves are typed as `Fin` and `Inf`. Internal nodes are typed as logical
   * connectives `And` and `Or`.
   */
  enum class TreeType {
    Fin,
    Inf,
    And,
    Or
  };

  enum class NodeContextType {
    NONE,
    SHARED_BREAKPOINT
  };

  enum class NodeContextState {
    GLOBAL_WAIT,
    RESAMPLE_LEAF,
    PROCESS_LEAF,
  };

  class check_macrostate;

  /**
   * @brief Collect IDs of `Inf` leaves under `And` nodes.
   *
   * Traverses the check tree @p t and appends each encountered
   * `inf_leaf::id` to @p out. For internal nodes, only `TreeType::And`
   * subtrees are traversed; `Or` subtrees are intentionally ignored.
   *
   * @param t   Check tree (subtree root) to traverse.
   * @param out Output vector to append leaf IDs into.
   */
  inline void collect_inf_leaf_ids(const check_macrostate& t, std::vector<unsigned>& out);

  struct NodeContext;
  std::ostream& operator<<(std::ostream& os, const NodeContext& ctx);

  struct NodeContext {
    NodeContextType type { NodeContextType::NONE };

    std::set<unsigned> breakpoint {};
    unsigned leaf_id {0};
    NodeContextState state {NodeContextState::GLOBAL_WAIT};

    // OR-FIN optimization: transient fields (not part of macrostate identity)
    bool collect_violating {false};                  // downward: tells FIN leaves to collect violations
    std::set<unsigned> violating_states {};          // upward: successors of violating states (for restrict_states_in_tree)
    std::set<unsigned> violating_predecessors {};    // upward: source states that fired Fin transitions (for right subtree check_states)

    // ----------------------------------------------------------------
    // Predicates — prefer these over direct field comparisons in callers
    // ----------------------------------------------------------------

    /// True when this context carries no special optimization (default state).
    bool is_none() const { return type == NodeContextType::NONE; }

    /// True when this context runs the shared-breakpoint optimization.
    bool is_shared_breakpoint() const { return type == NodeContextType::SHARED_BREAKPOINT; }

    /// True when this context is currently targeting the leaf with the given @p id.
    bool targets_leaf(unsigned id) const { return leaf_id == id; }

    /**
     * @brief True when the shared-breakpoint cycle for this And-node is complete.
     *
     * An And-node is cycle-complete when either the context is not a
     * shared-breakpoint context at all, or the breakpoint is empty and the
     * state machine has returned to `GLOBAL_WAIT`.
     */
    bool is_cycle_completed() const {
      return !is_shared_breakpoint() ||
             (breakpoint.empty() && state == NodeContextState::GLOBAL_WAIT);
    }

    // ----------------------------------------------------------------
    // State-machine
    // ----------------------------------------------------------------

    /**
     * @brief Compute the successor context for one transition step.
     *
     * Implements the small state machine used by the shared-breakpoint
     * optimization:
     * - If not a `SHARED_BREAKPOINT` context, returns itself unchanged.
     * - In `GLOBAL_WAIT`: a `resample` request advances to `RESAMPLE_LEAF`.
     * - In `PROCESS_LEAF`: advances the leaf index when the breakpoint empties;
     *   wrapping back to index 0 returns to `GLOBAL_WAIT`.
     *
     * @param resample Whether the caller requests a resampling step.
     * @return A copy of this context updated for the next step.
     */
    NodeContext get_succ_context(bool resample) const {
      NodeContext succ = *this;
      if (!this->is_shared_breakpoint()) {
        return succ;
      }

      if (this->state == NodeContextState::GLOBAL_WAIT) {
        if (resample) {
          succ.state = NodeContextState::RESAMPLE_LEAF;
        }
        // else stay in GLOBAL_WAIT
      } else if (this->state == NodeContextState::RESAMPLE_LEAF) {
        assert(false);
        return succ;
      } else {
        // PROCESS_LEAF: advance the cyclic leaf pointer when breakpoint drains
        if (succ.leaf_ids_.size() > 0) {
          succ.leaf_index_ = succ.leaf_index_ % succ.leaf_ids_.size();
          if (succ.breakpoint.empty()) {
            succ.leaf_index_ = (succ.leaf_index_ + 1) % succ.leaf_ids_.size();
            if (succ.leaf_index_ == 0) {
              succ.state = NodeContextState::GLOBAL_WAIT;
            } else {
              succ.state = NodeContextState::RESAMPLE_LEAF;
            }
          }
        }
        succ.leaf_id = succ.leaf_ids_[succ.leaf_index_];
      }

      return succ;
    }

    // ----------------------------------------------------------------
    // Factory
    // ----------------------------------------------------------------

    /**
     * @brief Build a shared-breakpoint context for a subtree.
     *
     * Collects all `Inf` leaf IDs from @p subtree_ and, if the subtree root
     * is an `And` node with at least one `Inf` leaf, initializes the returned
     * context as `SHARED_BREAKPOINT` and selects the first leaf ID.
     *
     * @param t        Type of the subtree root.
     * @param subtree_ Subtree to inspect for `Inf` leaf IDs.
     * @return A freshly initialized context for that subtree.
     */
    static NodeContext create_subtree_sh_context(TreeType t, const check_macrostate& subtree_) {
      NodeContext ctx;
      collect_inf_leaf_ids(subtree_, ctx.leaf_ids_);
      if (t == TreeType::And && ctx.leaf_ids_.size() > 0) {
        ctx.type = NodeContextType::SHARED_BREAKPOINT;
        ctx.leaf_id = ctx.leaf_ids_[ctx.leaf_index_ % ctx.leaf_ids_.size()];
      }
      return ctx;
    }

    // ----------------------------------------------------------------
    // Structural
    // ----------------------------------------------------------------

    /**
     * @brief True when this context opens a new shared-breakpoint scope.
     *
     * A context is a scope root when it carries a shared-breakpoint type
     * while its parent does not, i.e. this is the topmost And-node that
     * introduces the optimization.
     *
     * @param parent Context from the enclosing (parent) node.
     */
    bool is_scope_root(const NodeContext& parent) const {
      return is_shared_breakpoint() && !parent.is_shared_breakpoint();
    }

    /**
     * @brief Combine this context with a sibling context coming from the
     *        other child of an And/Or node.
     *
     * `NONE` is the identity element: if either side is `NONE`, the other
     * is returned.  When both sides carry a context, `*this` wins.
     *
     * @param other Context returned by the sibling subtree.
     * @return The combined context to propagate upward.
     */
    NodeContext union_contexts(const NodeContext& other) const {
      NodeContext result;
      if (this->is_none()) {
        result = other;
      } else {
        result = *this;
      }
      // Always union transient violation fields from both sides
      result.violating_states = get_set_union(this->violating_states, other.violating_states);
      result.violating_predecessors = get_set_union(this->violating_predecessors, other.violating_predecessors);
      return result;
    }

    void restrict_states(const std::set<unsigned>& forbidden) {
      this->breakpoint = get_set_difference(this->breakpoint, forbidden);
    }

    // ----------------------------------------------------------------
    // Comparison & serialization
    // ----------------------------------------------------------------

    /// Equality comparison (excludes transient fields: collect_violating, violating_states).
    bool operator==(const NodeContext& other) const {
      return type == other.type &&
             breakpoint == other.breakpoint &&
             leaf_id == other.leaf_id &&
             leaf_index_ == other.leaf_index_ &&
             leaf_ids_ == other.leaf_ids_ &&
             state == other.state;
    }

    /// Three-way comparison (excludes transient fields: collect_violating, violating_states).
    std::strong_ordering operator<=>(const NodeContext& other) const {
      if (this->type < other.type) return std::strong_ordering::less;
      if (other.type < this->type) return std::strong_ordering::greater;

      if (this->breakpoint < other.breakpoint) return std::strong_ordering::less;
      if (other.breakpoint < this->breakpoint) return std::strong_ordering::greater;

      if (this->leaf_id < other.leaf_id) return std::strong_ordering::less;
      if (other.leaf_id < this->leaf_id) return std::strong_ordering::greater;

      if (this->leaf_index_ < other.leaf_index_) return std::strong_ordering::less;
      if (other.leaf_index_ < this->leaf_index_) return std::strong_ordering::greater;

      if (this->leaf_ids_ < other.leaf_ids_) return std::strong_ordering::less;
      if (other.leaf_ids_ < this->leaf_ids_) return std::strong_ordering::greater;

      if (this->state < other.state) return std::strong_ordering::less;
      if (other.state < this->state) return std::strong_ordering::greater;

      return std::strong_ordering::equal;
    }

    std::string to_string() const {
      std::ostringstream os;
      os << *this;
      return os.str();
    }

  private:
    /// Index of the currently active leaf inside `leaf_ids_`.
    unsigned leaf_index_ {0};
    /// Ordered list of `Inf` leaf IDs under the enclosing `And` node.
    std::vector<unsigned> leaf_ids_ {};
  };

  inline std::ostream& operator<<(std::ostream& os, const NodeContext& ctx) {
    os << "NodeContext{";
    switch (ctx.type) {
      case NodeContextType::NONE:
        os << "type=NONE";
        break;
      case NodeContextType::SHARED_BREAKPOINT:
        os << "type=SHARED_BREAKPOINT";
        break;
    }
    os << ", leaf_id=" << ctx.leaf_id;
    os << ", breakpoint=" << std::to_string(ctx.breakpoint);
    os << ", state=" << (int)ctx.state;
    os << "}";
    return os;
  }

  /**
   * @brief Payload of an internal `And`/`Or` node in `check_macrostate`.
   */
  struct AndOrNode {
    TreeType type {TreeType::And};

    NodeContext context{};

    AndOrNode() = default;
    AndOrNode(TreeType t, NodeContext context_)
      : type(t), context(std::move(context_)) {}

    bool operator==(const AndOrNode& other) const;
    std::strong_ordering operator<=>(const AndOrNode& other) const;

    /// True when the shared-breakpoint cycle for this node has completed
    /// (or the context is not a shared-breakpoint context at all).
    bool is_satisfied() const {
      return this->context.is_cycle_completed();
    }

    NodeContext get_succ_context(bool resample) const {
      return this->context.get_succ_context(resample);
    }

    const NodeContext& get_context() const {
      return this->context;
    }

    void set_context(const NodeContext& ctx) {
      this->context = ctx;
    } 
  };

  /**
   * @brief Payload of a `Fin` leaf in `check_macrostate`.
   */
  struct fin_leaf {

    /** @brief Set of safe states relevant for the Fin-check. */
    std::set<unsigned> safe {};
    /** @brief Acceptance color associated with this Fin-check. */
    spot::acc_cond::mark_t color {}; 

    /** @brief Unique numeric id of this leaf within a check tree. */
    unsigned id {0};

    static std::string mark_to_string(const spot::acc_cond::mark_t& m) {
      std::ostringstream os;
      os << m;
      return os.str();
    }

    /**
     * @brief Convert this leaf payload into a human-readable string.
     *
     * @return String representation of this Fin leaf.
     */
    std::string to_string() const {
      return "safe=" + std::to_string(this->safe) + ", color=" + mark_to_string(this->color);
    }

    /**
     * @brief Defaulted three-way comparison.
     *
     * @param other The other leaf payload to compare with.
     * @return Ordering relation result.
     */
    bool operator==(const fin_leaf&) const = default;

    std::strong_ordering operator<=>(const fin_leaf& other) const {
      if (this->id < other.id)
        return std::strong_ordering::less;
      if (other.id < this->id)
        return std::strong_ordering::greater;
      if (this->safe < other.safe)
        return std::strong_ordering::less;
      if (other.safe < this->safe)
        return std::strong_ordering::greater;
      if (this->color < other.color)
        return std::strong_ordering::less;
      if (other.color < this->color)
        return std::strong_ordering::greater;
      return std::strong_ordering::equal;
    }

    std::vector<std::pair<check_macrostate, NodeContext>> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      options_ptr                        opts,
      const bdd&                        bdd,
      bool                              resample,
      NodeContext                       context) const;

    bool is_satisfied() const;
  };

  /**
   * @brief Payload of an `Inf` leaf in `check_macrostate`.
   */
  struct inf_leaf {
    /** @brief Set of tracked states relevant for the Inf-check. */
    std::set<unsigned> track {};
    /** @brief Breakpoint set used in the Inf-check. */
    std::set<unsigned> breakpoint {};
    /** @brief Acceptance color associated with this Inf-check. */
    spot::acc_cond::mark_t color {}; 

    /** @brief Unique numeric id of this leaf within a check tree. */
    unsigned id {0};

    /**
     * @brief Convert this leaf payload into a human-readable string.
     *
     * @return String representation of this Inf leaf.
     */
    std::string to_string() const {
      return "track=" + std::to_string(this->track) + ", breakpoint=" + std::to_string(this->breakpoint)
             + ", color=" + fin_leaf::mark_to_string(this->color);
    }

    /**
     * @brief Defaulted three-way comparison.
     *
     * @param other The other leaf payload to compare with.
     * @return Ordering relation result.
     */
    bool operator==(const inf_leaf&) const = default;

    std::strong_ordering operator<=>(const inf_leaf& other) const {
      if (this->id < other.id)
        return std::strong_ordering::less;
      if (other.id < this->id)
        return std::strong_ordering::greater;
      if (this->track < other.track)
        return std::strong_ordering::less;
      if (other.track < this->track)
        return std::strong_ordering::greater;
      if (this->breakpoint < other.breakpoint)
        return std::strong_ordering::less;
      if (other.breakpoint < this->breakpoint)
        return std::strong_ordering::greater;
      if (this->color < other.color)
        return std::strong_ordering::less;
      if (other.color < this->color)
        return std::strong_ordering::greater;
      return std::strong_ordering::equal;
    }

    std::vector<std::pair<check_macrostate, NodeContext>> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      options_ptr                        opts,
      const bdd&                        bdd,
      bool                              resample,
      NodeContext                       context) const;

    bool is_satisfied() const;

  };

  /**
   * @brief Stream output for `fin_leaf`.
   *
   * @param os Output stream.
   * @param l  Leaf payload.
   * @return Reference to @p os.
   */
  inline std::ostream& operator<<(std::ostream& os, const fin_leaf& l) {
    os << l.to_string();
    return os;
  }

  /**
   * @brief Stream output for `inf_leaf`.
   *
   * @param os Output stream.
   * @param l  Leaf payload.
   * @return Reference to @p os.
   */
  inline std::ostream& operator<<(std::ostream& os, const inf_leaf& l) {
    os << l.to_string();
    return os;
  }

  /**
   * @brief A binary tree representing a macrostate check formula.
   *
   * This is a thin wrapper over `kofola::types::binary_tree` with:
   * - node type: `TreeType`
   * - internal-node payload: `AndOrNode`
   * - leaf payload: `fin_leaf` or `inf_leaf`
   */
  class check_macrostate : public kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf> {
    /** @brief Base tree type used for representation. */
    using base_tree = kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf>;

  public:
    /** @brief Deleted default constructor (a check must be a leaf or a node). */
    check_macrostate() = delete;

    /**
     * @brief Copy constructor.
     *
     * @param other Object to copy.
     */
    check_macrostate(const check_macrostate&) = default;

    /**
     * @brief Move constructor.
     *
     * @param other Object to move from.
     */
    check_macrostate(check_macrostate&&) noexcept = default;

    /**
     * @brief Copy assignment.
     *
     * @param other Object to copy.
     * @return Reference to this object.
     */
    check_macrostate& operator=(const check_macrostate&) = default;

    /**
     * @brief Move assignment.
     *
     * @param other Object to move from.
     * @return Reference to this object.
     */
    check_macrostate& operator=(check_macrostate&&) noexcept = default;

    /**
     * @brief Construct a `check_macrostate` from an already-built base tree.
     *
     * @param opts Options shared by all nodes in this check tree.
     * @param tree Underlying tree representation.
     */
    explicit check_macrostate(options_ptr opts, base_tree tree)
      : base_tree(std::move(tree)),
        opts_(opts ? std::move(opts) : default_options()) {}

    static options_ptr default_options() {
      static const options_ptr opts = std::make_shared<options>();
      return opts;
    }

    const options& get_options() const { return *opts_; }
    const options_ptr& get_options_ptr() const { return opts_; }

    /**
     * @brief Build a `Fin` leaf.
     *
     * @param safe Set of safe states.
     * @return A `check_macrostate` leaf of type `TreeType::Fin`.
     */
    static check_macrostate fin(options_ptr opts, std::set<unsigned> safe, spot::acc_cond::mark_t color) {
      return check_macrostate(std::move(opts), base_tree::leaf(TreeType::Fin, fin_leaf{std::move(safe), color, 0}));
    }

    static check_macrostate fin(options_ptr opts, std::set<unsigned> safe, spot::acc_cond::mark_t color, unsigned id) {
      return check_macrostate(std::move(opts), base_tree::leaf(TreeType::Fin, fin_leaf{std::move(safe), color, id}));
    }

    static check_macrostate fin(std::set<unsigned> safe, spot::acc_cond::mark_t color) {
      return fin(default_options(), std::move(safe), color);
    }

    static check_macrostate fin(std::set<unsigned> safe, spot::acc_cond::mark_t color, unsigned id) {
      return fin(default_options(), std::move(safe), color, id);
    }

    /**
     * @brief Build an `Inf` leaf.
     *
     * @param track      Set of tracked states.
     * @param breakpoint Breakpoint set.
     * @return A `check_macrostate` leaf of type `TreeType::Inf`.
     */
    static check_macrostate inf(options_ptr opts, std::set<unsigned> track, std::set<unsigned> breakpoint, spot::acc_cond::mark_t color) {
      return check_macrostate(std::move(opts), base_tree::leaf(TreeType::Inf, inf_leaf{std::move(track), std::move(breakpoint), color, 0}));
    }

    static check_macrostate inf(options_ptr opts, std::set<unsigned> track, std::set<unsigned> breakpoint, spot::acc_cond::mark_t color, unsigned id) {
      return check_macrostate(std::move(opts), base_tree::leaf(TreeType::Inf, inf_leaf{std::move(track), std::move(breakpoint), color, id}));
    }

    static check_macrostate inf(std::set<unsigned> track, std::set<unsigned> breakpoint, spot::acc_cond::mark_t color) {
      return inf(default_options(), std::move(track), std::move(breakpoint), color);
    }

    static check_macrostate inf(std::set<unsigned> track, std::set<unsigned> breakpoint, spot::acc_cond::mark_t color, unsigned id) {
      return inf(default_options(), std::move(track), std::move(breakpoint), color, id);
    }

    /**
     * @brief Build an internal node with two children.
     *
     * @param type  Node type (typically `TreeType::And` or `TreeType::Or`).
     * @param left  Left subtree.
     * @param right Right subtree.
     * @return A `check_macrostate` internal node.
     */
    static check_macrostate make(options_ptr opts, TreeType type, check_macrostate left, check_macrostate right, NodeContext node_payload = NodeContext{}) {
      base_tree l(std::move(left));
      base_tree r(std::move(right));
      return check_macrostate(std::move(opts), base_tree::make_node(type, AndOrNode(type, std::move(node_payload)), std::move(l), std::move(r)));
    }

    static check_macrostate make(TreeType type, check_macrostate left, check_macrostate right, NodeContext node_payload = NodeContext{}) {
      return make(left.get_options_ptr(), type, std::move(left), std::move(right), std::move(node_payload));
    }

    /**
     * @brief Convert the whole check tree into a human-readable string.
     *
     * @return String representation of this macrostate check tree.
     */
    std::string to_string() const {
      return to_string_impl(static_cast<const base_tree&>(*this), opts_ && opts_->use_shared_breakpoint);
    }

    /**
     * @brief Stream output for `check_macrostate`.
     *
     * @param os Output stream.
     * @param ms Macrostate check tree.
     * @return Reference to @p os.
     */
    friend std::ostream& operator<<(std::ostream& os, const check_macrostate& ms) {
      os << ms.to_string();
      return os;
    }

    /// Build a `check_macrostate` tree from a Spot acceptance formula.
    static check_macrostate from_acc_code(options_ptr opts, const spot::acc_cond::acc_code& code);

    static check_macrostate from_acc_code(const spot::acc_cond::acc_code& code) {
      return from_acc_code(default_options(), code);
    }

    /// Re-initialize NodeContext for all internal nodes in this tree.
    ///
    /// This is intended to be called after leaf IDs have been assigned
    /// (e.g., after parsing acceptance and running an ID-renaming pass),
    /// because shared-breakpoint contexts depend on `inf_leaf::id`.
    check_macrostate init_contexts() const;

    /// Compute successor macrostate(s) for this check tree node.
    std::vector<std::pair<check_macrostate, NodeContext>> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd,
      bool                              resample,
      NodeContext                       parent_context) const;

    /// Check whether this macrostate check tree is satisfied.
    bool is_satisfied() const;

    /// Gather relevant automaton states from this check tree.
    /// - `Fin` leaf: returns `safe`
    /// - `Inf` leaf: returns `track`
    /// - `And`/`Or` node: returns union of recursively gathered states
    std::set<unsigned> gather_states() const;

    /// Reduce this check tree by enforcing disjointness in `Or` nodes.
    ///
    /// For every `Or(left, right)` node, states gathered from `left` are
    /// removed from all leaf sets in `right`:
    /// - `Fin` leaf: removes from `safe`
    /// - `Inf` leaf: removes from `track` and `breakpoint`
    check_macrostate reduce() const;

  private:

    static check_macrostate fold(options_ptr opts, TreeType op, const std::vector<spot::acc_cond::acc_code>& parts);

    static check_macrostate from_acc_code_impl(options_ptr opts, const spot::acc_cond::acc_code& code);

    /**
     * @brief Convert `TreeType` to a stable string name.
     *
     * @param t Tree node type.
     * @return String name of the type.
     */
    static std::string tree_type_to_string(TreeType t);

    /**
     * @brief Recursive helper for printing the base tree.
     *
     * @param tree Tree to print.
     * @return String representation of @p tree.
     */
    static std::string to_string_impl(const base_tree& tree, bool show_shared_breakpoint);

  private:
    options_ptr opts_;
  };

  inline void collect_inf_leaf_ids(const check_macrostate& t, std::vector<unsigned>& out) {
    using base_tree = kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf>;
    const base_tree& bt = static_cast<const base_tree&>(t);
    if (bt.is_leaf()) {
      if (bt.type() == TreeType::Inf) {
        out.push_back(std::get<inf_leaf>(bt.leaf_value()).id);
      }
      return;
    }
    if(bt.type() != TreeType::And) {
      return;
    }
    collect_inf_leaf_ids(check_macrostate(nullptr, base_tree(bt.left())), out);
    collect_inf_leaf_ids(check_macrostate(nullptr, base_tree(bt.right())), out);
  }

  /**
   * @brief Check whether a check tree contains only Fin leaves and has no internal And nodes.
   *
   * This predicate returns true if every leaf in @p t is a `Fin` leaf and all
   * internal nodes (if any) are `Or` nodes. In other words, internal `And`
   * nodes are forbidden.
   *
   * The function is used by the OR-FIN optimization to decide whether a subtree
   * is eligible: the optimization requires all leaves to be `Fin` and that the
   * subtree contains no `And` internals so that violation-collection semantics
   * through `Or` nodes remain straightforward.
   *
   * @param t Check tree to inspect.
   * @return true if every leaf in @p t is a `Fin` leaf and no internal node is an `And`.
   */
  inline bool has_only_fin_leaves_no_inner_or(const check_macrostate& t) {
    using base_tree = kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf>;
    const base_tree& bt = static_cast<const base_tree&>(t);
    if (bt.is_leaf()) {
      return bt.type() == TreeType::Fin;
    }
    if (bt.type() == TreeType::And) return false;
    const check_macrostate left(nullptr, base_tree(bt.left()));
    const check_macrostate right(nullptr, base_tree(bt.right()));
    return has_only_fin_leaves_no_inner_or(left) && has_only_fin_leaves_no_inner_or(right);
  }

  inline bool AndOrNode::operator==(const AndOrNode& other) const {
    if (this->type != other.type)
      return false;
    return this->context == other.context;
    // Structural child comparison is handled by binary_tree::operator==
    // via *left_ and *right_; no need to compare the subtree here.
  }

  inline std::strong_ordering AndOrNode::operator<=>(const AndOrNode& other) const {
    if (this->type < other.type) {
      return std::strong_ordering::less;
    }
    if (other.type < this->type) {
      return std::strong_ordering::greater;
    }

    return (this->context <=> other.context);
    // Structural child ordering is handled by binary_tree::operator<
    // via *left_ and *right_; no subtree comparison needed here.
  }

/// partial macrostate for the given component
class mstate_sd_inductive : public abstract_complement_alg::mstate
{ // {{{
public: // DATA MEMBERS

  std::set<unsigned> check_ {};       // states for runs that need to be checked
  check_macrostate check_tree_; // check macrostate tree
  

public: // METHODS

  /// constructor
  mstate_sd_inductive(
    const std::set<unsigned>&  check,
    const check_macrostate&  check_tree
  ) : check_(check),
    check_tree_(check_tree)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return true; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_sd_inductive() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override {
    throw std::runtime_error("mstate_sd_inductive: get_breakpoint() should not be called");
    return this->check_; // to suppress compiler warning
  }

  virtual void set_breakpoint(const std::set<unsigned>&) override {
    throw std::runtime_error("mstate_sd_inductive: set_breakpoint() should not be called");
  }

  virtual bool subsum_less_early(const mstate& rhs) override {
    (void)rhs; // suppress unused parameter warning
    // TODO: implement subsumption for mstate_sd_inductive
    return false;
  };
}; // mstate_sd_inductive }}}
} // namespace sd_inductive

class complement_sd_inductive : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_sd_inductive(const cmpl_info& info, unsigned part_index);

  virtual mstate_set get_init() override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&,
    const mstate*,
    const bdd&) override {
    
    throw std::runtime_error("complement_sd_inductive: get_succ_track() should not be called");
    return {};
  }

  virtual mstate_set lift_track_to_active(const mstate* ms) override { 
    const sd_inductive::mstate_sd_inductive* src_ms = dynamic_cast<const sd_inductive::mstate_sd_inductive*>(ms);
    assert(src_ms);

    std::shared_ptr<mstate> cp(new sd_inductive::mstate_sd_inductive(
      src_ms->check_,
      src_ms->check_tree_
    ));
    return {cp};
  };

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol,
    bool resample = true) override;

  virtual bool use_round_robin() const override { return false; }

  virtual bool use_shared_breakpoint() const override { return false; }

  virtual spot::acc_cond get_acc_cond() override
  { return spot::acc_cond(1, spot::acc_cond::inf({0})); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~complement_sd_inductive() override { };

// AUXILIARY METHODS
public:

// DATA MEMBERS
private:
  spot::acc_cond::acc_code acc_cond_ {};
  sd_inductive::options_ptr opts_ { nullptr };
}; // complement_sd_inductive }}}



} // namespace kofola }}}
