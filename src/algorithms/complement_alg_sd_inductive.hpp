// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

#include <compare>
#include <memory>
#include <ostream>
#include <sstream>
#include <set>
#include <string>
#include <vector>

#include "../types/binary_tree.hpp"
#include "../util/sets.hpp"

// SPOT
#include <spot/twa/acc.hh>

namespace kofola { // {{{

namespace sd_inductive {

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

  class check_macrostate;

  /**
   * @brief Payload of an internal `And`/`Or` node in `check_macrostate`.
   *
   * Stores a subtree (rooted at the corresponding internal node) as a
   * `check_macrostate`. This is stored indirectly to avoid recursive
   * by-value type definitions.
   */
  struct AndOrNode {
    TreeType type {TreeType::And};
    std::shared_ptr<check_macrostate> subtree {};

    AndOrNode() = default;
    AndOrNode(TreeType t, check_macrostate subtree_);

    bool operator==(const AndOrNode& other) const;
    std::strong_ordering operator<=>(const AndOrNode& other) const;
  };

  /**
   * @brief Payload of a `Fin` leaf in `check_macrostate`.
   */
  struct fin_leaf {

    /** @brief Set of safe states relevant for the Fin-check. */
    std::set<unsigned> safe {};
    /** @brief Acceptance color associated with this Fin-check. */
    spot::acc_cond::mark_t color {}; 

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

    std::vector<check_macrostate> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd,
      bool                              resample) const;

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

    std::vector<check_macrostate> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd,
      bool                              resample) const;

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
     * @param tree Underlying tree representation.
     */
    explicit check_macrostate(base_tree tree) : base_tree(std::move(tree)) {}

    /**
     * @brief Build a `Fin` leaf.
     *
     * @param safe Set of safe states.
     * @return A `check_macrostate` leaf of type `TreeType::Fin`.
     */
    static check_macrostate fin(std::set<unsigned> safe, spot::acc_cond::mark_t color) {
      return check_macrostate(base_tree::leaf(TreeType::Fin, fin_leaf{std::move(safe), color}));
    }

    /**
     * @brief Build an `Inf` leaf.
     *
     * @param track      Set of tracked states.
     * @param breakpoint Breakpoint set.
     * @return A `check_macrostate` leaf of type `TreeType::Inf`.
     */
    static check_macrostate inf(std::set<unsigned> track, std::set<unsigned> breakpoint, spot::acc_cond::mark_t color) {
      return check_macrostate(base_tree::leaf(TreeType::Inf, inf_leaf{std::move(track), std::move(breakpoint), color}));
    }

    /**
     * @brief Build an internal node with two children.
     *
     * @param type  Node type (typically `TreeType::And` or `TreeType::Or`).
     * @param left  Left subtree.
     * @param right Right subtree.
     * @return A `check_macrostate` internal node.
     */
    static check_macrostate make(TreeType type, check_macrostate left, check_macrostate right) {
      // Build a concrete subtree rooted at this internal node.
      // We intentionally construct this subtree using a default internal
      // payload, and store it in the node payload for now.
      check_macrostate subtree_root(base_tree::make_node(type, base_tree(left), base_tree(right)));
      base_tree l(std::move(left));
      base_tree r(std::move(right));
      return check_macrostate(base_tree::make_node(type, AndOrNode(type, std::move(subtree_root)), std::move(l), std::move(r)));
    }

    /**
     * @brief Convert the whole check tree into a human-readable string.
     *
     * @return String representation of this macrostate check tree.
     */
    std::string to_string() const {
      return to_string_impl(static_cast<const base_tree&>(*this));
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
    static check_macrostate from_acc_code(const spot::acc_cond::acc_code& code);

    /// Compute successor macrostate(s) for this check tree node.
    std::vector<check_macrostate> get_succ(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd,
      bool                              resample) const;

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

    static check_macrostate fold(TreeType op, const std::vector<spot::acc_cond::acc_code>& parts);

    static check_macrostate from_acc_code_impl(const spot::acc_cond::acc_code& code);

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
    static std::string to_string_impl(const base_tree& tree);
  };

  inline AndOrNode::AndOrNode(TreeType t, check_macrostate subtree_)
    : type(t),
      subtree(std::make_shared<check_macrostate>(std::move(subtree_))) {}

  inline bool AndOrNode::operator==(const AndOrNode& other) const {
    if (this->type != other.type) {
      return false;
    }
    if (!this->subtree && !other.subtree) {
      return true;
    }
    if (!this->subtree || !other.subtree) {
      return false;
    }
    return *this->subtree == *other.subtree;
  }

  inline std::strong_ordering AndOrNode::operator<=>(const AndOrNode& other) const {
    if (this->type < other.type) {
      return std::strong_ordering::less;
    }
    if (other.type < this->type) {
      return std::strong_ordering::greater;
    }

    // Same type: order by subtree presence then subtree structure.
    if (!this->subtree && !other.subtree) {
      return std::strong_ordering::equal;
    }
    if (!this->subtree) {
      return std::strong_ordering::less;
    }
    if (!other.subtree) {
      return std::strong_ordering::greater;
    }
    if (*this->subtree < *other.subtree) {
      return std::strong_ordering::less;
    }
    if (*other.subtree < *this->subtree) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equal;
  }

enum class mstate_type {
  GUESS,
  CHECK,
};

/// partial macrostate for the given component
class mstate_sd_inductive : public abstract_complement_alg::mstate
{ // {{{
public: // DATA MEMBERS

  std::set<unsigned> check_ {};       // states for runs that need to be checked
  check_macrostate check_tree_; // check macrostate tree
  mstate_type type_ {};          // type of the macrostate (GUESS / CHECK)
  

public: // METHODS

  /// constructor
  mstate_sd_inductive(
    const std::set<unsigned>&  check,
    const check_macrostate&  check_tree,
    const mstate_type&  type
  ) : check_(check),
    check_tree_(check_tree),
    type_(type)
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
      src_ms->check_tree_,
      src_ms->type_
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
}; // complement_sd_inductive }}}



} // namespace kofola }}}
