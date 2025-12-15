// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

#include <compare>
#include <ostream>
#include <sstream>
#include <stdexcept>

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

  inline std::string set_to_string(const std::set<unsigned>& s) {
    std::ostringstream os;
    os << "{";
    bool first = true;
    for (const auto& x : s) {
      if (!first) {
        os << ",";
      }
      first = false;
      os << x;
    }
    os << "}";
    return os.str();
  }

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
      return "safe=" + set_to_string(this->safe) + ", color=" + mark_to_string(this->color);
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

    std::vector<check_macrostate> getSucc(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd) const;

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
      return "track=" + set_to_string(this->track) + ", breakpoint=" + set_to_string(this->breakpoint)
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

    std::vector<check_macrostate> getSucc(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd) const;
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
   * - leaf payload: `fin_leaf` or `inf_leaf`
   */
  class check_macrostate : public kofola::types::binary_tree<TreeType, fin_leaf, inf_leaf> {
    /** @brief Base tree type used for representation. */
    using base_tree = kofola::types::binary_tree<TreeType, fin_leaf, inf_leaf>;

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
    static check_macrostate fin(std::set<unsigned> safe) {
      return check_macrostate(base_tree::leaf(TreeType::Fin, fin_leaf{std::move(safe)}));
    }

    /**
     * @brief Build an `Inf` leaf.
     *
     * @param track      Set of tracked states.
     * @param breakpoint Breakpoint set.
     * @return A `check_macrostate` leaf of type `TreeType::Inf`.
     */
    static check_macrostate inf(std::set<unsigned> track, std::set<unsigned> breakpoint) {
      return check_macrostate(base_tree::leaf(TreeType::Inf, inf_leaf{std::move(track), std::move(breakpoint)}));
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
      base_tree l(std::move(left));
      base_tree r(std::move(right));
      return check_macrostate(base_tree::make_node(type, std::move(l), std::move(r)));
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

    /**
     * @brief Build a `check_macrostate` tree from a Spot acceptance formula.
     *
     * Supported grammar (no negation expected):
     * - conjunction: `&`  -> `TreeType::And`
     * - disjunction: `|`  -> `TreeType::Or`
    * - leaves: `Fin(m)`  -> `TreeType::Fin` and stores the mark @p m in `fin_leaf::color`
    * - leaves: `Inf(m)`  -> `TreeType::Inf` and stores the mark @p m in `inf_leaf::color`
     *
     * Any other operator will throw `std::invalid_argument`.
     */
    static check_macrostate from_acc_code(const spot::acc_cond::acc_code& code) {
      return from_acc_code_impl(code);
    }

    /**
     * @brief Compute successor macrostate(s) for this check tree node.
     *
     * Behavior (follows the implementation below):
     * - If this `check_macrostate` is a leaf, the call is delegated to the
     *   leaf payload's `getSucc` implementation, which computes successor
     *   checks for that specific leaf type.
     * - If this node is an `And` node, the successors of the left and right
     *   children are computed and the Cartesian product of those successor
     *   sets is returned; each resulting pair is combined into a new `And`
     *   node.
     * - If this node is an `Or` node, the provided `check_states` set is
     *   nondeterministically split into two parts (binary partition) and
     *   successors are computed for the left child using the first part and
     *   for the right child using the second part; the combinations are
     *   assembled into `Or` nodes and deduplicated before returning.
     *
     * The function propagates `aut`, `scc_info`, and `bdd` down to leaves
     * which perform the actual transition computation. It will throw a
     * `std::logic_error` if an unexpected internal node type is encountered
     * or if the nondeterministic split routine returns a non-binary
     * partition.
     *
     * @param aut Spot automaton pointer used for successor computation.
     * @param scc_info SCC decomposition information for `aut`.
     * @param check_states Subset of automaton states assigned to this
     *                     subtree (used when splitting states for `Or`).
     * @param bdd Shared BDD structure passed to leaf computations.
     * @return A vector of successor `check_macrostate` trees.
     */
    std::vector<check_macrostate> getSucc(
      const spot::const_twa_graph_ptr&  aut,
      const spot::scc_info&             scc_info,
      const std::set<unsigned>&         check_states,
      const bdd&                        bdd) const {
      if (this->is_leaf()) {
        return std::visit(
          [&](const auto& leaf) {
            return leaf.getSucc(aut, scc_info, check_states, bdd);
          },
          this->leaf_value());
      }

      const auto node_type = this->type();
      const check_macrostate left_ms(base_tree(this->left()));
      const check_macrostate right_ms(base_tree(this->right()));

      if (node_type == TreeType::And) {
        const auto left_succ = left_ms.getSucc(aut, scc_info, check_states, bdd);
        const auto right_succ = right_ms.getSucc(aut, scc_info, check_states, bdd);
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
            throw std::logic_error("check_macrostate::getSucc: nondet_split_set returned non-binary partition");
          }
          const auto& left_states = part[0];
          const auto& right_states = part[1];
          const auto left_succ = left_ms.getSucc(aut, scc_info, left_states, bdd);
          const auto right_succ = right_ms.getSucc(aut, scc_info, right_states, bdd);

          const auto combined = cartesian_product<check_macrostate, check_macrostate>(
            left_succ,
            right_succ,
            [](const check_macrostate& l, const check_macrostate& r) {
              return check_macrostate::make(TreeType::Or, l, r);
            });
          out.insert(combined.begin(), combined.end());
        }
        return std::vector<check_macrostate>(out.begin(), out.end());
      }

      throw std::logic_error("check_macrostate::getSucc: unexpected internal node type");
    }

  private:

    static check_macrostate fold(TreeType op, const std::vector<spot::acc_cond::acc_code>& parts) {
      if (parts.empty()) {
        throw std::invalid_argument("check_macrostate: empty And/Or in acceptance formula");
      }
      check_macrostate acc = from_acc_code_impl(parts.front());
      for (size_t i = 1; i < parts.size(); ++i) {
        acc = check_macrostate::make(op, std::move(acc), from_acc_code_impl(parts[i]));
      }
      return acc;
    }

    static check_macrostate from_acc_code_impl(const spot::acc_cond::acc_code& code) {
      if (code.empty()) {
        throw std::invalid_argument("check_macrostate: empty acceptance formula");
      }

      // Leaf: [mark][op]
      if (code.size() == 2) {
        const auto op = code[1].sub.op;
        if (op == spot::acc_cond::acc_op::Fin) {
          return check_macrostate(base_tree::leaf(TreeType::Fin, fin_leaf{{}, code[0].mark}));
        }
        if (op == spot::acc_cond::acc_op::Inf) {
          return check_macrostate(base_tree::leaf(TreeType::Inf, inf_leaf{{}, {}, code[0].mark}));
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

      // Some forms may not be caught above (e.g., parenthesized singletons); try to unwrap once.
      if (!conjuncts.empty() && conjuncts.size() == 1 && conjuncts[0] != code) {
        return from_acc_code_impl(conjuncts[0]);
      }
      if (!disjuncts.empty() && disjuncts.size() == 1 && disjuncts[0] != code) {
        return from_acc_code_impl(disjuncts[0]);
      }

      throw std::invalid_argument("check_macrostate: unsupported acceptance formula operator");
    }

    /**
     * @brief Convert `TreeType` to a stable string name.
     *
     * @param t Tree node type.
     * @return String name of the type.
     */
    static std::string tree_type_to_string(TreeType t) {
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

    /**
     * @brief Recursive helper for printing the base tree.
     *
     * @param tree Tree to print.
     * @return String representation of @p tree.
     */
    static std::string to_string_impl(const base_tree& tree) {
      const std::string head = tree_type_to_string(tree.type());
      if (tree.is_leaf()) {
        const std::string payload = std::visit(
          [](const auto& leaf) { return leaf.to_string(); },
          tree.leaf_value());
        return head + "(" + payload + ")";
      }

      return head + "(" + to_string_impl(tree.left()) + ", " + to_string_impl(tree.right()) + ")";
    }
  };

  inline std::vector<check_macrostate> fin_leaf::getSucc(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const std::set<unsigned>&         check_states,
    const bdd&                        bdd) const {
    (void)aut;
    (void)scc_info;
    (void)check_states;
    (void)bdd;
    // TODO: fill with actual implementation
    return {};
  }

  inline std::vector<check_macrostate> inf_leaf::getSucc(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const std::set<unsigned>&         check_states,
    const bdd&                        bdd) const {
    (void)aut;
    (void)scc_info;
    (void)check_states;
    (void)bdd;
    // TODO: fill with actual implementation
    return {};
  }



} // namespace sd_inductive


} // namespace kofola }}}
