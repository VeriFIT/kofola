// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

#include <compare>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "../types/binary_tree.hpp"

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



} // namespace sd_inductive


} // namespace kofola }}}
