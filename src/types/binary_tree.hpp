#pragma once

#include <concepts>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <variant>

namespace kofola::types {
/** Concept for values that support total ordering (`==`, `<`, and friends). */
template <class T>
concept totally_ordered = std::totally_ordered<T>;

/**
 * Recursive binary tree with heterogeneous leaves and explicit per-node type.
 *
 * Representation:
 * - `kind_` distinguishes leaf vs internal node.
 * - Leaves:
 *   - `leaf_value_` stores the leaf payload (one of `LeafTypes...`).
 * - Internal nodes:
 *   - `type_` stores the node's type (e.g. operator/kind). (Also present on leaves.)
 *   - `left_`/`right_` store subtrees.
 *
 * Comparison:
 * - Ordered by (type_, kind_, payload/children)
 * - For leaves, payload is ordered by `std::variant` (index then value)
 * - For internal nodes, children are ordered lexicographically (left, right)
 */
template <class NodeType, class... LeafTypes>
class binary_tree {
  static_assert(sizeof...(LeafTypes) > 0, "binary_tree requires at least one leaf type");
  static_assert(totally_ordered<NodeType>, "NodeType must be totally ordered");
  static_assert((totally_ordered<LeafTypes> && ...), "All leaf types must be totally ordered");

public:
  /** Variant holding the payload of a leaf. */
  using leaf_variant = std::variant<LeafTypes...>;

  /** Discriminant of the node representation. */
  enum class kind { leaf, node };

private:
  kind kind_;
  NodeType type_;
  std::optional<leaf_variant> leaf_value_;
  std::unique_ptr<binary_tree> left_;
  std::unique_ptr<binary_tree> right_;

public:
  /** No default construction; a node must be a leaf or an internal node. */
  binary_tree() = delete;

  /**
   * Construct a leaf from a pre-built leaf variant.
   *
   * @param type Logical type associated with this leaf node.
   * @param v    Leaf payload.
   */
  static binary_tree leaf(NodeType type, leaf_variant v) {
    return binary_tree(std::move(type), std::move(v));
  }

  /**
   * Construct a leaf from a concrete payload type.
   *
   * This overload only participates when `T` is one of `LeafTypes...`.
   *
   * @param type Logical type associated with this leaf node.
   * @param v    Leaf payload.
   */
  template <class T>
    requires (std::same_as<std::remove_cvref_t<T>, LeafTypes> || ...)
  static binary_tree leaf(NodeType type, T&& v) {
    return binary_tree(std::move(type), leaf_variant{std::forward<T>(v)});
  }

  /**
   * Construct an internal node with two children.
   *
   * @param type  Logical type associated with this internal node.
   * @param left  Left subtree.
   * @param right Right subtree.
   */
  static binary_tree make_node(NodeType type, binary_tree left, binary_tree right) {
    return binary_tree(std::move(type), std::move(left), std::move(right));
  }

  /** Deep-copy constructor (recursively copies the subtree). */
  binary_tree(const binary_tree& other)
      : kind_(other.kind_),
        type_(other.type_),
        leaf_value_(other.leaf_value_),
        left_(other.left_ ? std::make_unique<binary_tree>(*other.left_) : nullptr),
        right_(other.right_ ? std::make_unique<binary_tree>(*other.right_) : nullptr) {}

  /** Move constructor (defaulted). */
  binary_tree(binary_tree&&) noexcept = default;

  /** Deep-copy assignment (recursively copies the subtree). */
  binary_tree& operator=(const binary_tree& other) {
    if (this == &other) {
      return *this;
    }
    *this = binary_tree(other);
    return *this;
  }

  /** Move assignment (defaulted). */
  binary_tree& operator=(binary_tree&&) noexcept = default;

  /** Returns whether this object is a leaf or an internal node. */
  kind get_kind() const { return kind_; }

  /** Returns true iff this node is a leaf. */
  bool is_leaf() const { return kind_ == kind::leaf; }

  /** Returns true iff this node is an internal node. */
  bool is_node() const { return kind_ == kind::node; }

  /** Returns the logical node type (present for both leaves and internal nodes). */
  const NodeType& type() const { return type_; }

  /** Returns the leaf payload (only valid when `is_leaf() == true`). */
  const leaf_variant& leaf_value() const { return *leaf_value_; }

  /** Returns the left subtree (only valid when `is_node() == true`). */
  const binary_tree& left() const { return *left_; }

  /** Returns the right subtree (only valid when `is_node() == true`). */
  const binary_tree& right() const { return *right_; }

  /** Structural equality comparison (compares type, kind, and payload/children). */
  friend bool operator==(const binary_tree& a, const binary_tree& b) {
    if (a.kind_ != b.kind_ || !(a.type_ == b.type_)) {
      return false;
    }
    if (a.is_leaf()) {
      return a.leaf_value_ == b.leaf_value_;
    }
    return *a.left_ == *b.left_ && *a.right_ == *b.right_;
  }

  /** Negation of `operator==`. */
  friend bool operator!=(const binary_tree& a, const binary_tree& b) { return !(a == b); }

  /**
   * Strict weak ordering.
   *
   * Orders first by `type()`, then by kind (leaf before node), then by leaf
   * payload or lexicographically by children.
   */
  friend bool operator<(const binary_tree& a, const binary_tree& b) {
    if (a.type_ < b.type_) {
      return true;
    }
    if (b.type_ < a.type_) {
      return false;
    }

    // Same type: order by kind (leaf before node).
    if (a.kind_ != b.kind_) {
      return a.kind_ == kind::leaf;
    }

    if (a.is_leaf() && b.is_leaf()) {
      // variant orders by index then value.
      return *a.leaf_value_ < *b.leaf_value_;
    }

    // Both are nodes (and same type).
    if (*a.left_ < *b.left_) {
      return true;
    }
    if (*b.left_ < *a.left_) {
      return false;
    }
    return *a.right_ < *b.right_;
  }

  /** Derived ordering relation. */
  friend bool operator<=(const binary_tree& a, const binary_tree& b) { return !(b < a); }

  /** Derived ordering relation. */
  friend bool operator>(const binary_tree& a, const binary_tree& b) { return b < a; }

  /** Derived ordering relation. */
  friend bool operator>=(const binary_tree& a, const binary_tree& b) { return !(a < b); }

private:
  /** Private constructor for building leaf nodes. */
  binary_tree(NodeType type, leaf_variant v)
      : kind_(kind::leaf),
        type_(std::move(type)),
        leaf_value_(std::move(v)),
        left_(nullptr),
        right_(nullptr) {}

  /** Private constructor for building internal nodes. */
  binary_tree(NodeType type, binary_tree left, binary_tree right)
      : kind_(kind::node),
        type_(std::move(type)),
        leaf_value_(std::nullopt),
        left_(std::make_unique<binary_tree>(std::move(left))),
        right_(std::make_unique<binary_tree>(std::move(right))) {}
};

} // namespace kofola::types
