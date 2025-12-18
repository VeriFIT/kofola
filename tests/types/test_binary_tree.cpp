#include <catch2/catch_test_macros.hpp>

#include "types/binary_tree.hpp"

#include <string>

namespace {

enum class op_t { and_, or_ };

using tree_t = kofola::types::binary_tree<op_t, int, std::string>;

} // namespace

TEST_CASE("binary_tree: leaf equality and ordering") {
  const auto a = tree_t::leaf(op_t::and_, 1);
  const auto b = tree_t::leaf(op_t::and_, 2);
  const auto c = tree_t::leaf(op_t::and_, std::string{"x"});

  REQUIRE(a == tree_t::leaf(op_t::and_, 1));
  REQUIRE(a != b);

  // Same type: leaf payload orders by variant (int before string here)
  REQUIRE(a < c);
  REQUIRE(b < c);
  REQUIRE(a < b);
}

TEST_CASE("binary_tree: node equality and ordering") {
  const auto l1 = tree_t::leaf(op_t::and_, 1);
  const auto l2 = tree_t::leaf(op_t::and_, 2);

  const auto n1 = tree_t::make_node(op_t::and_, l1, l2);
  const auto n2 = tree_t::make_node(op_t::and_, tree_t::leaf(op_t::and_, 1), tree_t::leaf(op_t::and_, 2));
  const auto n3 = tree_t::make_node(op_t::or_, tree_t::leaf(op_t::or_, 1), tree_t::leaf(op_t::or_, 2));
  const auto n4 = tree_t::make_node(op_t::and_, tree_t::leaf(op_t::and_, 2), tree_t::leaf(op_t::and_, 1));

  REQUIRE(n1 == n2);
  REQUIRE(n1 != n3);

  // order by type then kind then payload/children
  REQUIRE(n1 < n3);
  REQUIRE(n1 < n4);

  // leaf < node
  REQUIRE(l1 < n1);
}
