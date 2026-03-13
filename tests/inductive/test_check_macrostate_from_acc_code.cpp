#include <catch2/catch_test_macros.hpp>

#include "algorithms/complement_alg_sd_inductive.hpp"

// SPOT
#include <spot/twa/acc.hh>

namespace {
using kofola::sd_inductive::TreeType;
using kofola::sd_inductive::check_macrostate;
using kofola::sd_inductive::AndOrNode;
using kofola::sd_inductive::fin_leaf;
using kofola::sd_inductive::inf_leaf;
using kofola::sd_inductive::options;
using kofola::sd_inductive::options_ptr;
using kofola::sd_inductive::has_only_fin_leaves_no_inner_or;
using base_tree = kofola::types::binary_tree<TreeType, AndOrNode, fin_leaf, inf_leaf>;

static const auto& as_base(const check_macrostate& t) {
  return static_cast<const base_tree&>(t);
}

static check_macrostate fin_color(std::initializer_list<unsigned> idx, unsigned id = 0) {
  const spot::acc_cond::mark_t m(idx.begin(), idx.end());
  return check_macrostate::fin({}, m, id);
}

static check_macrostate inf_color(std::initializer_list<unsigned> idx, unsigned id = 0) {
  const spot::acc_cond::mark_t m(idx.begin(), idx.end());
  return check_macrostate::inf({}, {}, m, id);
}

struct leaf_sig {
  TreeType type;
  spot::acc_cond::mark_t color;
  bool operator==(const leaf_sig&) const = default;
};

static bool leaf_sig_less(const leaf_sig& a, const leaf_sig& b) {
  if (a.type != b.type)
    return a.type < b.type;
  return a.color < b.color;
}

static void collect(const base_tree& t,
                    std::vector<TreeType>& internal_types,
                    std::vector<leaf_sig>& leaves) {
  if (t.is_leaf()) {
    if (t.type() == TreeType::Fin) {
      const auto& f = std::get<fin_leaf>(t.leaf_value());
      REQUIRE(f.safe.empty());
      leaves.push_back(leaf_sig{TreeType::Fin, f.color});
      return;
    }
    if (t.type() == TreeType::Inf) {
      const auto& i = std::get<inf_leaf>(t.leaf_value());
      REQUIRE(i.track.empty());
      REQUIRE(i.breakpoint.empty());
      leaves.push_back(leaf_sig{TreeType::Inf, i.color});
      return;
    }
    FAIL("Unexpected leaf node type");
  }

  internal_types.push_back(t.type());
  collect(t.left(), internal_types, leaves);
  collect(t.right(), internal_types, leaves);
}
} // namespace

TEST_CASE("check_macrostate builds leaves from Spot acc_code", "[check_macrostate]") {
  SECTION("Inf leaf") {
    auto code = spot::acc_cond::acc_code("Inf(0)");
    auto got = check_macrostate::from_acc_code(code);
    auto exp = inf_color({0}, 0);
    REQUIRE(as_base(got) == as_base(exp));
  }

  SECTION("Fin leaf") {
    auto code = spot::acc_cond::acc_code("Fin(2)");
    auto got = check_macrostate::from_acc_code(code);
    auto exp = fin_color({2}, 0);
    REQUIRE(as_base(got) == as_base(exp));
  }
}

TEST_CASE("check_macrostate builds And/Or structure from Spot acc_code", "[check_macrostate]") {
  SECTION("Conjunction folds into And nodes") {
    auto code = spot::acc_cond::acc_code("Inf(0) & Fin(1) & Inf(2)");
    auto got = check_macrostate::from_acc_code(code);

    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);

    REQUIRE(internal_types.size() == 2);
    REQUIRE(std::all_of(internal_types.begin(), internal_types.end(), [](TreeType t) { return t == TreeType::And; }));

    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected {
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{1}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{0}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{2}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }

  SECTION("Disjunction folds into Or nodes") {
    auto code = spot::acc_cond::acc_code("Fin(0) | Inf(1) | Fin(3)");
    auto got = check_macrostate::from_acc_code(code);

    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);

    REQUIRE(internal_types.size() == 2);
    REQUIRE(std::all_of(internal_types.begin(), internal_types.end(), [](TreeType t) { return t == TreeType::Or; }));

    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected {
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{0}},
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{3}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{1}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }

  SECTION("Nested formulas") {
    auto code = spot::acc_cond::acc_code("(Fin(0) & Inf(1)) | (Fin(2) & Inf(3))");
    auto got = check_macrostate::from_acc_code(code);

    // Spot may reorder conjuncts/disjuncts depending on version/platform.
    // Verify the parsed tree shape and leaf multiset, ignoring IDs and ordering.
    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);

    REQUIRE(internal_types.size() == 3);
    REQUIRE(std::count(internal_types.begin(), internal_types.end(), TreeType::Or) == 1);
    REQUIRE(std::count(internal_types.begin(), internal_types.end(), TreeType::And) == 2);

    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected {
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{0}},
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{2}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{1}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{3}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }
}

// Helpers for reorganization tests.
namespace {
static check_macrostate child_left(const check_macrostate& t) {
  return check_macrostate(t.get_options_ptr(), base_tree(as_base(t).left()));
}
static check_macrostate child_right(const check_macrostate& t) {
  return check_macrostate(t.get_options_ptr(), base_tree(as_base(t).right()));
}
static options_ptr or_fin_opts() {
  return std::make_shared<const options>(options{.use_or_fin_opt = true});
}
} // namespace

TEST_CASE("reorganize_fins_left: FIN-only subtrees moved left in OR (use_or_fin_opt=true)",
          "[check_macrostate][or_fin_opt][reorganize_fins_left]") {

  SECTION("Inf(0) | Fin(1): FIN leaf moves to the left of OR") {
    auto code = spot::acc_cond::acc_code("Inf(0) | Fin(1)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);

    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Or);
    // Left subtree must be FIN-only after reorganization.
    REQUIRE(has_only_fin_leaves_no_inner_or(child_left(got)));

    // Leaf multiset unchanged.
    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);
    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected{
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{1}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{0}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }

  SECTION("Inf(0) | Fin(1) | Fin(2): both FINs grouped on left, INF on right (user example)") {
    // Corresponds to the user's example:
    //   (INF1 | FIN1) | FIN3  =>  (FIN1 | FIN3) | INF1
    auto code = spot::acc_cond::acc_code("Inf(0) | Fin(1) | Fin(2)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);

    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Or);

    // Top-level left child must be FIN-only (Or(Fin(1), Fin(2))).
    REQUIRE(has_only_fin_leaves_no_inner_or(child_left(got)));

    // Top-level right child must be the single INF leaf.
    auto right = child_right(got);
    REQUIRE(right.is_leaf());
    REQUIRE(right.type() == TreeType::Inf);

    // Leaf multiset unchanged.
    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);
    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected{
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{1}},
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{2}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{0}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }

  SECTION("Fin(0) | Inf(1) | Fin(2): both FINs again grouped on left") {
    auto code = spot::acc_cond::acc_code("Fin(0) | Inf(1) | Fin(2)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);

    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Or);

    REQUIRE(has_only_fin_leaves_no_inner_or(child_left(got)));
    auto right = child_right(got);
    REQUIRE(right.is_leaf());
    REQUIRE(right.type() == TreeType::Inf);

    std::vector<TreeType> internal_types;
    std::vector<leaf_sig> leaves;
    collect(as_base(got), internal_types, leaves);
    std::sort(leaves.begin(), leaves.end(), leaf_sig_less);
    std::vector<leaf_sig> expected{
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{0}},
      leaf_sig{TreeType::Fin, spot::acc_cond::mark_t{2}},
      leaf_sig{TreeType::Inf, spot::acc_cond::mark_t{1}},
    };
    std::sort(expected.begin(), expected.end(), leaf_sig_less);
    REQUIRE(leaves == expected);
  }

  SECTION("(Fin(0) & Inf(1)) | Fin(2): FIN leaf moves left, AND node stays right") {
    // And(Fin(0), Inf(1)) is NOT FIN-only; Fin(2) IS.
    auto code = spot::acc_cond::acc_code("(Fin(0) & Inf(1)) | Fin(2)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);

    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Or);

    // Left must be FIN-only (the plain Fin(2) leaf).
    auto left = child_left(got);
    REQUIRE(has_only_fin_leaves_no_inner_or(left));
    REQUIRE(left.is_leaf());
    REQUIRE(left.type() == TreeType::Fin);
    REQUIRE(std::get<fin_leaf>(as_base(left).leaf_value()).color == spot::acc_cond::mark_t{2});

    // Right must be the AND node.
    auto right = child_right(got);
    REQUIRE(!right.is_leaf());
    REQUIRE(right.type() == TreeType::And);
  }

  SECTION("Fin(0) | Fin(1): all-FIN OR leaves both children FIN-only") {
    auto code = spot::acc_cond::acc_code("Fin(0) | Fin(1)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);

    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Or);
    // Both children remain FIN-only (nothing to move).
    REQUIRE(has_only_fin_leaves_no_inner_or(child_left(got)));
    REQUIRE(has_only_fin_leaves_no_inner_or(child_right(got)));
  }

  SECTION("Inf(0): single INF leaf unaffected") {
    auto code = spot::acc_cond::acc_code("Inf(0)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);
    REQUIRE(as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Inf);
  }

  SECTION("Fin(0): single FIN leaf unaffected") {
    auto code = spot::acc_cond::acc_code("Fin(0)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);
    REQUIRE(as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::Fin);
  }

  SECTION("Fin(0) & Inf(1): AND at top level is preserved as AND") {
    auto code = spot::acc_cond::acc_code("Fin(0) & Inf(1)");
    auto got = check_macrostate::from_acc_code(or_fin_opts(), code);
    REQUIRE(!as_base(got).is_leaf());
    REQUIRE(as_base(got).type() == TreeType::And);
  }
}
