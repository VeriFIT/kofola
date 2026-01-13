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
