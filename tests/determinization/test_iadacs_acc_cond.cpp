// Soundness tests for the acceptance condition shared by the IADACs
// determinization (kofola::determinize_iadacs) and the IADACs complementation
// (kofola::complement_init_almost_det), see kofola::determinisation_acc_cond.
//
// The end-to-end cases force the IADACs algorithms (det_based_on_iadac=yes),
// check that they are really used, and compare the result with the input
// (determinization) or with Spot's complement (complementation).

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <string>

// Spot headers
#include <spot/twaalgos/contains.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/word.hh>

// Kofola headers
#include "algorithms/complement_alg_iadacs.hpp"
#include "complement/complement_tela.hpp"
#include "determinization/determinize_tela.hpp"

// Test utilities
#include "../utils/test_utils.hpp"

namespace {

/// sets the TELA options and forces the IADACs algorithms for the lifetime of
/// the object, so that a failing REQUIRE does not leak the option
struct force_iadacs
{
    force_iadacs()
    {
        test_utils::setup_tela_options();
        kofola::OPTIONS.params["det_based_on_iadac"] = "yes";
    }
    ~force_iadacs() { kofola::OPTIONS.params.clear(); }
};

/// whether some state name of 'aut' contains 'needle'
bool has_state_named(const spot::twa_graph_ptr& aut, const std::string& needle)
{
    auto names = aut->get_named_prop<std::vector<std::string>>("state-names");
    if (!names) { return false; }
    for (const std::string& name : *names) {
        if (name.find(needle) != std::string::npos) { return true; }
    }
    return false;
}

/// whether the IADACs determinization / complementation is used on 'aut'
/// (raw=yes skips the postprocessing, which would drop the state names)
bool iadacs_det_used(const spot::twa_graph_ptr& aut)
{
    kofola::OPTIONS.params["raw"] = "yes";
    bool res = has_state_named(kofola::determinize_tela(aut), "IADACs-det");
    kofola::OPTIONS.params.erase("raw");
    return res;
}

bool iadacs_complement_used(const spot::twa_graph_ptr& aut)
{
    kofola::OPTIONS.params["raw"] = "yes";
    bool res = has_state_named(kofola::complement_tela(aut), "INIT_ALMOST_DET");
    kofola::OPTIONS.params.erase("raw");
    return res;
}

/// An initial state looping on a|b that jumps on c into all the 'n' states of
/// a deterministic SCC over a|b (c is not used inside, so the SCCs are IADACs).
/// In the SCC, a is a self-loop and b moves to the next state; only the
/// a-loop of the last state carries all the colours of Inf(0)&...&Inf(k-1).
/// On c.a^w all the 'n' runs are alive forever and only the last one, i.e.
/// the one with the largest label, is accepting.
spot::twa_graph_ptr make_many_runs_aut(unsigned n, unsigned k)
{
    spot::bdd_dict_ptr dict = spot::make_bdd_dict();
    spot::twa_graph_ptr aut = spot::make_twa_graph(dict);
    bdd p0 = bdd_ithvar(aut->register_ap("p0"));
    bdd p1 = bdd_ithvar(aut->register_ap("p1"));
    bdd a = !p0 & !p1;
    bdd b = p0 & !p1;
    bdd c = !p0 & p1;

    spot::acc_cond::acc_code code = spot::acc_cond::acc_code::t();
    spot::acc_cond::mark_t all_cols = {};
    for (unsigned i = 0; i < k; ++i) {
        code &= spot::acc_cond::acc_code::inf({i});
        all_cols.set(i);
    }
    aut->set_acceptance(k, code);

    aut->new_states(n + 1);
    aut->set_init_state(0);
    aut->new_edge(0, 0, a | b);
    for (unsigned i = 1; i <= n; ++i) {
        aut->new_edge(0, i, c);
        aut->new_edge(i, i, a, (i == n) ? all_cols : spot::acc_cond::mark_t{});
        aut->new_edge(i, i % n + 1, b);
    }
    return aut;
}

} // anonymous namespace

TEST_CASE("IADACs acceptance condition is compacted onto the used colours", "[iadacs]")
{
    // Inf(1)&Inf(3): colours 0 and 2 are not mentioned
    auto code = spot::acc_cond::acc_code::inf({1}) & spot::acc_cond::acc_code::inf({3});
    kofola::determinisation_acc_cond cond(code, 2);

    // two used colours plus a discontinuation colour per label, starting at 0
    CHECK(cond.get_acc_cond().num_sets() == 6);
    CHECK(cond.get_acc_cond().get_acceptance().used_sets().min_set() == 1);   // i.e. colour 0
    CHECK(cond.get_fin_mark(0) == 2);
    CHECK(cond.get_fin_mark(1) == 5);

    std::set<unsigned> cols;
    cond.add_colours({0, 1, 2, 3}, 1, cols);   // 0 and 2 are dropped, 1 -> 3, 3 -> 4
    CHECK(cols == std::set<unsigned>{3, 4});

    CHECK(kofola::determinisation_acc_cond::colours_needed(code, 2) == 6);
}

TEST_CASE("IADACs acceptance condition refuses more colours than Spot supports", "[iadacs]")
{
    // Silently dropping the colours beyond the limit used to make the last
    // labels of the condition unsound; the callers now fall back instead.
    auto code = spot::acc_cond::acc_code::inf({0}) & spot::acc_cond::acc_code::inf({1}) &
                spot::acc_cond::acc_code::inf({2}) & spot::acc_cond::acc_code::inf({3});
    CHECK(kofola::determinisation_acc_cond::colours_needed(code, 7) == 35);
    CHECK_THROWS_AS(kofola::determinisation_acc_cond(code, 7), std::runtime_error);
    CHECK_NOTHROW(kofola::determinisation_acc_cond(code, 6));   // 30 colours
}

TEST_CASE("IADACs with an acceptance that does not mention colour 0", "[iadacs][determinize]")
{
    force_iadacs guard;

    // Regression test: the colours of the determinization were rebased by the
    // smallest colour of the condition (here 1), but the condition itself was
    // not, so Inf(1) ended up testing the discontinuation colour of the first
    // label.  The determinization then accepted a different language.
    spot::twa_graph_ptr aut =
        test_utils::load_automaton_from_file("tests/test_data/iadac_min_colour.hoa");
    REQUIRE(aut != nullptr);

    SECTION("determinization") {
        REQUIRE(iadacs_det_used(aut));
        CHECK(test_utils::test_determinize_equivalence(aut) ==
              test_utils::determinize_check::passed);
    }
    SECTION("complementation") {
        REQUIRE(iadacs_complement_used(aut));
        CHECK(test_utils::test_complement_equivalence(aut));
    }
}

TEST_CASE("IADACs with a transition colour the acceptance does not mention", "[iadacs][determinize]")
{
    force_iadacs guard;

    // Regression test: the acceptance is Inf(0), but some transitions carry
    // colour 1.  The discontinuation colour of the first label was 1 as well,
    // so an accepting run that also sees colour 1 infinitely often was taken
    // as discontinued infinitely often, and the word was rejected.
    spot::twa_graph_ptr aut =
        test_utils::load_automaton_from_file("tests/test_data/iadac_unused_colour.hoa");
    REQUIRE(aut != nullptr);

    SECTION("determinization") {
        REQUIRE(iadacs_det_used(aut));
        CHECK(test_utils::test_determinize_equivalence(aut) ==
              test_utils::determinize_check::passed);
    }
    SECTION("complementation") {
        REQUIRE(iadacs_complement_used(aut));
        CHECK(test_utils::test_complement_equivalence(aut));
    }
}

TEST_CASE("IADACs whose colours exceed Spot's limit", "[iadacs][determinize]")
{
    force_iadacs guard;

    // Regression test: 7 runs times (4 colours + 1 discontinuation colour)
    // needs 35 colours.  The shift of the last label (30) is below the limit,
    // so Spot did not complain and silently dropped the colours 32-34; the
    // determinization then rejected c.a^w, which is accepted by the last run
    // only.  Now IADACs is not used; the determinization may still fail for
    // lack of colours, but must never return a wrong automaton.
    spot::twa_graph_ptr aut = make_many_runs_aut(7, 4);
    const std::string word = "!p0&p1;cycle{!p0&!p1}";   // c.a^w
    REQUIRE(aut->intersects(spot::parse_word(word, aut->get_dict())->as_automaton()));

    SECTION("determinization") {
        spot::twa_graph_ptr det;
        try {
            det = kofola::determinize_tela(aut);
        } catch (const std::runtime_error&) {
            SUCCEED("too many colours for every algorithm");
        }
        if (det) {
            CHECK(spot::is_deterministic(det));
            CHECK(spot::are_equivalent(det, aut));
        }
    }
    SECTION("complementation") {
        CHECK(test_utils::test_complement_equivalence(aut));
    }
}
