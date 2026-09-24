// Tests for the modular determinization of TELA (kofola::determinize_tela).
//
// Inherently weak and deterministic partitions are supported, which together
// cover every elevator automaton; a nondeterministic accepting SCC is first
// removed by limit-determinization (kofola::Elevatorization).  The first test
// case runs over a handful of weak (nondeterministic) automata, the second one
// over every elevator automaton of test_data; both check that the result is
// deterministic and language equivalent to the input.  A third one covers a
// non-elevator automaton, i.e. the elevatorization path.

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <string>

// Spot headers
#include <spot/twaalgos/isdet.hh>

// Kofola headers
#include "determinization/determinize_tela.hpp"
#include "util/helpers.hpp"

// Test utilities
#include "../utils/test_utils.hpp"

TEST_CASE("determinization of weak automata", "[determinize]")
{
    test_utils::setup_tela_options();

    for (const std::string& filename : test_utils::WEAK_TEST_FILES) {
        SECTION("file: " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            INFO("determinizing " << filename);
            CHECK(test_utils::test_determinize_equivalence(aut) ==
                  test_utils::determinize_check::passed);
        }
    }
}

TEST_CASE("determinization of every elevator automaton of test_data", "[determinize]")
{
    test_utils::setup_tela_options();

    const std::vector<std::string> files = test_utils::list_test_data_automata();
    REQUIRE(!files.empty());

    unsigned checked = 0;
    unsigned unverifiable = 0;
    unsigned non_elevator = 0;

    for (const std::string& filename : files) {
        spot::twa_graph_ptr aut = test_utils::load_automaton_exact_path(filename);
        if (!aut) { continue; }

        // non-elevator automata go through elevatorization first and are
        // covered by the next test case
        if (!test_utils::is_elevator_automaton(aut)) { ++non_elevator; continue; }

        INFO("determinizing " << filename);
        test_utils::determinize_check res = test_utils::test_determinize_equivalence(aut);

        // Spot cannot complement a couple of the automata of the test data, so
        // their equivalence cannot be decided; the determinization itself and
        // the determinism of its result are still checked above
        if (test_utils::determinize_check::unverifiable == res) { ++unverifiable; continue; }

        CHECK(res == test_utils::determinize_check::passed);
        ++checked;
    }

    std::cout << "determinized " << checked << " elevator automata of test_data ("
              << unverifiable << " unverifiable, "
              << non_elevator << " non-elevator ones skipped)" << std::endl;
    CHECK(checked > 0);
}

TEST_CASE("IADACs determinization of every elevator automaton with an IADAC", "[determinize]")
{
    test_utils::setup_tela_options();
    // by default, the IADACs algorithm gives way to the DAC one when its
    // colours might not fit into Spot's budget; force it here
    kofola::OPTIONS.params["det_based_on_iadac"] = "yes";

    unsigned checked = 0;
    for (const std::string& filename : test_utils::list_test_data_automata()) {
        spot::twa_graph_ptr aut = test_utils::load_automaton_exact_path(filename);
        if (!aut || !test_utils::is_elevator_automaton(aut)) { continue; }

        spot::scc_info si(aut, spot::scc_info_options::ALL);
        si.determine_unknown_acceptance();
        std::string scc_types = helpers::get_scc_types(si);
        bool has_iadac = false;
        for (unsigned scc = 0; scc < si.scc_count(); ++scc) {
            has_iadac |= helpers::is_accepting_initial_almost_detscc(scc_types, scc);
        }
        if (!has_iadac) { continue; }

        INFO("determinizing " << filename);
        test_utils::determinize_check res = test_utils::test_determinize_equivalence(aut);
        if (test_utils::determinize_check::unverifiable == res) { continue; }

        CHECK(res == test_utils::determinize_check::passed);
        ++checked;
    }

    std::cout << "IADACs-determinized " << checked << " elevator automata of test_data" << std::endl;
    CHECK(checked > 0);
    kofola::OPTIONS.params.erase("det_based_on_iadac");
}

TEST_CASE("IADACs determinization with a gap in the colours", "[determinize]")
{
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["det_based_on_iadac"] = "yes";

    // Regression test: the acceptance is Inf(0)&Inf(2), colour 1 is unused.
    // determinisation_acc_cond used to size the colour block of a run by the
    // number of used colours (3) rather than by the largest one (4), so the
    // discontinuation colour of the first run coincided with colour 0 of the
    // second one.  On c(ab)^w the first run is accepting and the second one
    // sees colour 0 only, which then discontinued the first run infinitely
    // often and the word was rejected.
    spot::twa_graph_ptr aut =
        test_utils::load_automaton_from_file("tests/test_data/iadac_colour_gap.hoa");
    REQUIRE(aut != nullptr);

    CHECK(test_utils::test_determinize_equivalence(aut) ==
          test_utils::determinize_check::passed);
    kofola::OPTIONS.params.erase("det_based_on_iadac");
}

TEST_CASE("determinization of a non-elevator automaton limit-determinizes first", "[determinize]")
{
    test_utils::setup_tela_options();

    // Regression test for a heap-buffer-overflow in Elevatorization: the
    // acceptance condition of this automaton is a conjunction of eight clauses,
    // so its DNF has several disjuncts and limit_deter() calls
    // create_deter_part() once per disjunct.  Every call adds states to the
    // automaton, which lie outside the range of the scc_info built in the
    // constructor, and the later calls used to hand those fresh states to
    // scc_info::scc_of().  Only visible under a sanitizer or as a segfault.
    //
    // Only a single automaton is checked here: determinization blows up on
    // several of the other non-elevator automata of the test data, so a sweep
    // over all of them (the counterpart of the elevator test case above) would
    // not terminate.
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(
        "tests/test_data/tela_det/nac_multi_disjunct_acc.hoa");
    REQUIRE(aut != nullptr);
    REQUIRE(!test_utils::is_elevator_automaton(aut));   // has a nondet. accepting SCC

    spot::twa_graph_ptr det = kofola::determinize_tela(aut);
    REQUIRE(det != nullptr);
    CHECK(spot::is_deterministic(det));
}

TEST_CASE("determinization of a deterministic automaton is the identity", "[determinize]")
{
    test_utils::setup_tela_options();

    spot::twa_graph_ptr aut =
        test_utils::load_automaton_from_file("tests/test_data/inf_a.hoa");
    REQUIRE(aut != nullptr);
    REQUIRE(spot::is_deterministic(aut));

    spot::twa_graph_ptr det = kofola::determinize_tela(aut);
    REQUIRE(det != nullptr);
    CHECK(det == aut);
}

TEST_CASE("determinization of an automaton with no accepting SCC is empty", "[determinize]")
{
    test_utils::setup_tela_options();

    // an automaton with a single nonaccepting SCC
    spot::bdd_dict_ptr dict = spot::make_bdd_dict();
    spot::twa_graph_ptr aut = spot::make_twa_graph(dict);
    bdd a = bdd_ithvar(aut->register_ap("a"));
    aut->set_buchi();
    aut->new_states(2);
    aut->set_init_state(0);
    aut->new_edge(0, 0, a);
    aut->new_edge(0, 1, a);   // nondeterministic choice over 'a'
    aut->new_edge(1, 0, bddtrue);
    REQUIRE(!spot::is_deterministic(aut));

    spot::twa_graph_ptr det = kofola::determinize_tela(aut);
    REQUIRE(det != nullptr);
    CHECK(spot::is_deterministic(det));
    CHECK(det->is_empty());
}
