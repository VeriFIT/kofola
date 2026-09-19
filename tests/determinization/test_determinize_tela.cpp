// Tests for the modular determinization of TELA (kofola::determinize_tela).
//
// Inherently weak and deterministic partitions are supported, which together
// cover every elevator automaton; a nondeterministic accepting SCC is the one
// case left unimplemented.  The first test case runs over a handful of weak
// (nondeterministic) automata, the second one over every elevator automaton of
// test_data.  In both cases we check that the result is deterministic and
// language equivalent to the input.

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

        // a nondeterministic accepting SCC makes the automaton non-elevator;
        // there is no partial determinization algorithm for it yet
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
