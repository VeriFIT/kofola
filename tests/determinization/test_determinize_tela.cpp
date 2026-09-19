// Tests for the modular determinization of TELA (kofola::determinize_tela).
//
// Only inherently weak partitions are supported so far, so the test data are
// weak (nondeterministic) automata.  For each of them we check that the result
// is deterministic and language equivalent to the input.

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
            CHECK(test_utils::test_determinize_equivalence(aut));
        }
    }
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
