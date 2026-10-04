#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <sstream>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/misc/bddlt.hh>

#include "util/helpers.hpp"
#include "../utils/test_utils.hpp"

TEST_CASE("cola::get_scc_types - scc_det.hoa", "[scc_types]") {
    // Load the test automaton
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/scc_det.hoa");
    REQUIRE(aut != nullptr);
    
    // Create SCC info
    spot::scc_info scc_info(aut);
    
    // Get SCC types
    std::string scc_types = helpers::get_scc_types(scc_info);
    
    // Verify the automaton has 2 SCCs
    REQUIRE(scc_info.scc_count() == 2);
    REQUIRE(scc_types.size() == 2);
    
    // Analyze the automaton structure:
    // State 0: transitions [0] 0, [0] 1 (self-loop and to state 1)
    // State 1: transition [0] 1 {0} (self-loop with acceptance mark)
    // This creates 2 SCCs:
    // - SCC 0: contains only state 1 (accepting, deterministic, weak due to single acceptance set)
    // - SCC 1: contains only state 0 (non-accepting, deterministic)
    
    INFO("SCC 0 (state 1): " << static_cast<int>(scc_types[0]));
    INFO("SCC 1 (state 0): " << static_cast<int>(scc_types[1]));
    
    // SCC 0 (containing state 1) should be:
    // - Deterministic (both inside and fully deterministic)
    // - Accepting (has acceptance mark {0})
    // - Weak (inherently weak in Spot's analysis)
    REQUIRE((scc_types[0] & SCC_INSIDE_DET_TYPE) != 0); // Inside deterministic
    REQUIRE((scc_types[0] & SCC_DET_TYPE) != 0);        // Fully deterministic
    REQUIRE((scc_types[0] & SCC_ACC) != 0);             // Accepting
    REQUIRE((scc_types[0] & SCC_WEAK_TYPE) != 0);       // Inherently weak (single accepting state)
    
    // SCC 1 (containing state 0) should be:
    // - Deterministic (both inside and fully deterministic)
    // - Non-accepting (no acceptance marks)
    // - Not weak (no acceptance, cannot be weak)
    REQUIRE((scc_types[1] & SCC_INSIDE_DET_TYPE) != 0); // Inside deterministic
    REQUIRE((scc_types[1] & SCC_DET_TYPE) == 0);        // Fully deterministic
    REQUIRE((scc_types[1] & SCC_ACC) == 0);             // Not accepting
    REQUIRE((scc_types[1] & SCC_WEAK_TYPE) != 0);       // Not inherently weak
    
    // Verify individual SCC determinism using cola::is_deterministic_scc
    REQUIRE(helpers::is_deterministic_scc(0, scc_info) == true);  // SCC 0 is deterministic
    REQUIRE(helpers::is_deterministic_scc(1, scc_info) == true);  // SCC 1 is deterministic
    
    // Verify SCC acceptance
    REQUIRE(scc_info.is_accepting_scc(0) == true);   // SCC 0 (state 1) is accepting
    REQUIRE(scc_info.is_accepting_scc(1) == false);  // SCC 1 (state 0) is not accepting
}

TEST_CASE("cola::get_scc_types - utility functions", "[scc_types]") {
    // Load the test automaton
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/scc_det.hoa");
    REQUIRE(aut != nullptr);
    
    // Create SCC info
    spot::scc_info scc_info(aut);
    std::string scc_types = helpers::get_scc_types(scc_info);
    
    // Test utility functions for SCC type checking
    REQUIRE(helpers::is_accepting_scc(scc_types, 0) == true);   // SCC 0 is accepting
    REQUIRE(helpers::is_accepting_scc(scc_types, 1) == false);  // SCC 1 is not accepting
    
    // SCC 0 is weak and accepting, so it's not a "deterministic" SCC in the sense of is_accepting_detscc
    // (which excludes weak SCCs)
    REQUIRE(helpers::is_accepting_detscc(scc_types, 0) == false);  // SCC 0 is weak, not "det" in this context
    REQUIRE(helpers::is_accepting_detscc(scc_types, 1) == false);  // SCC 1 is not accepting
    
    REQUIRE(helpers::is_accepting_weakscc(scc_types, 0) == true);  // SCC 0 is weak and accepting
    REQUIRE(helpers::is_accepting_weakscc(scc_types, 1) == false); // SCC 1 is not accepting
    
    REQUIRE(helpers::is_weakscc(scc_types, 0) == true);   // SCC 0 is weak
    REQUIRE(helpers::is_weakscc(scc_types, 1) == true);  // SCC 1 is inherently weak
    
    REQUIRE(helpers::is_accepting_nondetscc(scc_types, 0) == false); // SCC 0 is deterministic and weak
    REQUIRE(helpers::is_accepting_nondetscc(scc_types, 1) == false); // SCC 1 is not accepting
}

TEST_CASE("cola::get_scc_types - automaton properties", "[scc_types]") {
    // Load the test automaton
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/scc_det.hoa");
    REQUIRE(aut != nullptr);
    
    // Create SCC info
    spot::scc_info scc_info(aut);
    std::string scc_types = helpers::get_scc_types(scc_info);
    
    // Test higher-level automaton properties
    REQUIRE(helpers::is_elevator_automaton(scc_info, scc_types) == true);  // Should be elevator (all SCCs det or weak)
    REQUIRE(helpers::is_weak_automaton(scc_info, scc_types) == true);      // Should be weak (SCC 0 is weak)
    REQUIRE(helpers::is_limit_deterministic_automaton(scc_info, scc_types) == true); // Should be limit deterministic
}

TEST_CASE("cola::get_scc_types - initial deterministic components", "[scc_types]") {
    SECTION("Test with simple HOA string - single deterministic SCC") {
        // Create a simple automaton with one deterministic SCC using HOA format
        std::string hoa_str = 
            "HOA: v1\n"
            "States: 1\n"
            "Start: 0\n"
            "AP: 1 \"a\"\n"
            "acc-name: Buchi\n"
            "Acceptance: 1 Inf(0)\n"
            "--BODY--\n"
            "State: 0\n"
            "[0] 0 {0}\n"
            "[!0] 0\n"
            "--END--\n";
        
        auto dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
        auto parsed_aut = parser.parse(dict);
        REQUIRE(!parsed_aut->format_errors(std::cerr));
        
        spot::twa_graph_ptr aut = parsed_aut->aut;
        REQUIRE(aut != nullptr);
        
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);
        
        REQUIRE(scc_info.scc_count() == 1);
        // Single deterministic SCC should be marked as almost initial deterministic
        REQUIRE((scc_types[0] & SCC_DET_TYPE) != 0);
        REQUIRE((scc_types[0] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) != 0);
    }
    
    SECTION("Test with HOA string - chain of deterministic SCCs") {
        // Create automaton: state 0 -> state 1 -> state 2 (all deterministic)
        std::string hoa_str = 
            "HOA: v1\n"
            "States: 3\n"
            "Start: 0\n"
            "AP: 1 \"a\"\n"
            "acc-name: Buchi\n"
            "Acceptance: 1 Inf(0)\n"
            "--BODY--\n"
            "State: 0\n"
            "[0] 1\n"
            "State: 1\n"
            "[0] 2\n"
            "State: 2\n"
            "[0] 2 {0}\n"
            "--END--\n";
        
        auto dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
        auto parsed_aut = parser.parse(dict);
        REQUIRE(!parsed_aut->format_errors(std::cerr));
        
        spot::twa_graph_ptr aut = parsed_aut->aut;
        REQUIRE(aut != nullptr);
        
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);
        
        // All deterministic SCCs in a chain should be almost initial deterministic
        for (unsigned sc = 0; sc < scc_info.scc_count(); ++sc) {
            if (scc_types[sc] & SCC_DET_TYPE) {
                REQUIRE((scc_types[sc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) != 0);
            }
        }
    }
    
    SECTION("Test nondeterministic choice leading to deterministic SCCs") {
        // Nondeterministic choice from state 0: 0 -> 1 and 0 -> 2, then deterministic loops
        std::string hoa_str = 
            "HOA: v1\n"
            "States: 3\n"
            "Start: 0\n"
            "AP: 1 \"a\"\n"
            "acc-name: Buchi\n"
            "Acceptance: 1 Inf(0)\n"
            "--BODY--\n"
            "State: 0\n"
            "[0] 1\n"
            "[0] 2\n"
            "State: 1\n"
            "[0] 1 {0}\n"
            "State: 2\n"
            "[0] 2 {0}\n"
            "--END--\n";
        
        auto dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
        auto parsed_aut = parser.parse(dict);
        REQUIRE(!parsed_aut->format_errors(std::cerr));
        
        spot::twa_graph_ptr aut = parsed_aut->aut;
        REQUIRE(aut != nullptr);
        
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);
        
        // Find the nondeterministic SCC (should contain state 0)
        // and deterministic SCCs (should contain states 1 and 2)
        bool found_nondet_scc = false;
        bool found_det_without_initial = false;
        
        for (unsigned sc = 0; sc < scc_info.scc_count(); ++sc) {
            auto states_in_scc = scc_info.states_of(sc);
            
            if (!(scc_types[sc] & SCC_DET_TYPE)) {
                // This should be the nondeterministic SCC containing state 0
                found_nondet_scc = true;
                REQUIRE((scc_types[sc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) != 0);
            } else {
                // These should be deterministic SCCs containing states 1 and 2
                // They should NOT be initial almost deterministic because they're reachable from nondet SCC
                if ((scc_types[sc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) == 0) {
                    found_det_without_initial = true;
                }
            }
        }
        
        REQUIRE(found_nondet_scc == true);
        REQUIRE(found_det_without_initial == false);
    }
    
    SECTION("Test mixed scenario - initial det, then nondet, then det") {
        // Initial deterministic: 0 -> 1, then nondet: 1 -> 2,3, then det: 2->2, 3->3
        std::string hoa_str = 
            "HOA: v1\n"
            "States: 4\n"
            "Start: 0\n"
            "AP: 2 \"a\" \"b\"\n"
            "acc-name: Buchi\n"
            "Acceptance: 1 Inf(0)\n"
            "--BODY--\n"
            "State: 0\n"
            "[0] 1\n"
            "State: 1\n"
            "[1] 2\n"
            "[1] 3\n"
            "State: 2\n"
            "[0] 2 {0}\n"
            "State: 3\n"
            "[0] 3 {0}\n"
            "--END--\n";
        
        auto dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
        auto parsed_aut = parser.parse(dict);
        REQUIRE(!parsed_aut->format_errors(std::cerr));
        
        spot::twa_graph_ptr aut = parsed_aut->aut;
        REQUIRE(aut != nullptr);
        
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);
        
        // Check that only initial deterministic SCCs are marked as such
        for (unsigned sc = 0; sc < scc_info.scc_count(); ++sc) {
            REQUIRE((scc_types[sc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) != 0);
        }
    }
}

TEST_CASE("kofola::get_scc_types - generalized initial deterministic", "[scc_types]") {

    SECTION("Nondet Border") {
        // Load the test automaton
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/det_sccs_with_nondet_between.hoa");
        REQUIRE(aut != nullptr);
            
        // Create SCC info
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);
     
        // Test utility functions for SCC type checking
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 0) == false);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 1) == false);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 2) == false);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 3) == false);
    }
    
    SECTION("Det Border") {
        // Load the test automaton
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/det_sccs_border.hoa");
        REQUIRE(aut != nullptr);
            
        // Create SCC info
        spot::scc_info scc_info(aut);
        std::string scc_types = helpers::get_scc_types(scc_info);

        // Test utility functions for SCC type checking
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 0) == true);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 1) == true);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 2) == false);
        REQUIRE(helpers::is_accepting_initial_almost_detscc(scc_types, 3) == false);
    }
}

TEST_CASE("cola::is_emerson_lei_elevator_automaton - SCC not using all sets", "[scc_types]") {
    // One nondeterministic SCC (state 0 has two "t" edges) that is not
    // inherently weak but whose acceptance, restricted to the SCC, is
    // generalized co-Buchi.  The first variant declares only the set the
    // SCC uses; the second declares an extra set the SCC never sees, so
    // the restricted condition is Fin(1) over 2 declared sets.
    auto acceptance = GENERATE(
        std::make_pair(std::string("1 Fin(0)"), std::string("{0}")),
        std::make_pair(std::string("2 Inf(0) | Fin(1)"), std::string("{1}")));

    std::string hoa_str =
        "HOA: v1\n"
        "States: 2\n"
        "Start: 0\n"
        "AP: 0\n"
        "Acceptance: " + acceptance.first + "\n"
        "--BODY--\n"
        "State: 0\n"
        "[t] 0\n"
        "[t] 1 " + acceptance.second + "\n"
        "State: 1\n"
        "[t] 0\n"
        "--END--\n";
    INFO("Acceptance: " << acceptance.first);

    auto dict = spot::make_bdd_dict();
    spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
    auto parsed_aut = parser.parse(dict);
    REQUIRE(!parsed_aut->format_errors(std::cerr));
    spot::twa_graph_ptr aut = parsed_aut->aut;
    REQUIRE(aut != nullptr);

    REQUIRE(helpers::is_elevator_automaton(aut) == false);
    REQUIRE(helpers::is_emerson_lei_elevator_automaton(aut) == true);

    spot::scc_info scc_info(aut);
    std::string scc_types = helpers::get_scc_types(scc_info);
    REQUIRE(helpers::is_elevator_automaton(scc_info, scc_types) == false);
    REQUIRE(helpers::is_emerson_lei_elevator_automaton(scc_info, scc_types) == true);
}

namespace {
    spot::twa_graph_ptr parse_hoa(const std::string& hoa_str) {
        auto dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(hoa_str.c_str(), "test_string");
        auto parsed_aut = parser.parse(dict);
        REQUIRE(!parsed_aut->format_errors(std::cerr));
        REQUIRE(parsed_aut->aut != nullptr);
        return parsed_aut->aut;
    }
}

TEST_CASE("kofola::get_scc_types - alternating automaton with Fin", "[scc_types]") {
    // Universal edge 0 -> 1&2.  Spot's determine_unknown_acceptance()
    // throws on alternating automata, and scc_of() on a universal
    // destination reads out of bounds.
    spot::twa_graph_ptr aut = parse_hoa(
        "HOA: v1\n"
        "States: 3\n"
        "Start: 0\n"
        "AP: 1 \"a\"\n"
        "Acceptance: 1 Fin(0)\n"
        "--BODY--\n"
        "State: 0\n"
        "[0] 1&2 {0}\n"
        "[!0] 1\n"
        "State: 1\n"
        "[t] 0\n"
        "State: 2\n"
        "[t] 2\n"
        "--END--\n");
    REQUIRE(!aut->is_existential());

    spot::scc_info scc_info(aut);
    std::string scc_types;
    REQUIRE_NOTHROW(scc_types = helpers::get_scc_types(scc_info));
    REQUIRE(scc_types.size() == scc_info.scc_count());

    REQUIRE(helpers::is_elevator_automaton(aut) == false);
    REQUIRE(helpers::is_emerson_lei_elevator_automaton(aut) == false);
    REQUIRE(helpers::is_elevator_automaton(scc_info, scc_types) == false);
    REQUIRE_NOTHROW(helpers::is_weak_automaton(aut));
}

TEST_CASE("kofola::get_scc_types - acceptance left untouched", "[scc_types]") {
    // Colour 2 is declared but not used by the formula.  Resolving the
    // unknown acceptance (and the weak-SCC check) makes Spot rebuild the
    // condition from the formula alone, shrinking num_sets() to 2.
    spot::twa_graph_ptr aut = parse_hoa(
        "HOA: v1\n"
        "States: 2\n"
        "Start: 0\n"
        "AP: 1 \"a\"\n"
        "Acceptance: 3 Fin(0) & Inf(1)\n"
        "--BODY--\n"
        "State: 0\n"
        "[0] 1 {0 2}\n"
        "[0] 1 {1}\n"
        "State: 1\n"
        "[t] 0\n"
        "--END--\n");
    const spot::acc_cond orig_acc = aut->acc();
    REQUIRE(orig_acc.num_sets() == 3);

    helpers::is_elevator_automaton(aut);
    CHECK(aut->acc() == orig_acc);
    helpers::is_emerson_lei_elevator_automaton(aut);
    CHECK(aut->acc() == orig_acc);
    helpers::is_weak_automaton(aut);
    CHECK(aut->acc() == orig_acc);
    spot::scc_info scc_info(aut);
    helpers::get_scc_types(scc_info);
    CHECK(aut->acc() == orig_acc);
}

TEST_CASE("kofola::get_scc_types - accepting SCC under mixed Fin/Inf", "[scc_types]") {
    // The cycle through the !a edge sees only colour 1, so the SCC is
    // accepting, but scc_info cannot tell without an emptiness check.
    spot::twa_graph_ptr aut = parse_hoa(
        "HOA: v1\n"
        "States: 2\n"
        "Start: 0\n"
        "AP: 1 \"a\"\n"
        "Acceptance: 2 Fin(0) & Inf(1)\n"
        "--BODY--\n"
        "State: 0\n"
        "[0] 1 {0}\n"
        "[!0] 1 {1}\n"
        "State: 1\n"
        "[t] 0\n"
        "--END--\n");

    spot::scc_info scc_info(aut);
    REQUIRE(scc_info.scc_count() == 1);
    std::string scc_types = helpers::get_scc_types(scc_info);
    REQUIRE((scc_types[0] & SCC_ACC) != 0);
    REQUIRE((scc_types[0] & SCC_WEAK_TYPE) == 0);
}

TEST_CASE("kofola::is_weak_automaton - rejecting SCC under mixed Fin/Inf", "[scc_types]") {
    // Every cycle is rejecting (the one through state 1 sees colour 0,
    // the self-loop never sees colour 1), so the SCC is inherently weak.
    spot::twa_graph_ptr aut = parse_hoa(
        "HOA: v1\n"
        "States: 2\n"
        "Start: 0\n"
        "AP: 1 \"a\"\n"
        "Acceptance: 2 Fin(0) & Inf(1)\n"
        "--BODY--\n"
        "State: 0\n"
        "[0] 0\n"
        "[!0] 1 {1}\n"
        "State: 1\n"
        "[t] 0 {0}\n"
        "--END--\n");

    REQUIRE(helpers::is_weak_automaton(aut) == true);
    REQUIRE(helpers::is_elevator_automaton(aut) == true);

    spot::scc_info scc_info(aut);
    std::string scc_types = helpers::get_scc_types(scc_info);
    REQUIRE((scc_types[0] & SCC_WEAK_TYPE) != 0);
    REQUIRE((scc_types[0] & SCC_ACC) == 0);
    REQUIRE(helpers::is_weak_automaton(scc_info, scc_types) == true);
}
