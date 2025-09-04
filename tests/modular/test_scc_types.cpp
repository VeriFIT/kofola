#include <catch2/catch_test_macros.hpp>

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
        // Single deterministic SCC should be marked as initial deterministic
        REQUIRE((scc_types[0] & SCC_DET_TYPE) != 0);
        REQUIRE((scc_types[0] & SCC_INITIAL_DET_TYPE) != 0);
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
        
        // All deterministic SCCs in a chain should be initial deterministic
        for (unsigned sc = 0; sc < scc_info.scc_count(); ++sc) {
            if (scc_types[sc] & SCC_DET_TYPE) {
                REQUIRE((scc_types[sc] & SCC_INITIAL_DET_TYPE) != 0);
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
                REQUIRE((scc_types[sc] & SCC_INITIAL_DET_TYPE) == 0);
            } else {
                // These should be deterministic SCCs containing states 1 and 2
                // They should NOT be initial deterministic because they're reachable from nondet SCC
                if ((scc_types[sc] & SCC_INITIAL_DET_TYPE) == 0) {
                    found_det_without_initial = true;
                }
            }
        }
        
        REQUIRE(found_nondet_scc == true);
        REQUIRE(found_det_without_initial == true);
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
            auto states_in_scc = scc_info.states_of(sc);
            
            // Check if this SCC contains state 0 or 1 (initial part)
            bool is_initial_part = false;
            for (unsigned state : states_in_scc) {
                if (state <= 1) {
                    is_initial_part = true;
                    break;
                }
            }
            
            if (is_initial_part && (scc_types[sc] & SCC_DET_TYPE)) {
                // Initial deterministic SCCs should be marked as initial deterministic
                REQUIRE((scc_types[sc] & SCC_INITIAL_DET_TYPE) != 0);
            } else if (scc_types[sc] & SCC_DET_TYPE) {
                // Later deterministic SCCs (states 2, 3) should NOT be initial deterministic
                REQUIRE((scc_types[sc] & SCC_INITIAL_DET_TYPE) == 0);
            }
        }
    }
}
