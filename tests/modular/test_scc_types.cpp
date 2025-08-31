#include <catch2/catch_test_macros.hpp>

#include <spot/parseaut/public.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/misc/bddlt.hh>

#include "kofola.hpp"
#include "../utils/test_utils.hpp"

TEST_CASE("cola::get_scc_types - scc_det.hoa", "[scc_types]") {
    // Load the test automaton
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/scc_det.hoa");
    REQUIRE(aut != nullptr);
    
    // Create SCC info
    spot::scc_info scc_info(aut);
    
    // Get SCC types
    std::string scc_types = cola::get_scc_types(scc_info);
    
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
    REQUIRE(cola::is_deterministic_scc(0, scc_info) == true);  // SCC 0 is deterministic
    REQUIRE(cola::is_deterministic_scc(1, scc_info) == true);  // SCC 1 is deterministic
    
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
    std::string scc_types = cola::get_scc_types(scc_info);
    
    // Test utility functions for SCC type checking
    REQUIRE(cola::is_accepting_scc(scc_types, 0) == true);   // SCC 0 is accepting
    REQUIRE(cola::is_accepting_scc(scc_types, 1) == false);  // SCC 1 is not accepting
    
    // SCC 0 is weak and accepting, so it's not a "deterministic" SCC in the sense of is_accepting_detscc
    // (which excludes weak SCCs)
    REQUIRE(cola::is_accepting_detscc(scc_types, 0) == false);  // SCC 0 is weak, not "det" in this context
    REQUIRE(cola::is_accepting_detscc(scc_types, 1) == false);  // SCC 1 is not accepting
    
    REQUIRE(cola::is_accepting_weakscc(scc_types, 0) == true);  // SCC 0 is weak and accepting
    REQUIRE(cola::is_accepting_weakscc(scc_types, 1) == false); // SCC 1 is not accepting
    
    REQUIRE(cola::is_weakscc(scc_types, 0) == true);   // SCC 0 is weak
    REQUIRE(cola::is_weakscc(scc_types, 1) == true);  // SCC 1 is inherently weak
    
    REQUIRE(cola::is_accepting_nondetscc(scc_types, 0) == false); // SCC 0 is deterministic and weak
    REQUIRE(cola::is_accepting_nondetscc(scc_types, 1) == false); // SCC 1 is not accepting
}

TEST_CASE("cola::get_scc_types - automaton properties", "[scc_types]") {
    // Load the test automaton
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("tests/test_data/scc_det.hoa");
    REQUIRE(aut != nullptr);
    
    // Create SCC info
    spot::scc_info scc_info(aut);
    std::string scc_types = cola::get_scc_types(scc_info);
    
    // Test higher-level automaton properties
    REQUIRE(cola::is_elevator_automaton(scc_info, scc_types) == true);  // Should be elevator (all SCCs det or weak)
    REQUIRE(cola::is_weak_automaton(scc_info, scc_types) == true);      // Should be weak (SCC 0 is weak)
    REQUIRE(cola::is_limit_deterministic_automaton(scc_info, scc_types) == true); // Should be limit deterministic
}
