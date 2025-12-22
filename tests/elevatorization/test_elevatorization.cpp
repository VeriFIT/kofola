// Test program for comparing complement_tela with Spot's complement
// 
// This program takes an automaton in HOA format as input, complements it using
// both complement_tela and Spot's complement, and checks if they produce
// language-equivalent results.

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <string>
#include <vector>

// Test utilities
#include "../utils/test_utils.hpp"

// Test case for complement_tela comparison
TEST_CASE("kofola::Elevatorization produces language-equivalent results", "[kofola::Elevatorization]") {
    // Set up TELA options
    test_utils::setup_tela_options();
    
    // Use the common test files from test utilities
    for (const std::string& filename : test_utils::ELEVATORIZATION_TEST_FILES) {
        SECTION("Testing file: " + filename) {
            // Load automaton from HOA file
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            // Test complement equivalence
            bool equivalent = test_utils::test_elevator_equivalence(aut, false);
            CHECK(equivalent);
        }
    }
}
