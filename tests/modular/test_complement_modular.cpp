// Test program for comparing complement_tela (with default modular settings) with Spot's complement
// 
// This program takes an automaton in HOA format as input, complements it using
// both complement_tela (with default modular settings) and Spot's complement, 
// and checks if they produce language-equivalent results.

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <string>
#include <vector>

// Test utilities
#include "../utils/test_utils.hpp"

// Test files specific to modular complementation testing
const std::vector<std::string> MODULAR_TEST_FILES = {
    "tests/test_data/simple_buchi.hoa",
};

/**
 * @brief Set up default modular options for kofola.
 * 
 * This function sets up kofola with default modular complementation options
 * (i.e., does NOT set tela="yes", allowing the modular approach to be used).
 */
void setup_modular_options() {
    // Use default options - specifically do NOT set tela="yes"
    // This allows kofola to use its modular complementation approach
    // instead of the TELA-specific optimizations
    
    // Clear any existing tela setting to ensure default modular behavior
    kofola::OPTIONS.params.erase("tela");
}

// Test case for modular complement comparison
TEST_CASE("complement_tela with default modular settings produces language-equivalent results to Spot", "[complement_modular]") {
    // Set up modular options (default settings)
    setup_modular_options();
    
    // Use the modular-specific test files
    for (const std::string& filename : MODULAR_TEST_FILES) {
        SECTION("Testing modular complement for file: " + filename) {
            // Load automaton from HOA file
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            // Test complement equivalence
            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }
}

// Manual test function that can be called from main
void run_modular_complement_test_on_file(const std::string& filename) {
    std::cout << "Testing modular complement equivalence for file: " << filename << std::endl;
    
    // Set up modular options
    setup_modular_options();
    
    bool result = test_utils::test_file_complement_equivalence(filename, true);
    
    if (result) {
        std::cout << "✓ Test PASSED: Modular complements are language equivalent" << std::endl;
    } else {
        std::cout << "✗ Test FAILED: Modular complements are NOT language equivalent" << std::endl;
    }
}
