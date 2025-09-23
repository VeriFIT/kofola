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
    "tests/test_data/det_sccs_with_nondet_between.hoa",
    "tests/test_data/det_sccs_border.hoa",
    "tests/test_data/All_Sturmian_words_contain_cubes_sub.autfilt",
    "tests/test_data/All_Sturmian_words_contain_cubes_sup.autfilt",
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

    kofola::OPTIONS.params["merge_iwa"] = "yes";
    kofola::OPTIONS.params["merge_det"] = "yes";
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

// Specific regression test with simulation-based pruning enabled
TEST_CASE("modular complement with sim-ms-prune on ostrowski_thms-heur-197-autfilt.hoa", "[complement_modular][sim]") {
    // Ensure modular pipeline (no tela) and enable simulation-based pruning
    setup_modular_options();
    kofola::OPTIONS.params["sim-ms-prune"] = "yes";

    const std::string filename = "tests/test_data/ostrowski_thms-heur-197-autfilt.hoa";

    SECTION("Testing with sim-ms-prune=yes on ostrowski_thms-heur-197-autfilt.hoa") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
        REQUIRE(aut != nullptr);

        bool equivalent = test_utils::test_complement_equivalence(aut, false);
        CHECK(equivalent);
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sim-ms-prune");
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
