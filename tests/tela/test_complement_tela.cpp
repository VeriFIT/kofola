// Test program for comparing complement_tela with Spot's complement
// 
// This program takes an automaton in HOA format as input, complements it using
// both complement_tela and Spot's complement, and checks if they produce
// language-equivalent results.

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Spot headers
#include <spot/parseaut/public.hh>
#include <spot/misc/bddlt.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/contains.hh>
#include <spot/twaalgos/hoa.hh>

// Kofola headers
#include "complement_tela.hpp"
#include "kofola.hpp"

namespace {

/**
 * @brief Load an automaton from a HOA file.
 * 
 * @param filename Path to the HOA file
 * @return spot::twa_graph_ptr The loaded automaton, or nullptr if loading failed
 */
spot::twa_graph_ptr load_automaton_from_file(const std::string& filename) {
    // Try different possible paths depending on where tests are run from
    std::vector<std::string> possible_paths = {
        filename,                    // Direct path
        "../../" + filename,         // From build/tests directory
        "../" + filename            // From build directory
    };
    
    for (const std::string& path : possible_paths) {
        try {
            std::ifstream file(path);
            if (file.good()) {
                spot::bdd_dict_ptr dict = spot::make_bdd_dict();
                spot::automaton_stream_parser parser(path);
                spot::parsed_aut_ptr parsed_aut = parser.parse(dict);
                
                if (parsed_aut->format_errors(std::cerr)) {
                    std::cerr << "Error parsing HOA file: " << path << std::endl;
                    continue;
                }
                
                return parsed_aut->aut;
            }
        } catch (const std::exception& ex) {
            // Try next path
            continue;
        }
    }
    
    std::cerr << "Failed to load automaton from any of the tried paths for: " << filename << std::endl;
    return nullptr;
}

/**
 * @brief Compare complement_tela with Spot's complement for a given automaton.
 * 
 * @param aut The input automaton to complement
 * @return true if both complements are language equivalent, false otherwise
 */
bool test_complement_equivalence(const spot::twa_graph_ptr& aut) {
    if (!aut) {
        std::cerr << "Input automaton is null" << std::endl;
        return false;
    }
    
    try {
        // Complement using kofola::complement_tela
        spot::twa_graph_ptr kofola_complement = kofola::complement_tela(aut);
        if (!kofola_complement) {
            std::cerr << "kofola::complement_tela returned null" << std::endl;
            return false;
        }
        
        // Complement using Spot's complement
        spot::twa_graph_ptr spot_complement = spot::complement(aut);
        if (!spot_complement) {
            std::cerr << "spot::complement returned null" << std::endl;
            return false;
        }
        
        // Check if they are language equivalent
        bool equivalent = spot::are_equivalent(kofola_complement, spot_complement);
        
        if (equivalent) {
            std::cout << "✓ Complements are language equivalent" << std::endl;
        } else {
            std::cout << "✗ Complements are NOT language equivalent" << std::endl;
            
            // Try to find a distinguishing word
            spot::twa_word_ptr distinguishing_word = kofola_complement->exclusive_word(spot_complement);
            if (distinguishing_word) {
                std::cout << "Distinguishing word found (cannot be printed directly)" << std::endl;
            }
        }
        
        return equivalent;
        
    } catch (const std::exception& ex) {
        std::cerr << "Exception during complement comparison: " << ex.what() << std::endl;
        return false;
    }
}

} // anonymous namespace

// Test case for complement_tela comparison
TEST_CASE("complement_tela produces language-equivalent results to Spot", "[complement_tela]") {
    // Set the tela parameter to yes for TELA simplifications
    kofola::OPTIONS.params["tela"] = "yes";
    
    // Array of HOA test files to use for testing
    std::vector<std::string> test_files = {
        "tests/test_data/inf_a.hoa",
        "tests/test_data/ndet_example.hoa",
        "tests/test_data/simple_buchi.hoa",
        "tests/test_data/random_sd_streett_001.hoa",
        "tests/test_data/random_sd_streett_002.hoa",
        "tests/test_data/random_sd_streett_003.hoa",
        "tests/test_data/random_sd_streett_004.hoa",
        "tests/test_data/random_sd_streett_005.hoa",
        "tests/test_data/random_sd_streett_006.hoa",
    };
    
    for (const std::string& filename : test_files) {
        SECTION("Testing file: " + filename) {
            // Load automaton from HOA file
            spot::twa_graph_ptr aut = load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            // Test complement equivalence
            bool equivalent = test_complement_equivalence(aut);
            CHECK(equivalent);
        }
    }
}

// Manual test function that can be called from main
void run_complement_test_on_file(const std::string& filename) {
    std::cout << "Testing complement equivalence for file: " << filename << std::endl;
    
    spot::twa_graph_ptr aut = load_automaton_from_file(filename);
    if (!aut) {
        std::cout << "Failed to load automaton from file: " << filename << std::endl;
        return;
    }
    
    std::cout << "Loaded automaton with " << aut->num_states() << " states" << std::endl;
    std::cout << "Acceptance condition: " << aut->get_acceptance() << std::endl;
    
    bool equivalent = test_complement_equivalence(aut);
    
    if (equivalent) {
        std::cout << "✓ Test PASSED: Complements are language equivalent" << std::endl;
    } else {
        std::cout << "✗ Test FAILED: Complements are NOT language equivalent" << std::endl;
    }
}
