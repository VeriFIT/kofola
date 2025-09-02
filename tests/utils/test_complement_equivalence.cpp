// Standalone test program for comparing complement_tela with Spot's complement
//
// Usage: ./test_complement_equivalence <hoa_file>
//
// This program takes an automaton in HOA format as input, complements it using
// both complement_tela and Spot's complement, and checks if they produce
// language-equivalent results.

#include <iostream>
#include <fstream>
#include <string>

// Spot headers
#include <spot/parseaut/public.hh>
#include <spot/misc/bddlt.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/contains.hh>
#include <spot/twaalgos/hoa.hh>

// Kofola headers
#include "complement_tela.hpp"
#include "util/kofola.hpp"

/**
 * @brief Load an automaton from a HOA file.
 * 
 * @param filename Path to the HOA file
 * @return spot::twa_graph_ptr The loaded automaton, or nullptr if loading failed
 */
spot::twa_graph_ptr load_automaton_from_file(const std::string& filename) {
    try {
        spot::bdd_dict_ptr dict = spot::make_bdd_dict();
        spot::automaton_stream_parser parser(filename);
        spot::parsed_aut_ptr parsed_aut = parser.parse(dict);
        
        if (parsed_aut->format_errors(std::cerr)) {
            std::cerr << "Error parsing HOA file: " << filename << std::endl;
            return nullptr;
        }
        
        return parsed_aut->aut;
    } catch (const std::exception& ex) {
        std::cerr << "Exception loading automaton from " << filename << ": " << ex.what() << std::endl;
        return nullptr;
    }
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
        std::cout << "Computing complement using kofola::complement_tela..." << std::endl;
        
        // Complement using kofola::complement_tela
        spot::twa_graph_ptr kofola_complement = kofola::complement_tela(aut);
        if (!kofola_complement) {
            std::cerr << "kofola::complement_tela returned null" << std::endl;
            return false;
        }
        
        std::cout << "Kofola complement has " << kofola_complement->num_states() << " states" << std::endl;
        
        std::cout << "Computing complement using Spot..." << std::endl;
        
        // Complement using Spot's complement
        spot::twa_graph_ptr spot_complement = spot::complement(aut);
        if (!spot_complement) {
            std::cerr << "spot::complement returned null" << std::endl;
            return false;
        }
        
        std::cout << "Spot complement has " << spot_complement->num_states() << " states" << std::endl;
        
        std::cout << "Checking language equivalence..." << std::endl;
        
        // Check if they are language equivalent
        bool equivalent = spot::are_equivalent(kofola_complement, spot_complement);
        
        if (equivalent) {
            std::cout << "✓ SUCCESS: Complements are language equivalent" << std::endl;
        } else {
            std::cout << "✗ FAILURE: Complements are NOT language equivalent" << std::endl;
            
            // Try to find a distinguishing word
            std::cout << "Searching for distinguishing word..." << std::endl;
            spot::twa_word_ptr distinguishing_word = kofola_complement->exclusive_word(spot_complement);
            if (distinguishing_word) {
                std::cout << "Distinguishing word found (cannot be printed directly)" << std::endl;
                std::cout << "This word is accepted by exactly one of the complements" << std::endl;
            } else {
                std::cout << "No distinguishing word found (this shouldn't happen if not equivalent)" << std::endl;
            }
        }
        
        return equivalent;
        
    } catch (const std::exception& ex) {
        std::cerr << "Exception during complement comparison: " << ex.what() << std::endl;
        return false;
    }
}

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <hoa_file>" << std::endl;
    std::cout << std::endl;
    std::cout << "Test program for comparing complement_tela with Spot's complement." << std::endl;
    std::cout << "Takes an automaton in HOA format as input and checks if both" << std::endl;
    std::cout << "complement algorithms produce language-equivalent results." << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  hoa_file    Path to HOA format automaton file" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string filename = argv[1];
    
    std::cout << "=== Complement Equivalence Test ===" << std::endl;
    std::cout << "Input file: " << filename << std::endl;
    std::cout << std::endl;
    
    // Check if file exists
    std::ifstream file_check(filename);
    if (!file_check) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }
    file_check.close();
    
    // Load the automaton
    std::cout << "Loading automaton from HOA file..." << std::endl;
    spot::twa_graph_ptr aut = load_automaton_from_file(filename);
    if (!aut) {
        std::cerr << "Failed to load automaton from file: " << filename << std::endl;
        return 1;
    }
    
    // Print automaton info
    std::cout << "Loaded automaton successfully:" << std::endl;
    std::cout << "  States: " << aut->num_states() << std::endl;
    std::cout << "  Edges: " << aut->num_edges() << std::endl;
    std::cout << "  Acceptance condition: " << aut->get_acceptance() << std::endl;
    std::cout << std::endl;
    
    // Run the complement equivalence test
    bool equivalent = test_complement_equivalence(aut);
    
    std::cout << std::endl;
    std::cout << "=== Test Result ===" << std::endl;
    if (equivalent) {
        std::cout << "✓ PASSED: complement_tela produces language-equivalent results to Spot" << std::endl;
        return 0;
    } else {
        std::cout << "✗ FAILED: complement_tela does NOT produce language-equivalent results to Spot" << std::endl;
        return 1;
    }
}
