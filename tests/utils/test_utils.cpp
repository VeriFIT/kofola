#include "test_utils.hpp"

namespace test_utils {

spot::twa_graph_ptr load_automaton_from_file(const std::string& filename) {
    // Try different possible paths depending on where tests are run from
    std::vector<std::string> possible_paths = {
        filename,                    // Direct path
        "../../" + filename,         // From build/tests directory
        "../" + filename,            // From build directory
        "../tests/" + filename,      // From build to tests
        "tests/" + filename          // From root directory
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

spot::twa_graph_ptr load_automaton_exact_path(const std::string& filename) {
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

bool test_elevator_equivalence(const spot::twa_graph_ptr& aut, bool verbose) {
    if (!aut) {
        std::cerr << "Input automaton is null" << std::endl;
        return false;
    }
    
    try {
        if (verbose) {
            std::cout << "Transforming to elevator using kofola::Elevatorization..." << std::endl;
        }
        
        kofola::Elevatorization elev(aut);
        spot::twa_graph_ptr elevatorized = elev.elevatorize();
        if (!elevatorized) {
            std::cerr << "kofola::Elevatorization::elevarize() returned null" << std::endl;
            return false;
        }
        
        if (verbose) {
            std::cout << "Kofola complement has " << elevatorized->num_states() << " states" << std::endl;
        }
        
        // Check if they are language equivalent
        bool equivalent = spot::are_equivalent(elevatorized, aut);
        
        if (equivalent) {
            if (verbose) {
                std::cout << "✓ SUCCESS: Transformation to elevator preserved the language" << std::endl;
            }
        } else {
            std::cout << "✗ FAILURE: Transformation to elevator did NOT preserve the language" << std::endl;
            
            if (verbose) {
                // Try to find a distinguishing word
                std::cout << "Searching for distinguishing word..." << std::endl;
                spot::twa_word_ptr distinguishing_word = elevatorized->exclusive_word(aut);
                if (distinguishing_word) {
                    std::cout << "Distinguishing word found (cannot be printed directly)" << std::endl;
                    std::cout << "This word is accepted by exactly one of the automata" << std::endl;
                } else {
                    std::cout << "No distinguishing word found (this shouldn't happen if not equivalent)" << std::endl;
                }
            }
        }
        
        return equivalent;
        
    } catch (const std::exception& ex) {
        std::cerr << "Exception during elevatorization language preservation check: " << ex.what() << std::endl;
        return false;
    }
}

bool test_complement_equivalence(const spot::twa_graph_ptr& aut, bool verbose) {
    if (!aut) {
        std::cerr << "Input automaton is null" << std::endl;
        return false;
    }
    
    try {
        if (verbose) {
            std::cout << "Computing complement using kofola::complement_tela..." << std::endl;
        }
        
        // Complement using kofola::complement_tela
        spot::twa_graph_ptr kofola_complement = kofola::complement_tela(aut);
        if (!kofola_complement) {
            std::cerr << "kofola::complement_tela returned null" << std::endl;
            return false;
        }
        
        if (verbose) {
            std::cout << "Kofola complement has " << kofola_complement->num_states() << " states" << std::endl;
            std::cout << "Computing complement using Spot..." << std::endl;
        }
        
        // Complement using Spot's complement
        spot::twa_graph_ptr spot_complement = spot::complement(aut);
        if (!spot_complement) {
            std::cerr << "spot::complement returned null" << std::endl;
            return false;
        }
        
        if (verbose) {
            std::cout << "Spot complement has " << spot_complement->num_states() << " states" << std::endl;
            std::cout << "Checking language equivalence..." << std::endl;
        }
        
        // Check if they are language equivalent
        bool equivalent = spot::are_equivalent(kofola_complement, spot_complement);
        
        if (equivalent) {
            if (verbose) {
                std::cout << "✓ SUCCESS: Complements are language equivalent" << std::endl;
            }
        } else {
            std::cout << "✗ FAILURE: Complements are NOT language equivalent" << std::endl;
            
            if (verbose) {
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
        }
        
        return equivalent;
        
    } catch (const std::exception& ex) {
        std::cerr << "Exception during complement comparison: " << ex.what() << std::endl;
        return false;
    }
}

bool test_file_complement_equivalence(const std::string& filename, bool verbose) {
    if (verbose) {
        std::cout << "Testing complement equivalence for file: " << filename << std::endl;
    }
    
    spot::twa_graph_ptr aut = load_automaton_from_file(filename);
    if (!aut) {
        std::cout << "Failed to load automaton from file: " << filename << std::endl;
        return false;
    }
    
    if (verbose) {
        print_automaton_info(aut, "Loaded automaton");
    }
    
    bool equivalent = test_complement_equivalence(aut, verbose);
    
    if (verbose) {
        if (equivalent) {
            std::cout << "✓ Test PASSED: Complements are language equivalent" << std::endl;
        } else {
            std::cout << "✗ Test FAILED: Complements are NOT language equivalent" << std::endl;
        }
    }
    
    return equivalent;
}

void setup_tela_options() {
    // Set the tela parameter to yes for TELA simplifications
    kofola::OPTIONS.params["tela"] = "yes";
}

void print_automaton_info(const spot::twa_graph_ptr& aut, const std::string& label) {
    if (!aut) {
        std::cout << label << ": null automaton" << std::endl;
        return;
    }
    
    std::string prefix = label.empty() ? "Automaton" : label;
    std::cout << prefix << " has " << aut->num_states() << " states" << std::endl;
    std::cout << prefix << " acceptance condition: " << aut->get_acceptance() << std::endl;
}

} // namespace test_utils
