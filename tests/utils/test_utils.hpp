#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// Spot headers
#include <spot/parseaut/public.hh>
#include <spot/misc/bddlt.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/contains.hh>
#include <spot/twaalgos/hoa.hh>

// Kofola headers
#include "complement_tela.hpp"
#include "kofola.hpp"

namespace test_utils {

/**
 * @brief Common test data files used across multiple test suites
 */
const std::vector<std::string> COMMON_TEST_FILES = {
    "tests/test_data/inf_a.hoa",
    "tests/test_data/ndet_example.hoa",
    "tests/test_data/simple_buchi.hoa",
    "tests/test_data/random_sd_streett_001.hoa",
    "tests/test_data/random_sd_streett_002.hoa",
    "tests/test_data/random_sd_streett_003.hoa",
    "tests/test_data/random_sd_streett_004.hoa",
    "tests/test_data/random_sd_streett_005.hoa",
    "tests/test_data/random_sd_streett_006.hoa",
    // TELA determinization automata
    "tests/test_data/tela_det/automaton0.hoa",
    "tests/test_data/tela_det/automaton1.hoa",
    "tests/test_data/tela_det/automaton102.hoa",
    "tests/test_data/tela_det/automaton103.hoa",
    "tests/test_data/tela_det/automaton108.hoa",
    "tests/test_data/tela_det/automaton130.hoa",
    "tests/test_data/tela_det/automaton162.hoa",
    // "tests/test_data/tela_det/automaton188.hoa",
    "tests/test_data/tela_det/automaton204.hoa",
    "tests/test_data/tela_det/automaton208.hoa",
    "tests/test_data/tela_det/automaton216.hoa",
    "tests/test_data/tela_det/automaton219.hoa",
    "tests/test_data/tela_det/automaton240.hoa",
    "tests/test_data/tela_det/automaton248.hoa",
    "tests/test_data/tela_det/automaton281.hoa",
    "tests/test_data/tela_det/automaton288.hoa",
    "tests/test_data/tela_det/automaton297.hoa",
    "tests/test_data/tela_det/automaton327.hoa",
    "tests/test_data/tela_det/automaton380.hoa",
    "tests/test_data/tela_det/automaton4.hoa",
    "tests/test_data/tela_det/automaton400.hoa",
    "tests/test_data/tela_det/automaton426.hoa",
    "tests/test_data/tela_det/automaton427.hoa",
    "tests/test_data/tela_det/automaton430.hoa",
    "tests/test_data/tela_det/automaton432.hoa",
    "tests/test_data/tela_det/automaton440.hoa",
    "tests/test_data/tela_det/automaton461.hoa",
    "tests/test_data/tela_det/automaton471.hoa",
    "tests/test_data/tela_det/automaton473.hoa",
    "tests/test_data/tela_det/automaton48.hoa",
    "tests/test_data/tela_det/automaton517.hoa",
    "tests/test_data/tela_det/automaton521.hoa",
    "tests/test_data/tela_det/automaton522.hoa",
    "tests/test_data/tela_det/automaton523.hoa",
    "tests/test_data/tela_det/automaton564.hoa",
    "tests/test_data/tela_det/automaton60.hoa",
    "tests/test_data/tela_det/automaton603.hoa",
    "tests/test_data/tela_det/automaton611.hoa",
    "tests/test_data/tela_det/automaton618.hoa",
    "tests/test_data/tela_det/automaton640.hoa",
    "tests/test_data/tela_det/automaton641.hoa",
    "tests/test_data/tela_det/automaton658.hoa",
    "tests/test_data/tela_det/automaton660.hoa",
    "tests/test_data/tela_det/automaton68.hoa",
    "tests/test_data/tela_det/automaton682.hoa",
    "tests/test_data/tela_det/automaton690.hoa",
    "tests/test_data/tela_det/automaton705.hoa",
    "tests/test_data/tela_det/automaton740.hoa",
    "tests/test_data/tela_det/automaton752.hoa",
    "tests/test_data/tela_det/automaton800.hoa",
    "tests/test_data/tela_det/automaton807.hoa",
    "tests/test_data/tela_det/automaton825.hoa",
    "tests/test_data/tela_det/automaton826.hoa",
    "tests/test_data/tela_det/automaton850.hoa",
    "tests/test_data/tela_det/automaton86.hoa",
    "tests/test_data/tela_det/automaton860.hoa",
    "tests/test_data/tela_det/automaton90.hoa",
    "tests/test_data/tela_det/automaton900.hoa",
    "tests/test_data/tela_det/automaton940.hoa"
};

/**
 * @brief Extended test data files for comprehensive testing
 */
const std::vector<std::string> EXTENDED_TEST_FILES = {
    "tests/test_data/inf_a.hoa",
    "tests/test_data/ndet_example.hoa",
    "tests/test_data/simple_buchi.hoa",
    "tests/test_data/multi_scc_buchi.hoa",
    "tests/test_data/mixed_scc_types.hoa",
    "tests/test_data/random_sd_streett_001.hoa",
    "tests/test_data/random_sd_streett_002.hoa",
    "tests/test_data/random_sd_streett_003.hoa",
    "tests/test_data/random_sd_streett_004.hoa",
    "tests/test_data/random_sd_streett_005.hoa",
    "tests/test_data/random_sd_streett_006.hoa"
};

/**
 * @brief Load an automaton from a HOA file, trying multiple possible paths.
 * 
 * This function handles the common case where tests may be run from different
 * directories (build/, build/tests/, etc.) by trying multiple relative paths.
 * 
 * @param filename Path to the HOA file
 * @return spot::twa_graph_ptr The loaded automaton, or nullptr if loading failed
 */
spot::twa_graph_ptr load_automaton_from_file(const std::string& filename);

/**
 * @brief Load an automaton from a HOA file using exact path (no path resolution).
 * 
 * This is a simpler version that doesn't try multiple paths, useful when the
 * exact path is known.
 * 
 * @param filename Exact path to the HOA file
 * @return spot::twa_graph_ptr The loaded automaton, or nullptr if loading failed
 */
spot::twa_graph_ptr load_automaton_exact_path(const std::string& filename);

/**
 * @brief Compare complement_tela with Spot's complement for a given automaton.
 * 
 * @param aut The input automaton to complement
 * @param verbose Whether to print detailed output during comparison
 * @return true if both complements are language equivalent, false otherwise
 */
bool test_complement_equivalence(const spot::twa_graph_ptr& aut, bool verbose = false);

/**
 * @brief Test complement equivalence for a file with comprehensive output.
 * 
 * This is a higher-level utility that loads an automaton from a file and
 * performs the complement equivalence test with detailed logging.
 * 
 * @param filename Path to the HOA file
 * @param verbose Whether to print detailed output
 * @return true if the test passes, false otherwise
 */
bool test_file_complement_equivalence(const std::string& filename, bool verbose = true);

/**
 * @brief Set up common kofola options for TELA complementation.
 * 
 * This function sets standard options that are commonly used across tests.
 */
void setup_tela_options();

/**
 * @brief Print automaton basic information.
 * 
 * @param aut The automaton to print info about
 * @param label Optional label to prefix the output
 */
void print_automaton_info(const spot::twa_graph_ptr& aut, const std::string& label = "");

} // namespace test_utils
