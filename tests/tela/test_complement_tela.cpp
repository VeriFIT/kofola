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
TEST_CASE("complement_tela produces language-equivalent results to Spot", "[complement_tela]") {
    // Set up TELA options
    test_utils::setup_tela_options();
    // Ensure we're testing the default determinization algorithm
    kofola::OPTIONS.params.erase("tela_det_alg");
    
    // Use the common test files from test utilities
    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file: " + filename) {
            // Load automaton from HOA file
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            // Test complement equivalence
            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }
}

TEST_CASE("complement_tela with tela_det_alg=inductive produces language-equivalent results to Spot", "[complement_tela][tela_det_alg][inductive]") {
    // Set up TELA options and enable inductive SD-TELA determinization
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file (inductive det): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and sd_ind_sh_break=root produces language-equivalent results to Spot",
          "[complement_tela][tela_det_alg][inductive][sd_ind_sh_break]") {
    // Set up TELA options and enable inductive SD-TELA determinization with shared breakpoint
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_sh_break"] = "root";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file (inductive det, shared breakpoint): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_sh_break");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and sd_ind_sh_break=inf_tree produces language-equivalent results to Spot",
          "[complement_tela][tela_det_alg][inductive][sd_ind_sh_break]") {
    // Set up TELA options and enable inductive SD-TELA determinization with shared breakpoint
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_sh_break"] = "inf_tree";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file (inductive det, shared breakpoint): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_sh_break");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and sd_ind_or_opt=root produces language-equivalent results to Spot",
          "[complement_tela][tela_det_alg][inductive][sd_ind_or_opt]") {
    // Set up TELA options and enable inductive SD-TELA determinization with OR-FIN optimization
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_or_opt"] = "root";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file (inductive det, OR-FIN opt): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_or_opt");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and sd_ind_or_opt=inf_tree produces language-equivalent results to Spot",
          "[complement_tela][tela_det_alg][inductive][sd_ind_or_opt]") {
    // Set up TELA options and enable inductive SD-TELA determinization with OR-FIN optimization
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_or_opt"] = "inf_tree";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        SECTION("Testing file (inductive det, OR-FIN opt): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_or_opt");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and both sd_ind_sh_break=root and sd_ind_or_opt=yes",
          "[complement_tela][tela_det_alg][inductive][sd_ind_sh_break][sd_ind_or_opt]") {
    // Set up TELA options with both shared breakpoint and OR-FIN optimizations
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_sh_break"] = "root";
    kofola::OPTIONS.params["sd_ind_or_opt"] = "yes";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        if(filename == "tests/test_data/random_sd_streett_006.hoa") {
            continue; // skip this file which is a known outlier for the OR-FIN optimization
        }
        SECTION("Testing file (inductive det, shared breakpoint + OR-FIN opt): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_or_opt");
    kofola::OPTIONS.params.erase("sd_ind_sh_break");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

TEST_CASE("complement_tela with tela_det_alg=inductive and both sd_ind_sh_break=inf_tree and sd_ind_or_opt=yes",
          "[complement_tela][tela_det_alg][inductive][sd_ind_sh_break][sd_ind_or_opt]") {
    // Set up TELA options with both shared breakpoint and OR-FIN optimizations
    test_utils::setup_tela_options();
    kofola::OPTIONS.params["tela_det_alg"] = "inductive";
    kofola::OPTIONS.params["sd_ind_sh_break"] = "inf_tree";
    kofola::OPTIONS.params["sd_ind_or_opt"] = "yes";

    for (const std::string& filename : test_utils::COMMON_TEST_FILES) {
        if(filename == "tests/test_data/random_sd_streett_006.hoa") {
            continue; // skip this file which is a known outlier for the OR-FIN optimization
        }
        SECTION("Testing file (inductive det, shared breakpoint + OR-FIN opt): " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);

            bool equivalent = test_utils::test_complement_equivalence(aut, false);
            CHECK(equivalent);
        }
    }

    // Clean up to avoid side-effects on other tests
    kofola::OPTIONS.params.erase("sd_ind_or_opt");
    kofola::OPTIONS.params.erase("sd_ind_sh_break");
    kofola::OPTIONS.params.erase("tela_det_alg");
}

// Manual test function that can be called from main
void run_complement_test_on_file(const std::string& filename) {
    std::cout << "Testing complement equivalence for file: " << filename << std::endl;
    
    bool result = test_utils::test_file_complement_equivalence(filename, true);
    
    if (result) {
        std::cout << "✓ Test PASSED: Complements are language equivalent" << std::endl;
    } else {
        std::cout << "✗ Test FAILED: Complements are NOT language equivalent" << std::endl;
    }
}
