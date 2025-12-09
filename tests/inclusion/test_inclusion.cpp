#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <functional>

// Spot
#include <spot/twa/twa.hh>

// Kofola
#include "inclusion/inclusion_check.hpp"

// Test utils
#include "../utils/test_utils.hpp"

namespace {
// Simple harness that checks A ⊆ B using kofola::inclusion_check
bool check_inclusion_kofola(const spot::twa_graph_ptr& A, const spot::twa_graph_ptr& B,
                            const std::function<void(kofola::options&)>& customise_options = {}) {
    REQUIRE(A != nullptr);
    REQUIRE(B != nullptr);

    // Preserve caller options to avoid cross-test contamination
    auto saved_options = kofola::OPTIONS;

    // Configure reasonable defaults for inclusion (non-intrusive)
    // Users can adjust via OPTIONS if needed in future tests
    // e.g., enable early simulation pruning and optional preprocessing
    // Note: these params are optional; inclusion works without them
    // kofola::OPTIONS.params["early_sim"] = "yes";
    // kofola::OPTIONS.params["preproc_incl_A"] = "low";
    kofola::OPTIONS.params["preproc_incl_B"] = "low";
    kofola::OPTIONS.params["nac-alg"] = "subs_tup";
    kofola::OPTIONS.operation = "inclusion";

    if (customise_options) {
        customise_options(kofola::OPTIONS);
    }

    kofola::inclusion_check checker(A, B);
    bool result = checker.inclusion();

    // Restore previous options
    kofola::OPTIONS = saved_options;
    return result;
}
}

TEST_CASE("E2E inclusion: vector of automata pairs checked via Kofola", "[inclusion][e2e]") {
    // Build a vector of triples (A path, B path, expected inclusion result)
    const std::vector<std::tuple<std::string, std::string, bool>> test_cases = {
        {"tests/test_data/NI_correct_NI_formula_A.hoa", "tests/test_data/NI_correct_NI_formula_B.hoa", true},
        {"tests/test_data/bakery_3procs_bakery_formula_sym1_3proc_A.hoa", "tests/test_data/bakery_3procs_bakery_formula_sym1_3proc_B.hoa", false},
        {"tests/test_data/All_Sturmian_words_contain_cubes_sup.autfilt", "tests/test_data/All_Sturmian_words_contain_cubes_sub.autfilt", true}
    };

    for (const auto& tc : test_cases) {
        const auto& [a_path, b_path, expected] = tc;

        SECTION(std::string("Checking inclusion for pair: ") + a_path + " ⊆ " + b_path) {
            auto A = test_utils::load_automaton_from_file(a_path);
            auto B = test_utils::load_automaton_from_file(b_path);

            // Kofola inclusion decides emptiness of A ∩ ¬B; returns true iff inclusion holds
            bool a_subset_b = check_inclusion_kofola(A, B);
            CHECK(a_subset_b == expected);
        }
    }
}

TEST_CASE("Inclusion with early simulation enabled", "[inclusion][early_sim]") {
    const std::string a_path = "tests/test_data/AliasDarteFeautrierGonnord-SAS2010-nestedLoop_true-termination_true-no-overflow.c_Iteration2_A.ba.hoa";
    const std::string b_path = "tests/test_data/AliasDarteFeautrierGonnord-SAS2010-nestedLoop_true-termination_true-no-overflow.c_Iteration2_B.ba.hoa";

    auto A = test_utils::load_automaton_from_file(a_path);
    auto B = test_utils::load_automaton_from_file(b_path);

    // Early simulation pruning should not change the expected inclusion outcome for this pair
    bool a_subset_b = check_inclusion_kofola(B, A, [](kofola::options&opts) {
        opts.params["early_sim"] = "yes";
    });

    CHECK(a_subset_b == false);
}