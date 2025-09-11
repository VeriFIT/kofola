#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>
#include <iostream>

// Spot
#include <spot/twa/twa.hh>

// Kofola
#include "inclusion/inclusion_check.hpp"

// Test utils
#include "../utils/test_utils.hpp"

namespace {
// Simple harness that checks A ⊆ B using kofola::inclusion_check
bool check_inclusion_kofola(const spot::twa_graph_ptr& A, const spot::twa_graph_ptr& B) {
    REQUIRE(A != nullptr);
    REQUIRE(B != nullptr);

    // Configure reasonable defaults for inclusion (non-intrusive)
    // Users can adjust via OPTIONS if needed in future tests
    // e.g., enable early simulation pruning and optional preprocessing
    // Note: these params are optional; inclusion works without them
    // kofola::OPTIONS.params["early_sim"] = "yes";
    // kofola::OPTIONS.params["preproc_incl_A"] = "low";
    kofola::OPTIONS.params["preproc_incl_B"] = "low";
    kofola::OPTIONS.params["merge_iwa"] = "yes";
    kofola::OPTIONS.params["merge_det"] = "yes";
    kofola::OPTIONS.params["nac-alg"] = "subs_tup";

    kofola::inclusion_check checker(A, B);
    return checker.inclusion();
}
}

TEST_CASE("E2E inclusion: vector of automata pairs checked via Kofola", "[inclusion][e2e]") {
    // Build a vector of pairs (A, B) as file paths
    const std::vector<std::pair<std::string, std::string>> test_pairs = {
        {"tests/test_data/NI_correct_NI_formula_A.hoa", "tests/test_data/NI_correct_NI_formula_B.hoa"},
    };

    for (const auto& paths : test_pairs) {
        const auto& a_path = paths.first;
        const auto& b_path = paths.second;

        SECTION(std::string("Checking inclusion for pair: ") + a_path + " ⊆ " + b_path) {
            auto A = test_utils::load_automaton_from_file(a_path);
            auto B = test_utils::load_automaton_from_file(b_path);

            // Kofola inclusion decides emptiness of A ∩ ¬B; returns true iff inclusion holds
            bool a_subset_b = check_inclusion_kofola(A, B);
            CHECK(a_subset_b);
        }
    }
}
