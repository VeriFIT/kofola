#include <catch2/catch_test_macros.hpp>
#include "complement_alg_sd_tela.hpp"
#include <set>
#include <vector>

using namespace kofola;
using namespace sd_tela;

TEST_CASE("guess_safe_models returns correct safe models for simple input", "[guess_safe_models]") {
    std::set<unsigned> check = {1, 2};
    std::vector<std::set<unsigned>> safe_models = {{}, {}};
    std::set<unsigned> breakpoint = {};
    unsigned model_index = 0;
    unsigned inf_index = 0;
    bool active = false;
    mstate_sd_tela init(check, safe_models, breakpoint, model_index, inf_index, active);

    std::set<unsigned> states = {3, 4};
    unsigned num_models = 2;
    auto result = guess_safe_models(init, states, num_models);

    // There should be 4 possible assignments (2^2 for 2 states, 2 models)
    REQUIRE(result.size() == 4);
    // Each result should have both states assigned to some model
    for (const auto& macro : result) {
        bool found_3 = false, found_4 = false;
        for (const auto& s : macro.safe_models_) {
            if (s.count(3)) found_3 = true;
            if (s.count(4)) found_4 = true;
        }
        REQUIRE(found_3);
        REQUIRE(found_4);
    }
}

TEST_CASE("guess_safe_models with empty states returns initial macrostate", "[guess_safe_models]") {
    std::set<unsigned> check = {1};
    std::vector<std::set<unsigned>> safe_models = {{}};
    std::set<unsigned> breakpoint = {};
    unsigned model_index = 0;
    unsigned inf_index = 0;
    bool active = true;
    mstate_sd_tela init(check, safe_models, breakpoint, model_index, inf_index, active);

    std::set<unsigned> states = {};
    unsigned num_models = 1;
    auto result = guess_safe_models(init, states, num_models);
    REQUIRE(result.size() == 1);
    // Should be identical to init
    REQUIRE(result[0].check_ == init.check_);
    REQUIRE(result[0].safe_models_ == init.safe_models_);
    REQUIRE(result[0].breakpoint_ == init.breakpoint_);
    REQUIRE(result[0].model_index_ == init.model_index_);
    REQUIRE(result[0].inf_index_ == init.inf_index_);
    REQUIRE(result[0].active_ == init.active_);
}
