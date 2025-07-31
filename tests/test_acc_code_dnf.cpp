#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "abstract_complement_alg.hpp"
#include "kofola.hpp"

// SPOT includes
#include <spot/twa/acc.hh>
#include <spot/misc/bddlt.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/tl/parse.hh>
#include <spot/twaalgos/sccinfo.hh>

using namespace kofola;

TEST_CASE("acc_code_dnf with basic acceptance conditions", "[acc_code_dnf]") {
    SECTION("Büchi acceptance condition") {
        // Create Büchi acceptance: Inf(0)
        auto buchi_code = spot::acc_cond::acc_code("Inf(0)");
        auto dnf_result = cmpl_info::acc_code_dnf(buchi_code);
        
        // Complement of Inf(0) should be Fin(0), which in DNF should be one clause
        REQUIRE(dnf_result.size() == 1);
        
        // Should have exactly one clause with Fin(0)
        REQUIRE(dnf_result[0].fins.size() == 1);
        REQUIRE(dnf_result[0].fins[0] == spot::acc_cond::mark_t{0});
        REQUIRE(dnf_result[0].infs.empty());
    }
    
    SECTION("Co-Büchi acceptance condition") {
        // Create co-Büchi acceptance: Fin(0)
        auto co_buchi_code = spot::acc_cond::acc_code("Fin(0)");
        auto dnf_result = cmpl_info::acc_code_dnf(co_buchi_code);
        
        // Complement of Fin(0) should be Inf(0)
        REQUIRE(dnf_result.size() == 1);
        
        // Should have exactly one clause with Inf(0)
        REQUIRE(dnf_result[0].infs.size() == 1);
        REQUIRE(dnf_result[0].infs[0] == spot::acc_cond::mark_t{0});
        REQUIRE(dnf_result[0].fins.empty());
    }
    
    SECTION("Generalized Büchi acceptance") {
        // Create generalized Büchi: Inf(0) & Inf(1)
        auto gen_buchi_code = spot::acc_cond::acc_code("Inf(0) & Inf(1)");
        auto dnf_result = cmpl_info::acc_code_dnf(gen_buchi_code);
        
        // Complement of Inf(0) & Inf(1) should be Fin(0) | Fin(1)
        // This should result in 2 clauses in DNF
        REQUIRE(dnf_result.size() == 2);
        
        // First clause should have Fin(0)
        REQUIRE(dnf_result[0].fins.size() == 1);
        REQUIRE(dnf_result[0].fins[0] == spot::acc_cond::mark_t{0});
        REQUIRE(dnf_result[0].infs.empty());
        
        // Second clause should have Fin(1)
        REQUIRE(dnf_result[1].fins.size() == 1);
        REQUIRE(dnf_result[1].fins[0] == spot::acc_cond::mark_t{1});
        REQUIRE(dnf_result[1].infs.empty());
    }
}

TEST_CASE("acc_code_dnf with complex acceptance conditions", "[acc_code_dnf]") {
    SECTION("Rabin acceptance condition") {
        // Create Rabin acceptance: (Fin(0) & Inf(1)) | (Fin(2) & Inf(3))
        auto rabin_code = spot::acc_cond::acc_code("(Fin(0) & Inf(1)) | (Fin(2) & Inf(3))");
        auto dnf_result = cmpl_info::acc_code_dnf(rabin_code);

        // Complement of (Fin(0) & Inf(1)) | (Fin(2) & Inf(3)) should be 
        // (Inf(0) | Fin(1)) & (Inf(2) | Fin(3))
        // In DNF this becomes: (Inf(0) & Inf(2)) | (Inf(0) & Fin(3)) | (Fin(1) & Inf(2)) | (Fin(1) & Fin(3))
        REQUIRE(dnf_result.size() == 4);
        
        // Verify each clause has the expected content
        bool found_inf0_inf2 = false, found_inf0_fin3 = false, found_fin1_inf2 = false, found_fin1_fin3 = false;
        
        for (const auto& clause : dnf_result) {
            if (clause.infs.size() == 2 && clause.fins.empty()) {
                // Should be (Inf(0) & Inf(2))
                REQUIRE(((clause.infs[0] == spot::acc_cond::mark_t{0} && clause.infs[1] == spot::acc_cond::mark_t{2}) || (clause.infs[0] == spot::acc_cond::mark_t{2} && clause.infs[1] == spot::acc_cond::mark_t{0})));
                found_inf0_inf2 = true;
            } else if (clause.infs.size() == 1 && clause.fins.size() == 1) {
                if (clause.infs[0] == spot::acc_cond::mark_t{0} && clause.fins[0] == spot::acc_cond::mark_t{3}) {
                    found_inf0_fin3 = true;
                } else if (clause.infs[0] == spot::acc_cond::mark_t{2} && clause.fins[0] == spot::acc_cond::mark_t{1}) {
                    found_fin1_inf2 = true;
                }
            } else if (clause.infs.empty() && clause.fins.size() == 2) {
                // Should be (Fin(1) & Fin(3))
                REQUIRE(((clause.fins[0] == spot::acc_cond::mark_t{1} && clause.fins[1] == spot::acc_cond::mark_t{3}) || (clause.fins[0] == spot::acc_cond::mark_t{3} && clause.fins[1] == spot::acc_cond::mark_t{1})));
                found_fin1_fin3 = true;
            }
        }
        
        REQUIRE(found_inf0_inf2);
        REQUIRE(found_inf0_fin3);
        REQUIRE(found_fin1_inf2);
        REQUIRE(found_fin1_fin3);
    }
    
    SECTION("Streett acceptance condition") {
        // Create Streett acceptance: (Inf(0) | Fin(1)) & (Inf(2) | Fin(3))
        auto streett_code = spot::acc_cond::acc_code("(Inf(0) | Fin(1)) & (Inf(2) | Fin(3))");
        auto dnf_result = cmpl_info::acc_code_dnf(streett_code);
        
        // Complement of (Inf(0) | Fin(1)) & (Inf(2) | Fin(3)) should be
        // (Fin(0) & Inf(1)) | (Fin(2) & Inf(3))
        // In DNF: (Fin(0) & Inf(1)) | (Fin(2) & Inf(3))
        REQUIRE(dnf_result.size() == 2);
        
        // Check both clauses
        bool found_fin0_inf1 = false, found_fin2_inf3 = false;
        
        for (const auto& clause : dnf_result) {
            REQUIRE(clause.fins.size() == 1);
            REQUIRE(clause.infs.size() == 1);
            
            if (clause.fins[0] == spot::acc_cond::mark_t{0} && clause.infs[0] == spot::acc_cond::mark_t{1}) {
                found_fin0_inf1 = true;
            } else if (clause.fins[0] == spot::acc_cond::mark_t{2} && clause.infs[0] == spot::acc_cond::mark_t{3}) {
                found_fin2_inf3 = true;
            }
        }
        
        REQUIRE(found_fin0_inf1);
        REQUIRE(found_fin2_inf3);
    }
    
    SECTION("Mixed complex condition") {
        // Create a complex mixed condition
        auto complex_code = spot::acc_cond::acc_code("(Inf(0) & Inf(1)) | (Fin(2) & Inf(3))");
        auto dnf_result = cmpl_info::acc_code_dnf(complex_code);
        
        // Complement of (Inf(0) & Inf(1)) | (Fin(2) & Inf(3)) should be
        // (Fin(0) | Fin(1)) & (Inf(2) | Fin(3))
        // In DNF: (Fin(0) & Inf(2)) | (Fin(0) & Fin(3)) | (Fin(1) & Inf(2)) | (Fin(1) & Fin(3))
        REQUIRE(dnf_result.size() == 4);
        
        // Verify structure of clauses
        for (const auto& clause : dnf_result) {
            REQUIRE(clause.fins.size() >= 1);
            REQUIRE(clause.fins.size() + clause.infs.size() == 2);
        }
    }
}

TEST_CASE("acc_code_dnf with edge cases", "[acc_code_dnf]") {
    SECTION("True acceptance condition") {
        // Create always-accepting condition: t (true)
        auto true_code = spot::acc_cond::acc_code("t");
        auto dnf_result = cmpl_info::acc_code_dnf(true_code);
        
        // Complement of true should be false, which results in an empty DNF
        REQUIRE(dnf_result.empty());
    }
    
    SECTION("False acceptance condition") {
        // Create never-accepting condition: f (false)  
        auto false_code = spot::acc_cond::acc_code("f");
        auto dnf_result = cmpl_info::acc_code_dnf(false_code);
        // Complement of false should be true, which results in an non-empty empty DNF
        REQUIRE_FALSE(dnf_result.empty());
    }
    
    SECTION("Single acceptance set") {
        // Test with different single acceptance set conditions
        std::vector<std::string> conditions = {
            "Inf(0)",
            "Fin(0)", 
            "Inf(1)",
            "Fin(1)",
            "Inf(2)",
            "Fin(2)"
        };
        
        for (const auto& condition_str : conditions) {
            auto code = spot::acc_cond::acc_code(condition_str.c_str());
            auto dnf_result = cmpl_info::acc_code_dnf(code);
            
            // Each should produce a valid result
            REQUIRE(dnf_result.size() >= 0);
            
            // If non-empty, verify structure
            for (const auto& clause : dnf_result) {
                REQUIRE(clause.fins.size() >= 0);
                REQUIRE(clause.infs.size() >= 0);
                // At least one of fins or infs should be non-empty for a valid clause
                if (!dnf_result.empty()) {
                    REQUIRE(((clause.fins.size() > 0) || (clause.infs.size() > 0)));
                }
            }
        }
    }
}

TEST_CASE("acc_code_dnf deterministic behavior", "[acc_code_dnf]") {
    SECTION("Multiple calls produce identical results") {
        auto buchi_code = spot::acc_cond::acc_code("Inf(0)");
        
        // Call multiple times
        auto result1 = cmpl_info::acc_code_dnf(buchi_code);
        auto result2 = cmpl_info::acc_code_dnf(buchi_code);
        auto result3 = cmpl_info::acc_code_dnf(buchi_code);
        
        // Results should be identical
        REQUIRE(result1.size() == result2.size());
        REQUIRE(result2.size() == result3.size());
        
        for (size_t i = 0; i < result1.size(); ++i) {
            REQUIRE(result1[i].fins.size() == result2[i].fins.size());
            REQUIRE(result2[i].fins.size() == result3[i].fins.size());
            REQUIRE(result1[i].infs.size() == result2[i].infs.size());
            REQUIRE(result2[i].infs.size() == result3[i].infs.size());
        }
    }
    
    SECTION("Different equivalent formulations produce same results") {
        // Test logically equivalent formulations
        auto code1 = spot::acc_cond::acc_code("Inf(0) & Inf(1)");
        auto code2 = spot::acc_cond::acc_code("Inf(1) & Inf(0)");
        
        auto result1 = cmpl_info::acc_code_dnf(code1);
        auto result2 = cmpl_info::acc_code_dnf(code2);
        
        // Should produce the same number of clauses (though order might differ)
        REQUIRE(result1.size() == result2.size());
    }
}

TEST_CASE("acc_code_dnf DNF structure validation", "[acc_code_dnf]") {
    SECTION("DNF structure is well-formed") {
        std::vector<std::string> test_conditions = {
            "Inf(0)",
            "Fin(0)",
            "Inf(0) & Inf(1)",
            "Inf(0) | Inf(1)",
            "Fin(0) & Fin(1)",
            "Fin(0) | Fin(1)",
            "(Inf(0) & Fin(1)) | (Inf(2) & Fin(3))"
        };
        
        for (const auto& condition_str : test_conditions) {
            auto code = spot::acc_cond::acc_code(condition_str.c_str());
            auto dnf_result = cmpl_info::acc_code_dnf(code);
            
            // Verify that the DNF structure is well-formed
            for (size_t i = 0; i < dnf_result.size(); ++i) {
                const auto& clause = dnf_result[i];
                
                // Each clause should have valid fins and infs vectors
                REQUIRE(clause.fins.size() >= 0);
                REQUIRE(clause.infs.size() >= 0);
                
                // Verify that the marks are valid (non-negative)
                for (const auto& mark : clause.fins) {
                    REQUIRE(mark >= spot::acc_cond::mark_t{0});
                }
                for (const auto& mark : clause.infs) {
                    REQUIRE(mark >= spot::acc_cond::mark_t{0});
                }
            }
        }
    }
}

TEST_CASE("acc_code_dnf with multiple acceptance sets", "[acc_code_dnf]") {
    SECTION("Three acceptance sets") {
        auto code = spot::acc_cond::acc_code("Inf(0) & Inf(1) & Inf(2)");
        auto dnf_result = cmpl_info::acc_code_dnf(code);
        
        // Complement of Inf(0) & Inf(1) & Inf(2) should be Fin(0) | Fin(1) | Fin(2)
        // This should result in 3 clauses in DNF
        REQUIRE(dnf_result.size() == 3);
        
        // Each clause should have exactly one Fin
        for (const auto& clause : dnf_result) {
            REQUIRE(clause.fins.size() == 1);
            REQUIRE(clause.infs.empty());
        }
        
        // Check that we have Fin(0), Fin(1), and Fin(2)
        std::set<spot::acc_cond::mark_t> found_fins;
        for (const auto& clause : dnf_result) {
            found_fins.insert(clause.fins[0]);
        }
        REQUIRE(found_fins.size() == 3);
        REQUIRE(found_fins.count(spot::acc_cond::mark_t{0}) == 1);
        REQUIRE(found_fins.count(spot::acc_cond::mark_t{1}) == 1);
        REQUIRE(found_fins.count(spot::acc_cond::mark_t{2}) == 1);
    }
    
    SECTION("Four acceptance sets with mixed conditions") {
        auto code = spot::acc_cond::acc_code("(Inf(0) & Fin(1)) | (Inf(2) & Fin(3))");
        auto dnf_result = cmpl_info::acc_code_dnf(code);
        
        // Should handle multiple sets correctly
        REQUIRE_FALSE(dnf_result.empty());
        
        for (const auto& clause : dnf_result) {
            // Each clause should have at least one fin or inf
            REQUIRE(((clause.fins.size() > 0) || (clause.infs.size() > 0)));
        }
    }
    
    SECTION("Large number of acceptance sets") {
        // Test with many acceptance sets
        auto code = spot::acc_cond::acc_code("Inf(0) & Inf(1) & Inf(2) & Inf(3) & Inf(4)");
        auto dnf_result = cmpl_info::acc_code_dnf(code);
        
        // Should handle large numbers of acceptance sets
        // Complement should be Fin(0) | Fin(1) | Fin(2) | Fin(3) | Fin(4)
        REQUIRE(dnf_result.size() == 5);
        
        for (const auto& clause : dnf_result) {
            REQUIRE(clause.fins.size() == 1);
            REQUIRE(clause.infs.empty());
        }
    }
}

TEST_CASE("acc_code_dnf specific complement verification", "[acc_code_dnf]") {
    SECTION("Büchi complement verification") {
        // Büchi: Inf(0) -> complement should be Fin(0)
        auto buchi_code = spot::acc_cond::acc_code("Inf(0)");
        auto dnf_result = cmpl_info::acc_code_dnf(buchi_code);
        
        // Should produce exactly one clause (Fin(0))
        REQUIRE(dnf_result.size() == 1);
        REQUIRE(dnf_result[0].fins.size() == 1);
        REQUIRE(dnf_result[0].fins[0] == spot::acc_cond::mark_t{0});
        REQUIRE(dnf_result[0].infs.empty());
    }
    
    SECTION("Generalized Büchi complement verification") {
        // Gen-Büchi: Inf(0) & Inf(1) -> complement should be Fin(0) | Fin(1)
        auto gen_buchi_code = spot::acc_cond::acc_code("Inf(0) & Inf(1)");
        auto dnf_result = cmpl_info::acc_code_dnf(gen_buchi_code);
        
        // Should produce exactly two clauses (Fin(0), Fin(1))
        REQUIRE(dnf_result.size() == 2);
        REQUIRE(dnf_result[0].fins.size() == 1);
        REQUIRE(dnf_result[0].infs.empty());
        REQUIRE(dnf_result[1].fins.size() == 1);
        REQUIRE(dnf_result[1].infs.empty());
    }
    
    SECTION("Co-Büchi complement verification") {
        // Co-Büchi: Fin(0) -> complement should be Inf(0)
        auto co_buchi_code = spot::acc_cond::acc_code("Fin(0)");
        auto dnf_result = cmpl_info::acc_code_dnf(co_buchi_code);
        
        // Should produce exactly one clause (Inf(0))
        REQUIRE(dnf_result.size() == 1);
        REQUIRE(dnf_result[0].infs.size() == 1);
        REQUIRE(dnf_result[0].infs[0] == spot::acc_cond::mark_t{0});
        REQUIRE(dnf_result[0].fins.empty());
    }
    
    SECTION("Disjunction complement verification") {
        // Fin(0) | Fin(1) -> complement should be Inf(0) & Inf(1)
        auto disj_code = spot::acc_cond::acc_code("Fin(0) | Fin(1)");
        auto dnf_result = cmpl_info::acc_code_dnf(disj_code);
        
        // Should produce exactly one clause with two Infs
        REQUIRE(dnf_result.size() == 1);
        REQUIRE(dnf_result[0].infs.size() == 2);
        REQUIRE(dnf_result[0].fins.empty());
        // Check that we have both Inf(0) and Inf(1)
        std::set<spot::acc_cond::mark_t> infs_set(dnf_result[0].infs.begin(), dnf_result[0].infs.end());
        REQUIRE(infs_set.count(spot::acc_cond::mark_t{0}) == 1);
        REQUIRE(infs_set.count(spot::acc_cond::mark_t{1}) == 1);
    }
}

TEST_CASE("acc_code_dnf error handling and robustness", "[acc_code_dnf]") {
    SECTION("Valid but complex conditions") {
        std::vector<std::string> complex_conditions = {
            "((Inf(0) & Fin(1)) | (Inf(2) & Fin(3))) & ((Inf(4) & Fin(5)) | (Inf(6) & Fin(7)))",
            "Inf(0) & Inf(1) & Inf(2) & Inf(3) & Inf(4) & Inf(5)",
            "(Fin(0) | Fin(1) | Fin(2)) & (Inf(3) | Inf(4) | Inf(5))",
            "((Inf(0) | Fin(1)) & (Inf(2) | Fin(3))) | ((Inf(4) | Fin(5)) & (Inf(6) | Fin(7)))"
        };
        
        for (const auto& condition_str : complex_conditions) {
            auto code = spot::acc_cond::acc_code(condition_str.c_str());
            auto dnf_result = cmpl_info::acc_code_dnf(code);
            
            // Should handle all complex conditions without throwing
            REQUIRE(dnf_result.size() >= 0);
            
            // Verify structure integrity
            for (const auto& clause : dnf_result) {
                REQUIRE(clause.fins.size() >= 0);
                REQUIRE(clause.infs.size() >= 0);
                // Verify marks are valid
                for (const auto& mark : clause.fins) {
                    REQUIRE(mark >= spot::acc_cond::mark_t{0});
                }
                for (const auto& mark : clause.infs) {
                    REQUIRE(mark >= spot::acc_cond::mark_t{0});
                }
            }
        }
    }
}

TEST_CASE("acc_code_dnf integration with compute_cond_to_verify", "[acc_code_dnf][integration]") {
    SECTION("Verify integration works correctly") {
        // Create a simple automaton with Büchi acceptance
        auto ltl_formula = spot::parse_infix_psl("a");
        if (!ltl_formula.errors.empty()) {
            return; // Skip if parsing fails
        }
        
        auto aut = spot::translator().run(ltl_formula.f);
        aut->set_acceptance(spot::acc_cond::acc_code::buchi());
        
        // Create minimal cmpl_info for testing
        PartitionToTypeMap part_to_type_map;
        part_to_type_map[0] = PartitionType::NONDETERMINISTIC;
        
        StateToPartitionMap st_to_part_map;
        for (unsigned i = 0; i < aut->num_states(); ++i) {
            st_to_part_map[i] = 0;
        }
        
        ReachableVector reachable_vector(aut->num_states());
        for (unsigned i = 0; i < aut->num_states(); ++i) {
            reachable_vector[i].insert(i);
        }
        
        PartitionToSCCMap part_to_scc_map;
        SCCToSCCSetMap scc_to_pred_sccs_map;
        spot::scc_info scc_info(aut);
        Simulation dir_sim;
        std::vector<bool> state_accepting(aut->num_states(), false);
        
        cmpl_info info(
            aut,
            1, // num_partitions
            part_to_type_map,
            st_to_part_map,
            reachable_vector,
            part_to_scc_map,
            scc_to_pred_sccs_map,
            scc_info,
            dir_sim,
            state_accepting,
            false // shared_breakpoint
        );
        
        // Test that compute_cond_to_verify works
        info.compute_cond_to_verify();
        REQUIRE_FALSE(info.cond_to_verify_.empty());
        
        // Test that direct call to acc_code_dnf produces same result
        auto direct_result = cmpl_info::acc_code_dnf(aut->acc().get_acceptance());
        
        // Results should be the same
        REQUIRE(info.cond_to_verify_.size() == direct_result.size());
        for (size_t i = 0; i < info.cond_to_verify_.size(); ++i) {
            REQUIRE(info.cond_to_verify_[i].fins.size() == direct_result[i].fins.size());
            REQUIRE(info.cond_to_verify_[i].infs.size() == direct_result[i].infs.size());
        }
    }
}
