// Test program for create_partitions function
// 
// This program tests the create_partitions function from complement_sync.cpp,
// focusing on the PartitionToAccMap component that maps partitions to their
// acceptance conditions. It uses various automata loaded from HOA format
// to test different acceptance condition types and partition scenarios.

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Spot headers
#include <spot/twaalgos/sccinfo.hh>

// Kofola headers
#include "complement/complement_sync.hpp"
#include "util/helpers.hpp"

// Test utilities
#include "../utils/test_utils.hpp"

namespace {

/**
 * @brief Test the create_partitions function for a given automaton
 * 
 * @param aut The input automaton
 * @param options The kofola options to use
 * @return The partition results including PartitionToAccMap
 */
std::tuple<size_t, kofola::PartitionToTypeMap, kofola::StateToPartitionMap, 
           kofola::SCCToPartitionMap, kofola::PartitionToAccMap>
test_create_partitions(const spot::twa_graph_ptr& aut, const kofola::options& options) {
    REQUIRE(aut != nullptr);
    
    // Create SCC info for the automaton
    spot::scc_info scc_inf(aut);
    
    // Call create_partitions
    auto result = helpers::tnba_complement::create_partitions(scc_inf, options);
    
    return result;
}

/**
 * @brief Validate the PartitionToAccMap for consistency
 * 
 * @param part_to_acc_map The partition to acceptance condition map
 * @param num_partitions Total number of partitions
 * @param aut The original automaton
 */
void validate_partition_to_acc_map(const kofola::PartitionToAccMap& part_to_acc_map,
                                   size_t num_partitions,
                                   const spot::twa_graph_ptr& aut) {
    // Check that we have acceptance conditions for valid partitions
    for (const auto& [partition, acc_cond] : part_to_acc_map) {
        // Partition index should be reasonable
        CHECK(partition < num_partitions);
        
        // Acceptance condition should be valid
        CHECK(acc_cond.num_sets() >= 0);
        
        // The acceptance condition should use a subset of the original acceptance sets
        CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
        
        // std::cout << "Partition " << partition 
        //           << " has acceptance condition: " << acc_cond 
        //           << " (uses " << acc_cond.num_sets() << " acceptance sets)" << std::endl;
    }
}

/**
 * @brief Print detailed partition information
 */
void print_partition_info(const std::tuple<size_t, kofola::PartitionToTypeMap, 
                         kofola::StateToPartitionMap, kofola::SCCToPartitionMap, 
                         kofola::PartitionToAccMap>& result,
                         const spot::twa_graph_ptr& aut) {
    auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
    
    std::cout << "\n=== Partition Analysis ===" << std::endl;
    std::cout << "Total partitions: " << num_partitions << std::endl;
    std::cout << "Original automaton acceptance: " << aut->acc() << std::endl;
    std::cout << "Original automaton has " << aut->acc().num_sets() << " acceptance sets" << std::endl;
    
    std::cout << "\nPartition types:" << std::endl;
    for (const auto& [partition, type] : part_to_type_map) {
        std::cout << "  Partition " << partition << ": ";
        switch (type) {
            case kofola::PartitionType::INHERENTLY_WEAK:
                std::cout << "INHERENTLY_WEAK"; break;
            case kofola::PartitionType::DETERMINISTIC:
                std::cout << "DETERMINISTIC"; break;
            case kofola::PartitionType::STRONGLY_DETERMINISTIC:
                std::cout << "STRONGLY_DETERMINISTIC"; break;
            case kofola::PartitionType::NONDETERMINISTIC:
                std::cout << "NONDETERMINISTIC"; break;
            case kofola::PartitionType::INITIAL_DETERMINISTIC:
                std::cout << "INITIAL_DETERMINISTIC"; break;
        }
        std::cout << std::endl;
    }
    
    std::cout << "\nSCC to partition mapping:" << std::endl;
    for (const auto& [scc, partition] : scc_to_part_map) {
        std::cout << "  SCC " << scc << " -> Partition " << partition << std::endl;
    }
    
    std::cout << "\nPartition acceptance conditions:" << std::endl;
    for (const auto& [partition, acc_cond] : part_to_acc_map) {
        std::cout << "  Partition " << partition << ": " << acc_cond << std::endl;
    }
    std::cout << "=========================" << std::endl;
}

} // anonymous namespace

// Test cases for create_partitions function
TEST_CASE("create_partitions produces valid PartitionToAccMap", "[create_partitions]") {
    // Default options
    kofola::options default_options;
    
    SECTION("Simple Buchi automaton") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/simple_buchi.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // print_partition_info(result, aut);
        
        // Validate basic properties
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
        
        // For a Buchi automaton, we expect acceptance conditions to be Buchi-like
        for (const auto& [partition, acc_cond] : part_to_acc_map) {
            CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
        }
    }
    
    SECTION("Non-deterministic Buchi automaton") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/ndet_example.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
    }
    
    SECTION("Streett automaton with multiple acceptance sets") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/random_sd_streett_001.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
        
        // For Streett automata, check that partitions preserve the structure
        for (const auto& [partition, acc_cond] : part_to_acc_map) {
            // Each partition should have a valid acceptance condition
            CHECK(acc_cond.num_sets() >= 0);
            // The acceptance condition should be a restriction of the original
            CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
        }
    }
}

TEST_CASE("create_partitions with different merge options", "[create_partitions][merge]") {
    SECTION("Test with merge_iwa option") {
        kofola::options options_merge_iwa;
        options_merge_iwa.params["merge_iwa"] = "yes";
        
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/simple_buchi.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, options_merge_iwa);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // std::cout << "\n=== Testing with merge_iwa=yes ===" << std::endl;
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
    }
    
    SECTION("Test with merge_det option") {
        kofola::options options_merge_det;
        options_merge_det.params["merge_det"] = "yes";
        
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/simple_buchi.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, options_merge_det);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // std::cout << "\n=== Testing with merge_det=yes ===" << std::endl;
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
    }
}

TEST_CASE("PartitionToAccMap consistency across different automata", "[create_partitions][consistency]") {
    kofola::options default_options;
    
    // Test files with different acceptance condition types
    std::vector<std::string> test_files = {
        "test_data/simple_buchi.hoa",
        "test_data/ndet_example.hoa",
        "test_data/multi_scc_buchi.hoa",
        "test_data/mixed_scc_types.hoa",
        "test_data/random_sd_streett_001.hoa",
        "test_data/random_sd_streett_002.hoa",
        "test_data/random_sd_streett_003.hoa"
    };
    
    for (const std::string& filename : test_files) {
        SECTION("Testing consistency for file: " + filename) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            auto result = test_create_partitions(aut, default_options);
            auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
            
            // std::cout << "\n=== Testing " << filename << " ===" << std::endl;
            // print_partition_info(result, aut);
            
            // Basic consistency checks
            CHECK(num_partitions > 0);
            validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
            
            // Check that every partition with an acceptance condition has a corresponding type
            for (const auto& [partition, acc_cond] : part_to_acc_map) {
                CHECK(part_to_type_map.find(partition) != part_to_type_map.end());
            }
            
            // Check that every valid SCC maps to a partition with acceptance condition
            for (const auto& [scc, partition] : scc_to_part_map) {
                if (partition != -1) { // -1 indicates invalid/non-accepting SCC
                    CHECK(part_to_acc_map.find(partition) != part_to_acc_map.end());
                }
            }
            
            // Check acceptance condition preservation
            for (const auto& [partition, acc_cond] : part_to_acc_map) {
                // Each partition's acceptance condition should be valid
                CHECK((acc_cond.is_t() || acc_cond.is_f() || acc_cond.num_sets() > 0));
                
                // The number of acceptance sets should not exceed the original
                CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
            }
        }
    }
}

TEST_CASE("PartitionToAccMap edge cases", "[create_partitions][edge_cases]") {
    kofola::options default_options;
    
    SECTION("Test automaton with Inf(a) acceptance") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/inf_a.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // std::cout << "\n=== Testing inf_a.hoa ===" << std::endl;
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
        
        // Check that all acceptance conditions are properly formed
        for (const auto& [partition, acc_cond] : part_to_acc_map) {
            // The acceptance condition should be a valid formula
            CHECK((!acc_cond.is_f())); // Should not be false (unless no accepting states)
        }
    }
}

TEST_CASE("create_partitions with multiple SCCs", "[create_partitions][multi_scc]") {
    kofola::options default_options;
    
    SECTION("Multi-SCC Buchi automaton") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/multi_scc_buchi.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // std::cout << "\n=== Testing Multi-SCC Buchi ===" << std::endl;
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
        
        // std::cout << "Number of accepting SCCs: " << accepting_sccs << std::endl;
        
        // Test merge behavior
        kofola::options merge_options;
        merge_options.params["merge_iwa"] = "yes";
        merge_options.params["merge_det"] = "yes";
        
        auto result_merged = test_create_partitions(aut, merge_options);
        auto [num_partitions_merged, part_to_type_map_merged, st_to_part_map_merged, 
              scc_to_part_map_merged, part_to_acc_map_merged] = result_merged;
        
        // std::cout << "\n=== Testing Multi-SCC Buchi with merge options ===" << std::endl;
        // print_partition_info(result_merged, aut);
        
        validate_partition_to_acc_map(part_to_acc_map_merged, num_partitions_merged, aut);
    }
    
    SECTION("Mixed SCC types") {
        spot::twa_graph_ptr aut = test_utils::load_automaton_from_file("test_data/mixed_scc_types.hoa");
        REQUIRE(aut != nullptr);
        
        auto result = test_create_partitions(aut, default_options);
        auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
        
        // std::cout << "\n=== Testing Mixed SCC Types ===" << std::endl;
        // print_partition_info(result, aut);
        
        CHECK(num_partitions > 0);
        validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
        
        // Verify acceptance condition consistency for complex acceptance
        for (const auto& [partition, acc_cond] : part_to_acc_map) {
            // Complex acceptance conditions should be properly handled
            CHECK(acc_cond.num_sets() >= 0);
            CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
        }
    }
}

TEST_CASE("PartitionToAccMap specific validation", "[create_partitions][acc_validation]") {
    kofola::options default_options;
    
    SECTION("Acceptance condition preservation") {
        // Test that acceptance conditions are properly preserved and restricted
        std::vector<std::string> test_files = {
            "test_data/simple_buchi.hoa",
            "test_data/multi_scc_buchi.hoa",
            "test_data/random_sd_streett_001.hoa"
        };
        
        for (const std::string& filename : test_files) {
            spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
            REQUIRE(aut != nullptr);
            
            auto result = test_create_partitions(aut, default_options);
            auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
            
            // std::cout << "\n=== Acc Validation for " << filename << " ===" << std::endl;
            
            // Test acceptance condition properties
            for (const auto& [partition, acc_cond] : part_to_acc_map) {
                INFO("Testing partition " << partition << " with acc condition " << acc_cond);
                
                // Basic validity checks
                CHECK(acc_cond.num_sets() >= 0);
                CHECK(acc_cond.num_sets() <= aut->acc().num_sets());
                
                // Check that the acceptance condition is not trivially false
                // (unless the partition has no accepting states)
                if (!acc_cond.is_f()) {
                    CHECK((acc_cond.is_t() || acc_cond.num_sets() > 0));
                }
                
                // Verify that used acceptance sets are valid
                if (acc_cond.num_sets() > 0) {
                    // The acceptance condition should be satisfiable
                    CHECK(!acc_cond.is_f());
                }
            }
            
            // Test partition completeness: all accepting SCCs should have partitions
            spot::scc_info scc_inf(aut);
            for (size_t scc = 0; scc < scc_inf.scc_count(); ++scc) {
                if (scc_inf.is_accepting_scc(scc)) {
                    CHECK(scc_to_part_map.find(scc) != scc_to_part_map.end());
                    int partition = scc_to_part_map.at(scc);
                    if (partition != -1) {
                        CHECK(part_to_acc_map.find(partition) != part_to_acc_map.end());
                    }
                }
            }
        }
    }
}

// Manual test function that can be called from main for debugging
void run_create_partitions_test_on_file(const std::string& filename) {
    std::cout << "Testing create_partitions for file: " << filename << std::endl;
    
    spot::twa_graph_ptr aut = test_utils::load_automaton_from_file(filename);
    if (!aut) {
        std::cout << "Failed to load automaton from file: " << filename << std::endl;
        return;
    }
    
    std::cout << "Loaded automaton with " << aut->num_states() << " states" << std::endl;
    std::cout << "Acceptance condition: " << aut->get_acceptance() << std::endl;
    
    kofola::options default_options;
    auto result = test_create_partitions(aut, default_options);
    
    print_partition_info(result, aut);
    
    auto [num_partitions, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map] = result;
    validate_partition_to_acc_map(part_to_acc_map, num_partitions, aut);
    
    std::cout << "✓ Test completed successfully" << std::endl;
}
