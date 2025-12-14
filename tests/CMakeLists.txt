cmake_minimum_required(VERSION 3.16)

project(kofola_tests
    LANGUAGES CXX
    VERSION 1.0
)

if (APPLE) # we need to set correct path for Spot on Mac OS/homebrew
    execute_process(COMMAND brew --prefix
        OUTPUT_VARIABLE BREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("-- Brew prefix: ${BREW_PREFIX}")
    include_directories("${BREW_PREFIX}/include")
    link_directories("${BREW_PREFIX}/lib")
endif()

# Test configuration
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for tests
set(CXX_FLAGS)
list(APPEND CXX_FLAGS "-Wall")
list(APPEND CXX_FLAGS "-g")
list(APPEND CXX_FLAGS "-pedantic-errors")
list(APPEND CXX_FLAGS "-Wextra")
list(APPEND CXX_FLAGS "-fPIC")
list(APPEND CXX_FLAGS "-fno-strict-aliasing")

# Find or fetch Catch2
find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found, fetching from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0  # or any newer version
    )
    FetchContent_MakeAvailable(Catch2)
endif()

# Find required libraries (same as main project)
find_library(BDDX bddx)
find_library(SPOT spot)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# Test source files
set(TEST_SOURCES
    test_main.cpp
    types/test_binary_tree.cpp
    inductive/test_check_macrostate_from_acc_code.cpp
    tela/test_acc_code_dnf.cpp
    tela/test_complement_tela.cpp
    tela/test_guess_safe_models.cpp
    tela/test_contains_transition_color.cpp
    tela/test_create_partitions.cpp
    modular/test_complement_modular.cpp
    modular/test_scc_types.cpp
    inclusion/test_inclusion.cpp
)

# Create test sources that need to be compiled from main project
set(KOFOLA_SOURCES
    ${CMAKE_SOURCE_DIR}/src/abstract_complement_alg.cpp
    ${CMAKE_SOURCE_DIR}/src/util/helpers.cpp
    ${CMAKE_SOURCE_DIR}/src/util.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_tela.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_sync.cpp
    ${CMAKE_SOURCE_DIR}/src/decomposer.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_init_det.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_mh.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_ncsb.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_ncsb_delay.cpp
    # Excluding rank algorithms due to compilation issues
    # ${CMAKE_SOURCE_DIR}/src/complement_alg_rank.cpp
    # ${CMAKE_SOURCE_DIR}/src/complement_alg_rank2.cpp
    # ${CMAKE_SOURCE_DIR}/src/ranking.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_safra.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_sd_tela.cpp
    ${CMAKE_SOURCE_DIR}/src/complement_alg_subs_tuple.cpp
    ${CMAKE_SOURCE_DIR}/src/emptiness_check.cpp
    ${CMAKE_SOURCE_DIR}/src/inclusion_check.cpp
)

# Add the utils subdirectory
add_subdirectory(utils)

# Create test executable
add_executable(kofola_tests
    ${TEST_SOURCES}
)

target_compile_options(kofola_tests PRIVATE ${CXX_FLAGS})

# Link libraries including test utilities
target_link_libraries(kofola_tests 
    PRIVATE 
    Catch2::Catch2WithMain
    test_utils
    ${SPOT}
    ${BDDX}
)

# Include catch2
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
if(Catch2_FOUND)
    include(CTest)
    include(Catch)
    catch_discover_tests(kofola_tests)
endif()
