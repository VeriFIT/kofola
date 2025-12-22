# Utils subdirectory CMakeLists.txt

# Ensure the same compiler settings as parent
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for utils
set(CXX_FLAGS)
list(APPEND CXX_FLAGS "-Wall")
list(APPEND CXX_FLAGS "-g")
list(APPEND CXX_FLAGS "-pedantic-errors")
list(APPEND CXX_FLAGS "-Wextra")
list(APPEND CXX_FLAGS "-fPIC")
list(APPEND CXX_FLAGS "-fno-strict-aliasing")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/algorithms)
include_directories(${CMAKE_SOURCE_DIR}/src/complement)
include_directories(${CMAKE_SOURCE_DIR}/src/inclusion)
include_directories(${CMAKE_SOURCE_DIR}/src/types)
include_directories(${CMAKE_SOURCE_DIR}/src/util)

# Define the same KOFOLA_SOURCES as parent with new paths
set(KOFOLA_SOURCES
    ${CMAKE_SOURCE_DIR}/src/algorithms/abstract_complement_alg.cpp
    ${CMAKE_SOURCE_DIR}/src/util/helpers.cpp
    ${CMAKE_SOURCE_DIR}/src/util/util.cpp
    ${CMAKE_SOURCE_DIR}/src/util/sets.cpp
    ${CMAKE_SOURCE_DIR}/src/complement/complement_tela.cpp
    ${CMAKE_SOURCE_DIR}/src/complement/complement_sync.cpp
    ${CMAKE_SOURCE_DIR}/src/complement/elevatorization.cpp
    ${CMAKE_SOURCE_DIR}/src/complement/decomposer.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_init_det.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_mh.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_ncsb.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_ncsb_delay.cpp
    # Excluding rank algorithms due to compilation issues
    # ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_rank.cpp
    # ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_rank2.cpp
    # ${CMAKE_SOURCE_DIR}/src/types/ranking.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_safra.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_sd_tela.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_sd_inductive.cpp
    ${CMAKE_SOURCE_DIR}/src/algorithms/complement_alg_subs_tuple.cpp
    ${CMAKE_SOURCE_DIR}/src/inclusion/emptiness_check.cpp
    ${CMAKE_SOURCE_DIR}/src/inclusion/inclusion_check.cpp
)

# Create test utilities library
add_library(test_utils STATIC
    test_utils.cpp
    ${KOFOLA_SOURCES}
)

target_compile_options(test_utils PRIVATE ${CXX_FLAGS})
target_link_libraries(test_utils 
    PRIVATE 
    ${SPOT}
    ${BDDX}
)

# Create standalone test executable for complement equivalence
add_executable(test_complement_equivalence
    test_complement_equivalence.cpp
)

target_compile_options(test_complement_equivalence PRIVATE ${CXX_FLAGS})

# Link with test utilities library
target_link_libraries(test_complement_equivalence 
    PRIVATE 
    test_utils
    ${SPOT}
    ${BDDX}
)
