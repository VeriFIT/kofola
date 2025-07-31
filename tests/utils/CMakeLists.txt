# Utils subdirectory CMakeLists.txt

# Ensure the same compiler settings as parent
set(CMAKE_CXX_STANDARD 17)
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

# Define the same KOFOLA_SOURCES as parent
set(KOFOLA_SOURCES
    ${CMAKE_SOURCE_DIR}/src/abstract_complement_alg.cpp
    ${CMAKE_SOURCE_DIR}/src/kofola.cpp
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

# Create standalone test executable for complement equivalence
add_executable(test_complement_equivalence
    test_complement_equivalence.cpp
    ${KOFOLA_SOURCES}
)

target_compile_options(test_complement_equivalence PRIVATE ${CXX_FLAGS})

# Link libraries for standalone test
target_link_libraries(test_complement_equivalence 
    PRIVATE 
    ${SPOT}
    ${BDDX}
)
