# Test configuration
set(CMAKE_CXX_STANDARD 17)
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
    test_acc_code_dnf.cpp
)

# Create test sources that need to be compiled from main project
set(KOFOLA_SOURCES
    ${CMAKE_SOURCE_DIR}/src/abstract_complement_alg.cpp
    ${CMAKE_SOURCE_DIR}/src/kofola.cpp
    ${CMAKE_SOURCE_DIR}/src/util.cpp
)

# Create test executable
add_executable(kofola_tests
    ${TEST_SOURCES}
    ${KOFOLA_SOURCES}
)

target_compile_options(kofola_tests PRIVATE ${CXX_FLAGS})

# Link libraries
target_link_libraries(kofola_tests 
    PRIVATE 
    Catch2::Catch2WithMain
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
