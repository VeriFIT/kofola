cmake_minimum_required(VERSION 3.16)

set(CMAKE_COLOR_MAKEFILE ON)

project(kofola
	LANGUAGES CXX
	VERSION 1.0
)

# Create a custom target to generate git hash at build time
# This runs every time make is called to ensure fresh git hash
add_custom_target(generate_version ALL
	COMMAND ${CMAKE_COMMAND}
		-DSOURCE_DIR=${CMAKE_SOURCE_DIR}
		-DBINARY_DIR=${CMAKE_BINARY_DIR}
		-P ${CMAKE_SOURCE_DIR}/cmake/GenerateVersion.cmake
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	COMMENT "Generating version header with current git hash"
	VERBATIM
)

if (APPLE) # we need to set correct path for Spot on Mac OS/homebrew
	execute_process(COMMAND brew --prefix
		OUTPUT_VARIABLE BREW_PREFIX
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	message("-- Brew prefix: ${BREW_PREFIX}")
	include_directories("${BREW_PREFIX}/include")
	link_directories("${BREW_PREFIX}/lib")
endif()



add_subdirectory(src)

# Option to build tests
option(BUILD_TESTS "Build the tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
