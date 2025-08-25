cmake_minimum_required(VERSION 3.16)

set(CMAKE_COLOR_MAKEFILE ON)

project(kofola
	LANGUAGES CXX
	VERSION 1.0
)

# Get git hash at compile time
execute_process(
	COMMAND git rev-parse HEAD
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE GIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE
	ERROR_QUIET
)

# If git command failed, set a default value
if(NOT GIT_HASH)
	set(GIT_HASH "unknown")
endif()

# Configure version header
configure_file(
	${CMAKE_SOURCE_DIR}/src/version.hpp.in
	${CMAKE_BINARY_DIR}/src/version.hpp
	@ONLY
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
