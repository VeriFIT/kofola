cmake_minimum_required(VERSION 3.16)

set(CMAKE_COLOR_MAKEFILE ON)

project(kofola
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



add_subdirectory(src)
