# Script to generate version.hpp with current git hash at build time

# Set CMake policy to avoid warnings
cmake_policy(SET CMP0053 NEW)

# Get git hash at build time
execute_process(
	COMMAND git rev-parse HEAD
	WORKING_DIRECTORY ${SOURCE_DIR}
	OUTPUT_VARIABLE GIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE
	ERROR_QUIET
)

# If git command failed, set a default value
if(NOT GIT_HASH)
	set(GIT_HASH "unknown")
endif()

# Generate the version header content directly
set(VERSION_CONTENT "#ifndef KOFOLA_VERSION_HPP
#define KOFOLA_VERSION_HPP

#define KOFOLA_GIT_HASH \"${GIT_HASH}\"

#endif // KOFOLA_VERSION_HPP
")

# Check if the file already exists and has the same content
set(VERSION_FILE "${BINARY_DIR}/src/version.hpp")
if(EXISTS "${VERSION_FILE}")
	file(READ "${VERSION_FILE}" EXISTING_CONTENT)
	if("${VERSION_CONTENT}" STREQUAL "${EXISTING_CONTENT}")
		# File content is the same, don't write to avoid unnecessary rebuilds
		return()
	endif()
endif()

# Write the new version file
file(WRITE "${VERSION_FILE}" "${VERSION_CONTENT}")
message(STATUS "Generated version.hpp with git hash: ${GIT_HASH}")
