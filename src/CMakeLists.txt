set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags: common warnings and per-build-type optimizations
set(CXX_FLAGS "")
list(APPEND CXX_FLAGS "-Wall")
list(APPEND CXX_FLAGS "-pedantic-errors")
list(APPEND CXX_FLAGS "-Wextra")
list(APPEND CXX_FLAGS "-Wfloat-equal")
list(APPEND CXX_FLAGS "-fdiagnostics-show-option")
list(APPEND CXX_FLAGS "-Wctor-dtor-privacy")
list(APPEND CXX_FLAGS "-Woverloaded-virtual")
# list(APPEND CXX_FLAGS "-Wold-style-cast")
list(APPEND CXX_FLAGS "-fPIC")
list(APPEND CXX_FLAGS "-fno-strict-aliasing")
list(APPEND CXX_FLAGS "-Weffc++")

# for counting states in inclusion
option(ENABLE_CNT "Enable counter debugging" OFF)

if(ENABLE_CNT)
	add_definitions(-DENABLE_COUNTER)
	set(BINARY_SUFFIX "_cnt") 
else()
	set(BINARY_SUFFIX "") 
endif()
# end for counting states in inclusion

# Include subdirectories
add_subdirectory(algorithms)
add_subdirectory(complement)
add_subdirectory(inclusion)
add_subdirectory(types)
add_subdirectory(util)

# Collect all source files from subdirectories
set(ALL_SOURCES
	main.cpp
)

# Add sources from subdirectories with proper paths
foreach(src ${ALGORITHMS_SOURCES})
	list(APPEND ALL_SOURCES "algorithms/${src}")
endforeach()

foreach(src ${COMPLEMENT_SOURCES})
	list(APPEND ALL_SOURCES "complement/${src}")
endforeach()

foreach(src ${INCLUSION_SOURCES})
	list(APPEND ALL_SOURCES "inclusion/${src}")
endforeach()

foreach(src ${TYPES_SOURCES})
	list(APPEND ALL_SOURCES "types/${src}")
endforeach()

foreach(src ${UTIL_SOURCES})
	list(APPEND ALL_SOURCES "util/${src}")
endforeach()

add_executable(kofola ${ALL_SOURCES})

# Include the directory with the generated version header
target_include_directories(kofola PRIVATE ${CMAKE_BINARY_DIR}/src)

# Include all subdirectories for headers
target_include_directories(kofola PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/algorithms
	${CMAKE_CURRENT_SOURCE_DIR}/complement
	${CMAKE_CURRENT_SOURCE_DIR}/inclusion
	${CMAKE_CURRENT_SOURCE_DIR}/types
	${CMAKE_CURRENT_SOURCE_DIR}/util
)

# Ensure version.hpp is generated before building kofola
add_dependencies(kofola generate_version)

target_compile_options(kofola PRIVATE ${CXX_FLAGS})

# target_compile_features(kofola PUBLIC cxx_std_17)
find_library(BDDX bddx)
find_library(SPOT spot)

# kofola's libraries
target_link_libraries(kofola
	${SPOT}
	${BDDX}
)

set_target_properties(kofola PROPERTIES
	OUTPUT_NAME "kofola${BINARY_SUFFIX}"  # Append _cnt if ENABLE_CNT is ON
)