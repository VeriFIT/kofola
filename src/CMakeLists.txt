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

add_executable(kofola
	main.cpp
	kofola.cpp
	abstract_complement_alg.cpp
	decomposer.cpp
	complement_tela.cpp
	complement_sync.cpp
	complement_alg_mh.cpp
	complement_alg_safra.cpp
	complement_alg_ncsb.cpp
	complement_alg_ncsb_delay.cpp
	# complement_alg_rank.cpp
	complement_alg_rank2.cpp
	complement_alg_init_det.cpp
	complement_alg_subs_tuple.cpp
	ranking.cpp
	util.cpp
	emptiness_check.cpp
	inclusion_check.cpp
	complement_alg_sd_tela.cpp
	)

# Include the directory with the generated version header
target_include_directories(kofola PRIVATE ${CMAKE_BINARY_DIR}/src)

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