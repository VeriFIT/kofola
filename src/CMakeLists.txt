set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# https://dane-bulat.medium.com/cmake-how-to-inspect-and-configure-the-compiler-877e6cb0317f
set(CXX_FLAGS)
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
	# complement_alg_rank.cpp
	complement_alg_rank2.cpp
	ranking.cpp
	util.cpp
)

target_compile_options(kofola PRIVATE ${CXX_FLAGS})

# target_compile_features(kofola PUBLIC cxx_std_17)
find_library(BDDX bddx)
find_library(SPOT spot)

# kofola's libraries
target_link_libraries(kofola
	${SPOT}
	${BDDX}
)
