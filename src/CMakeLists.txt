

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
	complement_alg_rank.cpp
	rankings.cpp
	util.cpp
)

target_compile_features(kofola PUBLIC cxx_std_17)
find_library(BDDX bddx)
find_library(SPOT spot)

# kofola's libraries
target_link_libraries(kofola
	${SPOT}
	${BDDX}
)
