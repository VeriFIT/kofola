# Algorithms module
set(ALGORITHMS_SOURCES
    abstract_complement_alg.cpp
    complement_alg_mh.cpp
    complement_alg_ncsb.cpp
    complement_alg_ncsb_delay.cpp
    # complement_alg_rank.cpp  # Temporarily excluded due to compilation issues
    complement_alg_rank2.cpp
    complement_alg_safra.cpp
    complement_alg_sd_inductive.cpp
    complement_alg_sd_tela.cpp
    complement_alg_subs_tuple.cpp
    complement_alg_iadacs.cpp
)

set(ALGORITHMS_HEADERS
    abstract_complement_alg.hpp
    complement_alg_mh.hpp
    complement_alg_ncsb.hpp
    complement_alg_ncsb_delay.hpp
    complement_alg_rank.hpp
    complement_alg_rank2.hpp
    complement_alg_safra.hpp
    complement_alg_sd_inductive.hpp
    complement_alg_sd_tela.hpp
    complement_alg_subs_tuple.hpp
    complement_mh.hpp
    complement_alg_iadacs.hpp
)

# Make the sources available to parent scope
set(ALGORITHMS_SOURCES ${ALGORITHMS_SOURCES} PARENT_SCOPE)
set(ALGORITHMS_HEADERS ${ALGORITHMS_HEADERS} PARENT_SCOPE)
