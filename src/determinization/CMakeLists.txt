# Determinization module
set(DETERMINIZATION_SOURCES
    abstract_determinize_alg.cpp
    determinize_alg_dac.cpp
    determinize_alg_mh.cpp
    determinize_tela.cpp
)

set(DETERMINIZATION_HEADERS
    abstract_determinize_alg.hpp
    determinize_alg_dac.hpp
    determinize_alg_mh.hpp
    determinize_tela.hpp
)

# Make the sources available to parent scope
set(DETERMINIZATION_SOURCES ${DETERMINIZATION_SOURCES} PARENT_SCOPE)
set(DETERMINIZATION_HEADERS ${DETERMINIZATION_HEADERS} PARENT_SCOPE)
