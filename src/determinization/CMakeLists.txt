# Determinization module
set(DETERMINIZATION_SOURCES
    determinize_tela.cpp
)

set(DETERMINIZATION_HEADERS
    determinize_tela.hpp
)

# Make the sources available to parent scope
set(DETERMINIZATION_SOURCES ${DETERMINIZATION_SOURCES} PARENT_SCOPE)
set(DETERMINIZATION_HEADERS ${DETERMINIZATION_HEADERS} PARENT_SCOPE)
