# Types module
set(TYPES_SOURCES
    ranking.cpp
)

set(TYPES_HEADERS
    ranking.hpp
    safra_tree.hpp
    types.hpp
)

# Make the sources available to parent scope
set(TYPES_SOURCES ${TYPES_SOURCES} PARENT_SCOPE)
set(TYPES_HEADERS ${TYPES_HEADERS} PARENT_SCOPE)
