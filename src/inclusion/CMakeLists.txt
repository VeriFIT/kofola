# Inclusion module
set(INCLUSION_SOURCES
    emptiness_check.cpp
    inclusion_check.cpp
)

set(INCLUSION_HEADERS
    emptiness_check.hpp
    inclusion_check.hpp
)

# Make the sources available to parent scope
set(INCLUSION_SOURCES ${INCLUSION_SOURCES} PARENT_SCOPE)
set(INCLUSION_HEADERS ${INCLUSION_HEADERS} PARENT_SCOPE)
