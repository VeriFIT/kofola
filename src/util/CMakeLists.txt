# Util module
set(UTIL_SOURCES
    util.cpp
    kofola.cpp
)

set(UTIL_HEADERS
    dataflow.hpp
    util.hpp
    kofola.hpp
)

# Make the sources available to parent scope
set(UTIL_SOURCES ${UTIL_SOURCES} PARENT_SCOPE)
set(UTIL_HEADERS ${UTIL_HEADERS} PARENT_SCOPE)
