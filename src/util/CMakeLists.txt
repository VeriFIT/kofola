# Util module
set(UTIL_SOURCES
    util.cpp
    helpers.cpp
)

set(UTIL_HEADERS
    dataflow.hpp
    util.hpp
    helpers.hpp
)

# Make the sources available to parent scope
set(UTIL_SOURCES ${UTIL_SOURCES} PARENT_SCOPE)
set(UTIL_HEADERS ${UTIL_HEADERS} PARENT_SCOPE)
