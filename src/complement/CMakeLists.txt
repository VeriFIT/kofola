# Complement module
set(COMPLEMENT_SOURCES
    complement_sync.cpp
    complement_tela.cpp
    decomposer.cpp
)

set(COMPLEMENT_HEADERS
    complement_sync.hpp
    complement_tela.hpp
    decomposer.hpp
)

# Make the sources available to parent scope
set(COMPLEMENT_SOURCES ${COMPLEMENT_SOURCES} PARENT_SCOPE)
set(COMPLEMENT_HEADERS ${COMPLEMENT_HEADERS} PARENT_SCOPE)
