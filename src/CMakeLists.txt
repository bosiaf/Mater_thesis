set(CMAKE_CXX_STANDARD 14)

set(HIVEVO_VERSION_MAJOR 0)
set(HIVEVO_VERSION_MINOR 1)
set(HIVEVO_VERSION_PATCH 2)
set(HIVEVO_VERSION ${HIVEVO_VERSION_MAJOR}.${HIVEVO_VERSION_MINOR}.${HIVEVO_VERSION_PATCH})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CompilerSpecificSetup)
include(Utils)

add_subdirectory(HIVEvo)

hivevo_install_cmake_target_file(hivevo)
hivevo_install_cmake_config(hivevo)
