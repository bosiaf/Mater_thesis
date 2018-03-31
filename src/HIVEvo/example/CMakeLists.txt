

project(HIVEvo-example) # Hard-code location of the parameter sets into the test by configuring a header file for it.
if(WIN32)
  string(REPLACE "/" "\\\\" current_dir "${current_dir}")
endif()
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
add_executable(HIVEvo-example main.cpp)
target_link_libraries(HIVEvo-example hivevo::HIVEvo)
target_include_directories(HIVEvo-example PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
hivevo_install_target(TARGETS HIVEvo-example EXPORT "")

