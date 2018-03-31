project(HIVEvoTests)
set(TEST_NAME ${PROJECT_NAME}) # Hard-code location of the parameter sets into the test by configuring a header file for it.
if(WIN32)
  string(REPLACE "/" "\\\\" current_dir "${current_dir}")
endif()
include_directories(${CMAKE_CURRENT_BINARY_DIR})
set(test_sources Tests.cpp)
add_executable(${TEST_NAME} ${test_sources})
target_link_libraries(${TEST_NAME} hivevo::HIVEvo )
add_google_test(${TEST_NAME})

