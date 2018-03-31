cmake_minimum_required(VERSION 3.5)

project(HIVEvo)

# necessary at least as long as the minimum cmake required version is less than 3.10
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 OLD)
endif()

enable_testing()
add_subdirectory(src)
