# Compiler-specific CMake commands
# They are included from this file in order to avoid unnecessarily making the example CMakeLists.txt too complex.

if(UNIX)
  set(Boost_NO_BOOST_CMAKE On)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif(UNIX)

if(MSVC)
  #Avoid some dumb warnings from Visual studio
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  #Compilation on different cores (faster)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  #For recursive template depth
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _VARIADIC_MAX=10")
endif()

