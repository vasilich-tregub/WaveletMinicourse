# CMakeList.txt : CMake project for exercise, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (legall53moments legall53moments.cpp)
add_executable (cdf97moments "cdf97moments.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET legall531D PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
