cmake_minimum_required(VERSION 3.1.0)
project(DatasetGen)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Compilation flags: adapt to your needs
if(MSVC)
  # Enable parallel compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  # Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

# Add your project files

file(GLOB header *.h)
file(GLOB source *.cpp)

include_directories(
  # Eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen
)

add_executable(DatasetGen ${source} ${header})
target_link_libraries(${PROJECT_NAME})