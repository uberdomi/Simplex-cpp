cmake_minimum_required(VERSION 3.10)
project(Simplex LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directory for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Collect all source files
file(GLOB SRC_FILES src/*.cpp)

# Create a static library from the source files (excluding main.cpp)
list(REMOVE_ITEM SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
add_library(LinAlg STATIC ${SRC_FILES})

# Add the executable target
add_executable(main src/main.cpp)

# Link the library to the executable
target_link_libraries(main PRIVATE LinAlg)
