# Define the target library, include current directory in it and make it public
add_library(simplex_lib STATIC algorithm.cpp)

target_include_directories(simplex_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})