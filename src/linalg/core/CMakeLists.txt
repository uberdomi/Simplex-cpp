# Define the target library, include current directory in it and make it public
add_library(linalg_core STATIC matrix.cpp)

target_include_directories(linalg_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})