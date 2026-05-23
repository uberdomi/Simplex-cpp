# Define the 'linalg_core' target library
add_library(linalg_core matrix.cpp)

# Any target linking to 'linalg_core' will automatically know to look in this directory for headers
target_include_directories(linalg_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})