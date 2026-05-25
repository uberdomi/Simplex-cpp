# --- Quick testing ---
add_executable(example_app example.cpp)

# Import necessary libraries to the executables
target_link_libraries(example_app PRIVATE linalg_lib)
target_include_directories(example_app PUBLIC ${PROJECT_SOURCE_DIR}/src)

# --- Main executable ---
add_executable(main_app main.cpp)

# Import necessary libraries to the executables
target_link_libraries(main_app PRIVATE linalg_lib)