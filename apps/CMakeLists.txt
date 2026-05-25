# --- Quick testing ---
add_executable(example_app example.cpp)

# Import necessary libraries to the executables
target_link_libraries(example_app PRIVATE linalg_lib)

# --- Main executable ---
add_executable(main_app main.cpp)

# Import necessary libraries to the executables
target_link_libraries(main_app PRIVATE linalg_lib)