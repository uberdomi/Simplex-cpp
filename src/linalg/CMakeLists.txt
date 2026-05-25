add_subdirectory(core)
add_subdirectory(ops)

# Umbrella 'linalg_lib' interface library
add_library(linalg_lib INTERFACE)
target_link_libraries(linalg_lib INTERFACE linalg_core linalg_ops)