# Define the target library, include current directory in it and make it public
add_library(linalg_ops STATIC basic.cpp gauss.cpp matmul.cpp)
target_include_directories(linalg_ops PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link the public other library to this library ~ import it
target_link_libraries(linalg_ops PUBLIC linalg_core)