cmake_minimum_required(VERSION 3.20)
project(SimplexLinAlg VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Help IntelliSense find headers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enable testing globally
enable_testing()

# Orchestrate the subdirectories
add_subdirectory(src/linalg)
add_subdirectory(src/simplex)
add_subdirectory(apps)
add_subdirectory(test)

# --- Deprecated ---

# # Fetch Catch2
# include(FetchContent)
# FetchContent_Declare(
# Catch2
# GIT_REPOSITORY https://github.com/catchorg/Catch2.git
# GIT_TAG v3.4.0 # or any tag you want
# )
# FetchContent_MakeAvailable(Catch2)

# # Add test executable
# add_executable(unit_tests test/unit_testing.cpp)
# target_link_libraries(unit_tests PRIVATE LinAlg Catch2::Catch2WithMain)

# # Enable testing
# enable_testing()
# add_test(NAME LinAlgTests COMMAND unit_tests)
