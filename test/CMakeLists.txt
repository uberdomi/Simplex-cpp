include(FetchContent)

# Declare GTest source Git URL
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)

# Download and make the targets available
FetchContent_MakeAvailable(googletest)

# Create test executables
add_executable(sample_test sample_test.cpp)

# Link libraries
target_link_libraries(sample_test PRIVATE gtest_main)