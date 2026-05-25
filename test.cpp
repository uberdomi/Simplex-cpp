#include <iostream>
#include <vector>

// g++ -o test_script test.cpp && ./test_script

int main() {
  std::vector<int> myVector;

  if (myVector.empty()) {
    std::cout << "The vector is empty." << std::endl;
  }

  myVector.push_back(10);

  if (!myVector.empty()) {
    std::cout << "The vector now has elements." << std::endl;
  }

  return 0;
}
