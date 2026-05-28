#include <iostream>
#include <memory>
#include <vector>

// g++ -o test_script test.cpp && ./test_script

template <typename T> class Array {
private:
  std::size_t size;
  std::shared_ptr<const T[]> data;

public:
  Array(std::size_t size)
      : size{size}, data(std::shared_ptr<T[]>(new T[size]{})) {}

  void print() {
    std::cout << "Printing the array contents!" << std::endl;
    for (int idx = 0; idx < size; idx++) {
      std::cout << data[idx] << ", ";
    }
    std::cout << std::endl;
  }
};

class CustomResource {
public:
  CustomResource() { std::cout << "Resource acquired!\n"; }
  ~CustomResource() { std::cout << "Resource destroyed and memory freed!\n"; }
  void do_something() const { std::cout << "Resource is being used.\n"; }
};

int main() {

  std::shared_ptr<CustomResource> shared_pointer =
      std::make_shared<CustomResource>();
  shared_pointer->do_something();

  std::unique_ptr<CustomResource> unique_pointer =
      std::make_unique<CustomResource>();
  unique_pointer->do_something();

  std::size_t n{5};

  auto arr = Array<int>(n);

  arr.print();

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
