#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>
#include <span>
#include <vector>

// g++ -std=c++23 -o test_script test.cpp && ./test_script

int main() {

  int data[]{1, 2, 3, 6, 8};

  auto data_view = std::span{data};

  int agg = std::accumulate(std::begin(data_view), std::end(data_view), 0);

  std::cout << "Data sum: " << std::to_string(agg) << std::endl;

  return 0;
}
