#include <numeric>
#include <iostream>
#include <vector>

// g++ -o hello hello.cpp && ./hello

int main(int argc, char **argv)
{
    std::cout << "Hello world!" << std::endl;

    // Basic accesses
    int n = 5;

    std::vector<double> v(10);
    std::iota(v.begin(), v.end(), 0);

    std::cout << n << "-th element: " << v.at(n) << std::endl;

    return 0;
}