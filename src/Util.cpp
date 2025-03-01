#include "../include/Util.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

bool util::isEqual(const double& x, const double& y, const double& tol) {
    return std::abs(x-y) < tol;
}

bool util::isEqual(const vec& v, const vec& w, const double& tol=1e-6){
    int n = v.size();
    if(w.size() != n) {
        return false;
    }
    for (int i=0; i<n; i++) {
        if (!util::isEqual(v.at(i), w.at(i), tol)) {
            return false;
        }
    }
    return true;
}

bool util::isEqual(const matrix& A, const matrix& B, const double& tol=1e-6){
    int n = A.size();
    if(B.size() != n) {
        return false;
    }
    for (int i=0; i<n; i++) {
        if (!util::isEqual(A.at(i), B.at(i), tol)) {
            return false;
        }
    }
    return true;
}

// v * sf
void util::scaleRow(vec& v, double sf) {
    for(double& val : v){
        val *= sf;
    }
}

// v - w * sf
void util::subtractRow(vec& v, vec& w, double sf) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for(int i=0; i<v.size(); i++) {
        v.at(i) -= w.at(i)*sf;
    }
}

void util::print(const vec& v) {
    std::cout << "{ ";
    for (const double& val : v) {
        std::cout << val << ", ";
    }
    std::cout << "}" << std::endl;
}


void util::print(const matrix& A) {
    std::cout << "{ ";
    for (const vec& v : A) {
        print(v);
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}