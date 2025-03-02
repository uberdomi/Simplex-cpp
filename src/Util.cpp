#include "Util.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

// Comparisons

bool util::isEqual(const double& x, const double& y, const double& tol) {
    return std::abs(x-y) < tol;
}

bool util::isEqual(const vec& v, const vec& w, const double& tol){
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

bool util::isEqual(const matrix& A, const matrix& B, const double& tol){
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

bool util::checkRectangular(const matrix& A){
    size_t n = A.at(0).size();
    for(const util::vec& row : A) {
        if(row.size() != n) {
            return false;
        }
    }
    return true;
}

// Vector manipulations

// v * sf
void util::scaleRow(vec& v, const double& sf) {
    for(double& val : v){
        val *= sf;
    }
}

// v - w * sf
void util::subtractRow(vec& v, const vec& w, const double& sf) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for(int i=0; i<v.size(); i++) {
        v.at(i) -= w.at(i)*sf;
    }
}

// Printing

void util::print(const vec& v) {
    std::cout << "[ ";
    for (const double& val : v) {
        std::cout << val << ", ";
    }
    std::cout << "]" << std::endl;
}


void util::print(const matrix& A) {
    std::cout << "{" << std::endl;
    for (const vec& v : A) {
        print(v);
    }
    std::cout << "}" << std::endl;
}

// Conversions

util::matrix util::vecAsMatrix(const vec& v){
    util::matrix A(v.size(), util::vec(1, 0.0));
    for (int i=0; i<v.size(); i++) {
        A.at(i).at(0) = v.at(i);
    }

    return A;
}

util::vec util::matrixAsVec(const matrix& A) {
    vec v(A.size(), 0.0);

    for (int i=0; i<A.size(); i++) {
        if(A.at(i).size() != 1) {
            throw std::invalid_argument("Matrix not convertible to a vector!");
        }
        else {
            v.at(i) = A.at(i).at(0);
        }
    }

    return v;
}