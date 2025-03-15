#include "Utils.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

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

bool util::isRectangular(const matrix& A){
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

void util::print(const std::vector<int>& v) {
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

// Common matrices

util::matrix util::eye(size_t n) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }

    util::matrix I(n, util::vec(n, 0.0));

    for(int i=0; i<n; i++) {
        I.at(i).at(i) = 1;
    }

    return I;
}

util::matrix util::ones(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    return util::matrix(n, util::vec(m, 1.0));
}

util::matrix util::zeros(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    return util::matrix(n, util::vec(m, 0.0));
}

util::matrix util::t(const matrix& A) {
    util::matrix B(A.at(0).size(), util::vec(A.size(), 0.0));

    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A.at(0).size(); j++) {
            B.at(j).at(i) = A.at(i).at(j);
        }
    }

    return B;
}

void util::append(matrix& A, matrix&& B) {
    A.reserve(A.size() + B.size());
    // Append B to the end of A
    std::transform(B.begin(), B.end(), A.end(), [](util::vec& row) {
        return std::move(row);
    });
}

void util::append(vec& v, vec&& w) {
    v.reserve(v.size() + w.size());
    // Append w to the end of v
    std::transform(w.begin(), w.end(), v.end(), [](double& val) {
        return std::move(val);
    });
}

// Vector operations

util::vec util::add(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    vec result(v.size(), 0.0);

    for(int i=0; i<v.size(); i++) {
        result.at(i) = v.at(i) + w.at(i);
    }

    return result;
}

util::vec util::sub(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    vec result(v.size(), 0.0);

    for(int i=0; i<v.size(); i++) {
        result.at(i) = v.at(i) - w.at(i);
    }

    return result;
}

util::vec util::scale(const vec& v, const double& sf) {
    vec result(v.size(), 0.0);
    result.reserve(v.size());

    for(int i=0; i<v.size(); i++) {
        result.at(i) = v.at(i)*sf;
    }

    return result;
}

double util::dot(const vec& v, const vec& w) {
    // v^t * w

    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    double sum{0.0};
    for(int i=0; i<v.size(); i++) {
        sum += v.at(i)*w.at(i);
    }

    return sum;
}

util::matrix util::dyadic(const vec& v, const vec& w) {
    // v * w^t
    matrix vw(v.size(), vec(w.size(), 0.0));
    for(int i=0; i<v.size(); i++) {
        for(int j=0; j<w.size(); j++) {
            vw.at(i).at(j) = v.at(i)*w.at(j);
        }
    }

    return vw;
}