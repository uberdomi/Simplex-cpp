#include "LinAlg.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

using vec = la::vec;
using matrix = la::matrix;

// Comparisons

bool la::isEqual(const double& x, const double& y, const double& tol) {
    return std::abs(x-y) < tol;
}

bool la::isEqual(const vec& v, const vec& w, const double& tol){
    int n = v.size();
    if(w.size() != n) {
        return false;
    }
    for (int i=0; i<n; i++) {
        if (!la::isEqual(v.at(i), w.at(i), tol)) {
            return false;
        }
    }
    return true;
}

bool la::isEqual(const matrix& A, const matrix& B, const double& tol){
    int n = A.size();
    if(B.size() != n) {
        return false;
    }
    for (int i=0; i<n; i++) {
        if (!la::isEqual(A.at(i), B.at(i), tol)) {
            return false;
        }
    }
    return true;
}

bool la::isRectangular(const matrix& A){
    size_t n = A.at(0).size();
    for(const vec& row : A) {
        if(row.size() != n) {
            return false;
        }
    }
    return true;
}

// Printing

void la::print(const vec& v){
    std::cout << "[ ";
    for (const double& val : v) {
        std::cout << val << ", ";
    }
    std::cout << "]" << std::endl;
}

void la::print(const std::vector<int>& v){
    std::cout << "[ ";
    for (const double& val : v) {
        std::cout << val << ", ";
    }
    std::cout << "]" << std::endl;
}

void la::print(const matrix& A){
    std::cout << "{" << std::endl;
    for (const vec& v : A) {
        print(v);
    }
    std::cout << "}" << std::endl;
}

// Conversions

matrix la::vecAsMatrix(const vec& v){
    matrix A(v.size(), vec(1, 0.0));
    for (int i=0; i<v.size(); i++) {
        A.at(i).at(0) = v.at(i);
    }

    return A;
}

vec la::matrixAsVec(const matrix& A){
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

matrix la::eye(size_t n) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }

    matrix I(n, vec(n, 0.0));

    for(int i=0; i<n; i++) {
        I.at(i).at(i) = 1;
    }

    return I;
}

matrix la::ones(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    return matrix(n, vec(m, 1.0));
}

matrix la::zeros(size_t n, size_t m) {
    if(n <= 0){
        throw std::invalid_argument("Matrix size must be positive");
    }
    if(m <= 0) {
        m = n;
    }

    return matrix(n, vec(m, 0.0));
}

matrix la::t(const matrix& A) {
    matrix B(A.at(0).size(), vec(A.size(), 0.0));

    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A.at(0).size(); j++) {
            B.at(j).at(i) = A.at(i).at(j);
        }
    }

    return B;
}

// Glueing matrices
void la::append(matrix& A, const matrix& B) {
    if(!la::isRectangular(B) || !la::isRectangular(A)) {
        throw std::invalid_argument("Matrix not rectangular!");
    }
    if(B.at(0).size() != A.at(0).size()) {
        throw std::invalid_argument("Incompatible dimensions for appending matrices!");
    }

    A.reserve(A.size() + B.size());
    // Append B to the end of A
    std::transform(B.begin(), B.end(), A.end(), [](const vec& row) {
        return row;
    });
}

matrix& la::operator^=(matrix& A, const matrix& B) {
    la::append(A,B);
}

matrix la::operator^(const matrix& A, const matrix& B) {
    matrix C = A;
    C ^= B;
    return C;
}

void la::append(matrix& A, const vec& v) {
    if(v.size() != A.at(0).size()) {
        throw std::invalid_argument("Incompatible dimensions for appending a vector to a matrix!");
    }

    A.push_back(v);
}

matrix& la::operator^=(matrix& A, const vec& v) {
    la::append(A,v);
    return A;
}

matrix la::operator^(const matrix& A, const vec& v) {
    matrix C = A;
    C ^= v;
    return C;
}

void la::conc(vec& v, const vec& w) {
    v.reserve(v.size() + w.size());
    // Append w to the end of v
    std::transform(w.begin(), w.end(), v.end(), [](double& val) {
        return val;
    });
}

vec& la::operator|=(vec& v, const vec& w) {
    la::conc(v,w);
    return v;
}

vec la::operator|(const vec& v, const vec& w) {
    vec x = v;
    x |= w;
    return x;
}

void la::conc(matrix& A,const matrix& B) {
    if(B.size() != A.size()) {
        throw std::invalid_argument("Incompatible dimensions for matrix concatenation!");
    }
    if(!la::isRectangular(B) || !la::isRectangular(A)) {
        throw std::invalid_argument("Matrix not rectangular!");
    }

    for(int i=0; i<A.size(); i++) {
        A.at(i).insert(A.at(i).end(), B.at(i).begin(), B.at(i).end());  // Concatenate in-place
    }
}

matrix& la::operator|=(matrix& A, const matrix& B) {
    la::conc(A,B);
    return A;
}

matrix la::operator|(const matrix& A, const matrix& B) {
    matrix C = A;
    C |= B;
    return C;
}

void la::conc(matrix& A, const vec& v) {
    int i=0;
    for (vec& row : A) {
        row.push_back(v.at(i++));
    }
}

matrix& la::operator|=(matrix& A, const vec& v) {
    la::conc(A,v);
    return A;
}

matrix la::operator|(const matrix& A, const vec& v){
    matrix C = A;
    C |= v;
    return C;
}


// Vector operations
vec la::add(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    vec result(v.size(), 0.0);

    for(int i=0; i<v.size(); i++) {
        result.at(i) = v.at(i) + w.at(i);
    }

    return result;
}

vec la::operator+(const vec& v, const vec& w) {
    return la::add(v,w);
}

void la::add_inplace(vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    for(int i=0; i<v.size(); i++) {
        v.at(i) = v.at(i) + w.at(i);
    }
}

vec& la::operator+=(vec& v, const vec& w){
    la::add_inplace(v,w);
    return v;
}

vec la::sub(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    vec result(v.size(), 0.0);

    for(int i=0; i<v.size(); i++) {
        result.at(i) = v.at(i) - w.at(i);
    }

    return result;
}

vec la::operator-(const vec& v, const vec& w) {
    return la::sub(v,w);
}

void la::sub_inplace(vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    for(int i=0; i<v.size(); i++) {
        v.at(i) = v.at(i) - w.at(i);
    }
}

vec& la::operator-=(vec& v, const vec& w){
    la::sub_inplace(v,w);
    return v;
}

vec la::scale(const vec& v, const double& sf) {
    vec result = v;

    for(int i=0; i<result.size(); i++) {
        result.at(i) *= sf;
    }

    return result;
}

vec la::operator*(const vec& v, const double& sf) {
    return la::scale(v,sf);
}


void la::scale_inplace(vec& v, const double& sf){
    vec result = v;

    for(int i=0; i<v.size(); i++) {
        v.at(i) *= sf;
    }
}

vec& la::operator*=(vec& v, const double& sf){
    la::scale_inplace(v,sf);
    return v;
}

double la::dot(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    double sum{0.0};
    for(int i=0; i<v.size(); i++) {
        sum += v.at(i)*w.at(i);
    }

    return sum;
}
double la::operator*(const vec& v, const vec& w){
    return la::dot(v,w);
}

vec la::vecmult(const vec& v, const vec& w) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    vec result = v;
    for(int i=0; i<v.size(); i++) {
        result.at(i) *= w.at(i);
    }

    return result;
}
vec la::operator&(const vec& v, const vec& w){
    return la::vecmult(v,w);
}

matrix la::dyadic(const vec& v, const vec& w) {
    matrix vw(v.size(), vec(w.size(), 0.0));
    for(int i=0; i<v.size(); i++) {
        for(int j=0; j<w.size(); j++) {
            vw.at(i).at(j) = v.at(i)*w.at(j);
        }
    }

    return vw;
}
matrix la::operator%(const vec& v, const vec& w) {
    return la::dyadic(v,w);
}

vec la::subvector(const vec& v, int start, int end) {
    if(start < 0 || end > v.size()) {
        throw std::invalid_argument("The provided indices are out of range!");
    }

    if(end <= start) {
        throw std::invalid_argument("The start iterator must be before the end!");
    }

    vec result;
    result.reserve(end - start);
    std::transform(v.begin() + start, v.begin() + end, result.begin(), [](const double& val){
        return val;
    });

    return result;
}   

// --- Matrix operations
// Basic
matrix la::add(const matrix& A, const matrix& B);
matrix la::operator+(const matrix& A, const matrix& B);
void la::add_inplace(matrix& A, const matrix& B);
matrix& la::operator+=(matrix& A, const matrix& B);

matrix la::sub(const matrix& A, const matrix& B);
matrix la::operator-(const matrix& A, const matrix& B);
void la::sub_inplace(matrix& A, const matrix& B);
matrix& la::operator-=(matrix& A, const matrix& B);

matrix la::scale(const matrix& A, const double& sf);
matrix la::operator*(const matrix& A, const double& sf);
matrix& la::scale_inplace(matrix& A, const double& sf);
matrix& la::operator*=(matrix& A, const double& sf);

// Multiplications
vec la::mult(const matrix& A, const vec& v);
vec la::operator*(const matrix& A, const vec& v);

matrix la::mult(const matrix& A, const matrix& B);
matrix la::operator*(const matrix& A, const matrix& B);

// Solving LSE's
matrix la::inv(const matrix& A);
matrix la::operator~(const matrix& A);

double la::det(const matrix& A);

vec la::solve(const matrix& A, const vec& v);
vec la::operator/(const matrix& A, const vec& v);

std::pair<double, matrix> la::gauss(const matrix& LHS, matrix&& RHS); // perform the gauss algorithm, getting det. and inv. together
std::pair<double, matrix> la::gauss(const matrix& LHS, const vec& RHS);

matrix la::ShermanMorrison(const matrix& A, const vec& u, const vec& v);

// Helper
void la::swapRows(matrix& A, matrix& B, const int& i_A, const int& i_B);