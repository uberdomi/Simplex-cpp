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
    A.insert(A.end(), B.begin(), B.end());
}

matrix& la::operator^=(matrix& A, const matrix& B) {
    la::append(A,B);
    return A;
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

void la::append(vec& v, const vec& w) {
    v.reserve(v.size() + w.size());
    // Append w to the end of v
    v.insert(v.end(), w.begin(), w.end());
}

vec& la::operator|=(vec& v, const vec& w) {
    la::append(v,w);
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

    vec result(end - start);
    std::copy(v.begin() + start, v.begin() + end, result.begin());

    return result;
}   

// Vector manipulations
void la::scaleVec(vec& v, const double& sf) {
    for(double& val : v) {
        val *= sf;
    }
}

void la::subtractVec(vec& v, const vec& w, const double& sf) {
    if(v.size() != w.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }
    for(int i = 0; i < v.size(); i++) {
        v[i] -= sf * w[i];
    }
}

// --- Matrix operations
// Sums along axes
vec la::rowSum(const matrix& A) {
    if(A.empty()) {
        throw std::invalid_argument("Matrix is empty");
    }
    
    vec result(A.size(), 0.0);
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            result[i] += A[i][j];
        }
    }
    return result;
}

vec la::colSum(const matrix& A) {
    if(A.empty() || A[0].empty()) {
        throw std::invalid_argument("Matrix is empty");
    }
    
    vec result(A[0].size(), 0.0);
    for(const vec& row : A) {
        if(row.size() != result.size()) {
            throw std::invalid_argument("Incompatible dimensions for column sum");
        }
        for(int j = 0; j < row.size(); j++) {
            result[j] += row[j];
        }
    }
    return result;
}

// Basic
matrix la::add(const matrix& A, const matrix& B) {
    if(A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices have incompatible dimensions");
    }
    matrix result(A.size(), vec(A[0].size(), 0.0));
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

matrix la::operator+(const matrix& A, const matrix& B) {
    return la::add(A, B);
}

void la::add_inplace(matrix& A, const matrix& B) {
    if(A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices have incompatible dimensions");
    }
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            A[i][j] += B[i][j];
        }
    }
}

matrix& la::operator+=(matrix& A, const matrix& B) {
    la::add_inplace(A, B);
    return A;
}

matrix la::sub(const matrix& A, const matrix& B) {
    if(A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices have incompatible dimensions");
    }
    matrix result(A.size(), vec(A[0].size(), 0.0));
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

matrix la::operator-(const matrix& A, const matrix& B) {
    return la::sub(A, B);
}

void la::sub_inplace(matrix& A, const matrix& B) {
    if(A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices have incompatible dimensions");
    }
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            A[i][j] -= B[i][j];
        }
    }
}

matrix& la::operator-=(matrix& A, const matrix& B) {
    la::sub_inplace(A, B);
    return A;
}

matrix la::scale(const matrix& A, const double& sf) {
    matrix result = A;
    for(int i = 0; i < result.size(); i++) {
        for(int j = 0; j < result[0].size(); j++) {
            result[i][j] *= sf;
        }
    }
    return result;
}

matrix la::operator*(const matrix& A, const double& sf) {
    return la::scale(A, sf);
}

matrix& la::scale_inplace(matrix& A, const double& sf) {
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            A[i][j] *= sf;
        }
    }
    return A;
}

matrix& la::operator*=(matrix& A, const double& sf) {
    return la::scale_inplace(A, sf);
}

// Multiplications
vec la::mult(const matrix& A, const vec& v) {
    if(A[0].size() != v.size()) {
        throw std::invalid_argument("Matrix and vector dimensions are incompatible");
    }
    vec result(A.size(), 0.0);
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < A[0].size(); j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

vec la::operator*(const matrix& A, const vec& v) {
    return la::mult(A, v);
}

matrix la::mult(const matrix& A, const matrix& B) {
    if(A[0].size() != B.size()) {
        throw std::invalid_argument("Matrix dimensions are incompatible for multiplication");
    }
    matrix result(A.size(), vec(B[0].size(), 0.0));
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j < B[0].size(); j++) {
            for(int k = 0; k < A[0].size(); k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

matrix la::operator*(const matrix& A, const matrix& B) {
    return la::mult(A, B);
}

// Solving LSE's
std::pair<double, matrix> la::gauss(const matrix& LHS, matrix&& RHS) {
    if(LHS.size() != LHS[0].size()) {
        throw std::invalid_argument("LHS matrix must be square for Gaussian elimination");
    }
    
    int n = LHS.size();
    matrix A = LHS;  // Copy LHS
    matrix B = std::move(RHS);  // Move RHS
    double det = 1.0;
    
    // Forward elimination
    for(int i = 0; i < n; i++) {
        // Find pivot
        int pivot = i;
        for(int k = i + 1; k < n; k++) {
            if(std::abs(A[k][i]) > std::abs(A[pivot][i])) {
                pivot = k;
            }
        }
        
        // Swap rows if needed
        if(pivot != i) {
            std::swap(A[i], A[pivot]);
            std::swap(B[i], B[pivot]);
            det *= -1;
        }
        
        // Check for zero pivot
        if(std::abs(A[i][i]) < 1e-12) {
            throw std::runtime_error("Matrix is singular");
        }
        
        det *= A[i][i];
        
        // Eliminate column
        for(int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for(int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            for(int j = 0; j < B[0].size(); j++) {
                B[k][j] -= factor * B[i][j];
            }
        }
    }
    
    // Back substitution
    for(int i = n - 1; i >= 0; i--) {
        for(int j = 0; j < B[0].size(); j++) {
            for(int k = i + 1; k < n; k++) {
                B[i][j] -= A[i][k] * B[k][j];
            }
            B[i][j] /= A[i][i];
        }
    }
    
    return std::make_pair(det, B);
}

std::pair<double, matrix> la::gauss(const matrix& LHS, const vec& RHS) {
    matrix rhs_matrix = la::vecAsMatrix(RHS);
    return la::gauss(LHS, std::move(rhs_matrix));
}

matrix la::inv(const matrix& A) {
    if(A.size() != A[0].size()) {
        throw std::invalid_argument("Matrix must be square for inversion");
    }
    
    matrix I = la::eye(A.size());
    auto result = la::gauss(A, std::move(I));
    return result.second;
}

matrix la::operator~(const matrix& A) {
    return la::inv(A);
}

double la::det(const matrix& A) {
    if(A.size() != A[0].size()) {
        throw std::invalid_argument("Matrix must be square for determinant calculation");
    }
    
    matrix I = la::eye(A.size());
    auto result = la::gauss(A, std::move(I));
    return result.first;
}

vec la::solve(const matrix& A, const vec& v) {
    auto result = la::gauss(A, v);
    return la::matrixAsVec(result.second);
}

vec la::operator/(const matrix& A, const vec& v) {
    return la::solve(A, v);
}

matrix la::ShermanMorrisonFull(const matrix& A_inv, const vec& A_inv_u, const vec& v_T_A_inv, const double denominator) {
    // Sherman-Morrison formula: (A + uv^T)^-1 = A^-1 - (A^-1 u v^T A^-1) / (1 + v^T A^-1 u) = A_inv - (A_inv_u * v_T_A_inv) / denominator

    if(std::abs(denominator) < 1e-12) {
        throw std::runtime_error("Sherman-Morrison formula is not applicable (denominator is zero)");
    }
    
    matrix correction = la::dyadic(A_inv_u, v_T_A_inv);
    la::scale_inplace(correction, 1.0 / denominator);
    
    return A_inv - correction;
}

matrix la::ShermanMorrison(const matrix& A, const vec& u, const vec& v) {
    // Sherman-Morrison formula: (A + uv^T)^-1 = A^-1 - (A^-1 u v^T A^-1) / (1 + v^T A^-1 u)
    matrix A_inv = la::inv(A);
    vec A_inv_u = A_inv * u;
    vec v_T_A_inv = la::t(A_inv) * v;
    double denominator = 1.0 + la::dot(v, A_inv_u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

matrix la::ShermanMorrisonInv(const matrix& A_inv, const vec& u, const vec& v) {
    vec A_inv_u = la::mult(A_inv, u);
    vec v_T_A_inv = la::mult(la::t(A_inv), v);
    double denominator = 1.0 + la::dot(v, A_inv_u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

matrix la::ShermanMorrisonRow(const matrix& A, const vec& v) {
    matrix A_inv = la::inv(A);
    vec A_inv_u = la::colSum(A_inv);
    vec v_T_A_inv = la::mult(la::t(A_inv), v);
    double denominator = 1.0 + la::dot(v, A_inv_u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

matrix la::ShermanMorrisonCol(const matrix& A, const vec& u) {
    matrix A_inv = la::inv(A);
    vec A_inv_u = la::mult(A_inv, u);
    vec v_T_A_inv = la::rowSum(A_inv);
    double denominator = 1.0 + la::dot(v_T_A_inv, u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

matrix la::ShermanMorrisonRowInv(const matrix& A_inv, const vec& v) {
    vec A_inv_u = la::colSum(A_inv);
    vec v_T_A_inv = la::mult(la::t(A_inv), v);
    double denominator = 1.0 + la::dot(v, A_inv_u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

matrix la::ShermanMorrisonColInv(const matrix& A_inv, const vec& u) {
    vec A_inv_u = la::mult(A_inv, u);
    vec v_T_A_inv = la::rowSum(A_inv);
    double denominator = 1.0 + la::dot(v_T_A_inv, u);

    return la::ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
}

// Helper
void la::swapRows(matrix& A, matrix& B, const int& i_A, const int& i_B) {
    if(i_A < 0 || i_A >= A.size() || i_B < 0 || i_B >= B.size()) {
        throw std::invalid_argument("Row indices are out of bounds");
    }
    std::swap(A[i_A], B[i_B]);
}

vec la::getRow(const matrix& A, const int& i) {
    if(i < 0 || i >= A.size()) {
        throw std::invalid_argument("Row index is out of bounds");
    }
    return A[i];
}

vec la::getCol(const matrix& A, const int& i) {
    if(i < 0 || i >= A[0].size()) {
        throw std::invalid_argument("Column index is out of bounds");
    }
    vec col(A.size(), 0.0);
    for(int j = 0; j < A.size(); j++) {
        col[j] = A[j][i];
    }
    return col;
}