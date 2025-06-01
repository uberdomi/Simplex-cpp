#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "LinAlg.h"
#include <vector>
#include <cmath>

using namespace la;
using Catch::Approx;

// Helper function to check if matrices are approximately equal
bool approxEqual(const matrix& A, const matrix& B, double tolerance = 1e-10) {
    if (A.size() != B.size()) return false;
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i].size() != B[i].size()) return false;
        for (size_t j = 0; j < A[i].size(); ++j) {
            if (std::abs(A[i][j] - B[i][j]) > tolerance) return false;
        }
    }
    return true;
}

// Helper function to check if vectors are approximately equal
bool approxEqual(const vec& v, const vec& w, double tolerance = 1e-10) {
    if (v.size() != w.size()) return false;
    for (size_t i = 0; i < v.size(); ++i) {
        if (std::abs(v[i] - w[i]) > tolerance) return false;
    }
    return true;
}

TEST_CASE("Vector Operations", "[vector]") {
    vec v1 = {1.0, 2.0, 3.0};
    vec v2 = {4.0, 5.0, 6.0};
    
    SECTION("Vector Addition") {
        vec result = add(v1, v2);
        vec expected = {5.0, 7.0, 9.0};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        vec result2 = v1 + v2;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Vector Subtraction") {
        vec result = sub(v2, v1);
        vec expected = {3.0, 3.0, 3.0};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        vec result2 = v2 - v1;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Vector Scaling") {
        vec result = scale(v1, 2.0);
        vec expected = {2.0, 4.0, 6.0};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        vec result2 = v1 * 2.0;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Dot Product") {
        double result = dot(v1, v2);
        double expected = 1.0*4.0 + 2.0*5.0 + 3.0*6.0; // 32.0
        REQUIRE(result == Approx(expected));
        
        // Test operator overload
        double result2 = v1 * v2;
        REQUIRE(result2 == Approx(expected));
    }
    
    SECTION("Element-wise Multiplication") {
        vec result = vecmult(v1, v2);
        vec expected = {4.0, 10.0, 18.0};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        vec result2 = v1 & v2;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Dyadic Product") {
        matrix result = dyadic(v1, v2);
        matrix expected = {{4.0, 5.0, 6.0}, {8.0, 10.0, 12.0}, {12.0, 15.0, 18.0}};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        matrix result2 = v1 % v2;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Subvector") {
        vec v = {1.0, 2.0, 3.0, 4.0, 5.0};
        vec result = subvector(v, 1, 4);
        vec expected = {2.0, 3.0, 4.0};
        REQUIRE(approxEqual(result, expected));
    }
    
    SECTION("Vector Concatenation") {
        vec v = v1;
        conc(v, v2);
        vec expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        REQUIRE(approxEqual(v, expected));
        
        // Test operator overload
        vec result2 = v1 | v2;
        REQUIRE(approxEqual(result2, expected));
    }
}

TEST_CASE("Matrix Operations", "[matrix]") {
    matrix A = {{1.0, 2.0}, {3.0, 4.0}};
    matrix B = {{5.0, 6.0}, {7.0, 8.0}};
    
    SECTION("Matrix Addition") {
        matrix result = add(A, B);
        matrix expected = {{6.0, 8.0}, {10.0, 12.0}};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        matrix result2 = A + B;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix Subtraction") {
        matrix result = sub(B, A);
        matrix expected = {{4.0, 4.0}, {4.0, 4.0}};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        matrix result2 = B - A;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix Scaling") {
        matrix result = scale(A, 2.0);
        matrix expected = {{2.0, 4.0}, {6.0, 8.0}};
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        matrix result2 = A * 2.0;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix-Vector Multiplication") {
        vec v = {1.0, 2.0};
        vec result = mult(A, v);
        vec expected = {5.0, 11.0}; // [1*1+2*2, 3*1+4*2] = [5, 11]
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        vec result2 = A * v;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix-Matrix Multiplication") {
        matrix result = mult(A, B);
        matrix expected = {{19.0, 22.0}, {43.0, 50.0}};
        // [1*5+2*7, 1*6+2*8] = [19, 22]
        // [3*5+4*7, 3*6+4*8] = [43, 50]
        REQUIRE(approxEqual(result, expected));
        
        // Test operator overload
        matrix result2 = A * B;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix Transpose") {
        matrix result = t(A);
        matrix expected = {{1.0, 3.0}, {2.0, 4.0}};
        REQUIRE(approxEqual(result, expected));
    }
    
    SECTION("Identity Matrix") {
        matrix I = eye(3);
        matrix expected = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        REQUIRE(approxEqual(I, expected));
    }
    
    SECTION("Ones Matrix") {
        matrix ones_mat = ones(2, 3);
        matrix expected = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
        REQUIRE(approxEqual(ones_mat, expected));
    }
    
    SECTION("Zeros Matrix") {
        matrix zeros_mat = zeros(2, 3);
        matrix expected = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        REQUIRE(approxEqual(zeros_mat, expected));
    }
    
    SECTION("Matrix Concatenation Horizontal") {
        matrix C = A;
        conc(C, B);
        matrix expected = {{1.0, 2.0, 5.0, 6.0}, {3.0, 4.0, 7.0, 8.0}};
        REQUIRE(approxEqual(C, expected));
        
        // Test operator overload
        matrix result2 = A | B;
        REQUIRE(approxEqual(result2, expected));
    }
    
    SECTION("Matrix Concatenation Vertical") {
        matrix C = A;
        append(C, B);
        matrix expected = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}};
        REQUIRE(approxEqual(C, expected));
        
        // Test operator overload
        matrix result2 = A ^ B;
        REQUIRE(approxEqual(result2, expected));
    }
}

TEST_CASE("Linear System of Equations", "[lse]") {
    SECTION("2x2 System Solution") {
        matrix A = {{2.0, 1.0}, {1.0, 3.0}};
        vec b = {5.0, 6.0};
        
        vec x = solve(A, b);
        
        // Verify solution by substitution
        vec result = A * x;
        REQUIRE(approxEqual(result, b));
        
        // Test operator overload
        vec x2 = A / b;
        REQUIRE(approxEqual(x2, x));
    }
    
    SECTION("3x3 System Solution") {
        matrix A = {{2.0, 1.0, -1.0}, {-3.0, -1.0, 2.0}, {-2.0, 1.0, 2.0}};
        vec b = {8.0, -11.0, -3.0};
        
        vec x = solve(A, b);
        
        // Verify solution by substitution
        vec result = A * x;
        REQUIRE(approxEqual(result, b, 1e-9));
    }
    
    SECTION("Matrix Determinant") {
        matrix A = {{2.0, 1.0}, {1.0, 3.0}};
        double det_A = det(A);
        double expected = 2.0 * 3.0 - 1.0 * 1.0; // 5.0
        REQUIRE(det_A == Approx(expected));
        
        // 3x3 determinant
        matrix B = {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};
        double det_B = det(B);
        double expected_B = 1.0 * (1.0 * 0.0 - 4.0 * 6.0) - 2.0 * (0.0 * 0.0 - 4.0 * 5.0) + 3.0 * (0.0 * 6.0 - 1.0 * 5.0);
        REQUIRE(det_B == Approx(expected_B));
    }
    
    SECTION("Matrix Inversion") {
        matrix A = {{2.0, 1.0}, {1.0, 3.0}};
        matrix A_inv = inv(A);
        
        // Test that A * A_inv = I
        matrix I = A * A_inv;
        matrix expected_I = eye(2);
        REQUIRE(approxEqual(I, expected_I, 1e-10));
        
        // Test operator overload
        matrix A_inv2 = ~A;
        REQUIRE(approxEqual(A_inv2, A_inv));
    }
    
    SECTION("3x3 Matrix Inversion") {
        matrix A = {{2.0, 1.0, 0.0}, {1.0, 2.0, 1.0}, {0.0, 1.0, 2.0}};
        matrix A_inv = inv(A);
        
        // Test that A * A_inv = I
        matrix I = A * A_inv;
        matrix expected_I = eye(3);
        REQUIRE(approxEqual(I, expected_I, 1e-10));
        
        // Test that A_inv * A = I
        matrix I2 = A_inv * A;
        REQUIRE(approxEqual(I2, expected_I, 1e-10));
    }
    
    SECTION("Gaussian Elimination") {
        matrix A = {{2.0, 1.0}, {1.0, 3.0}};
        matrix rhs = {{1.0}, {1.0}};
        
        auto [det_val, solution] = gauss(A, std::move(rhs));
        
        REQUIRE(det_val == Approx(5.0)); // determinant of A
        
        // Verify the solution
        vec b = {1.0, 1.0};
        vec x = matrixAsVec(solution);
        vec result = A * x;
        REQUIRE(approxEqual(result, b, 1e-10));
    }
}

TEST_CASE("Sherman-Morrison Formula", "[sherman-morrison]") {
    SECTION("Basic Sherman-Morrison Update") {
        // Start with a simple invertible matrix
        matrix A = {{2.0, 1.0}, {1.0, 2.0}};
        vec u = {1.0, 0.0};
        vec v = {0.0, 1.0};
        
        // Compute (A + u*v^T)^(-1) using Sherman-Morrison
        matrix result = ShermanMorrison(A, u, v);
        
        // Compute the updated matrix A + u*v^T
        matrix uv_T = dyadic(u, v);
        matrix A_updated = A + uv_T;
        
        // Compute the inverse directly
        matrix A_updated_inv = inv(A_updated);
        
        // Compare results
        REQUIRE(approxEqual(result, A_updated_inv, 1e-10));
    }
    
    SECTION("Sherman-Morrison with Known Inverse") {
        matrix A = {{4.0, 1.0}, {1.0, 3.0}};
        matrix A_inv = inv(A);
        vec u = {1.0, 1.0};
        vec v = {1.0, -1.0};
        
        // Use Sherman-Morrison with precomputed inverse
        vec A_inv_u = A_inv * u;
        vec v_T_A_inv = t(A_inv) * v;
        double denominator = 1.0 + dot(v, A_inv_u);
        
        matrix result = ShermanMorrisonFull(A_inv, A_inv_u, v_T_A_inv, denominator);
        
        // Verify by computing the updated matrix and its inverse
        matrix uv_T = dyadic(u, v);
        matrix A_updated = A + uv_T;
        matrix A_updated_inv = inv(A_updated);
        
        REQUIRE(approxEqual(result, A_updated_inv, 1e-10));
    }
    
    SECTION("Sherman-Morrison Row Update") {
        matrix A = {{3.0, 1.0}, {1.0, 2.0}};
        vec v = {1.0, 1.0};
        
        matrix result = ShermanMorrisonRow(A, v);
        
        // For row update, u is the ones vector
        vec u = {1.0, 1.0};
        matrix uv_T = dyadic(u, v);
        matrix A_updated = A + uv_T;
        matrix A_updated_inv = inv(A_updated);
        
        REQUIRE(approxEqual(result, A_updated_inv, 1e-10));
    }
    
    SECTION("Sherman-Morrison Column Update") {
        matrix A = {{3.0, 1.0}, {1.0, 2.0}};
        vec u = {1.0, 1.0};
        
        matrix result = ShermanMorrisonCol(A, u);
        
        // For column update, v is the ones vector
        vec v = {1.0, 1.0};
        matrix uv_T = dyadic(u, v);
        matrix A_updated = A + uv_T;
        matrix A_updated_inv = inv(A_updated);
        
        REQUIRE(approxEqual(result, A_updated_inv, 1e-10));
    }
    
    SECTION("3x3 Sherman-Morrison") {
        matrix A = {{3.0, 1.0, 0.0}, {1.0, 3.0, 1.0}, {0.0, 1.0, 3.0}};
        vec u = {1.0, 0.0, 1.0};
        vec v = {0.0, 1.0, 0.0};
        
        matrix result = ShermanMorrison(A, u, v);
        
        // Verify by direct computation
        matrix uv_T = dyadic(u, v);
        matrix A_updated = A + uv_T;
        matrix A_updated_inv = inv(A_updated);
        
        REQUIRE(approxEqual(result, A_updated_inv, 1e-10));
    }
}

TEST_CASE("Utility Functions", "[utility]") {
    SECTION("Matrix Equality Check") {
        matrix A = {{1.0, 2.0}, {3.0, 4.0}};
        matrix B = {{1.0, 2.0}, {3.0, 4.0}};
        matrix C = {{1.0, 2.0}, {3.0, 4.1}};
        
        REQUIRE(isEqual(A, B));
        REQUIRE_FALSE(isEqual(A, C));
        REQUIRE(isEqual(A, C, 0.2)); // With larger tolerance
    }
    
    SECTION("Vector Equality Check") {
        vec v1 = {1.0, 2.0, 3.0};
        vec v2 = {1.0, 2.0, 3.0};
        vec v3 = {1.0, 2.0, 3.1};
        
        REQUIRE(isEqual(v1, v2));
        REQUIRE_FALSE(isEqual(v1, v3));
        REQUIRE(isEqual(v1, v3, 0.2)); // With larger tolerance
    }
    
    SECTION("Scalar Equality Check") {
        REQUIRE(isEqual(1.0, 1.0));
        REQUIRE_FALSE(isEqual(1.0, 1.1));
        REQUIRE(isEqual(1.0, 1.1, 0.2)); // With larger tolerance
    }
    
    SECTION("Matrix Rectangularity Check") {
        matrix rect = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
        matrix non_rect = {{1.0, 2.0}, {3.0, 4.0, 5.0}};
        
        REQUIRE(isRectangular(rect));
        REQUIRE_FALSE(isRectangular(non_rect));
    }
    
    SECTION("Matrix-Vector Conversions") {
        vec v = {1.0, 2.0, 3.0};
        matrix m = vecAsMatrix(v);
        vec v_back = matrixAsVec(m);
        
        REQUIRE(approxEqual(v, v_back));
        REQUIRE(m.size() == 3);
        REQUIRE(m[0].size() == 1);
    }
    
    SECTION("Row Swapping") {
        matrix A = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        matrix B = {{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};
        
        swapRows(A, B, 0, 1);
        
        // A[0] should now be B[1] and B[1] should be A[0]
        vec expected_A0 = {9.0, 10.0};
        vec expected_B1 = {1.0, 2.0};
        
        REQUIRE(approxEqual(A[0], expected_A0));
        REQUIRE(approxEqual(B[1], expected_B1));
    }
}

TEST_CASE("Edge Cases and Error Handling", "[edge-cases]") {
    SECTION("Vector Size Mismatch") {
        vec v1 = {1.0, 2.0};
        vec v2 = {1.0, 2.0, 3.0};
        
        REQUIRE_THROWS_AS(add(v1, v2), std::invalid_argument);
        REQUIRE_THROWS_AS(dot(v1, v2), std::invalid_argument);
    }
    
    SECTION("Matrix Size Mismatch") {
        matrix A = {{1.0, 2.0}, {3.0, 4.0}};
        matrix B = {{1.0}, {2.0}};
        
        REQUIRE_THROWS_AS(add(A, B), std::invalid_argument);
    }
    
    SECTION("Invalid Subvector Indices") {
        vec v = {1.0, 2.0, 3.0};
        
        REQUIRE_THROWS_AS(subvector(v, -1, 2), std::invalid_argument);
        REQUIRE_THROWS_AS(subvector(v, 0, 5), std::invalid_argument);
        REQUIRE_THROWS_AS(subvector(v, 2, 1), std::invalid_argument);
    }
    
    SECTION("Singular Matrix Inversion") {
        matrix singular = {{1.0, 2.0}, {2.0, 4.0}}; // Rank 1 matrix
        
        // Should throw or return a matrix with very large values
        // The exact behavior depends on implementation
        REQUIRE_THROWS_AS(inv(singular), std::runtime_error);
    }
    
    SECTION("Zero Size Matrix") {
        REQUIRE_THROWS_AS(eye(0), std::invalid_argument);
        REQUIRE_THROWS_AS(ones(0), std::invalid_argument);
        REQUIRE_THROWS_AS(zeros(0), std::invalid_argument);
    }
}