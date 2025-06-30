#ifndef _MATRIXUTILS_INCLUDE_
#define _MATRIXUTILS_INCLUDE_

#include "CSRMatrix/CSRMatrix.hpp"
#include "MSRMatrix/MSRMatrix.hpp"
#include <iostream>
#include <vector>
#include <string>

std::vector <double> operator*(const std::vector <double> &x ,const std::vector <double> &y);
std::vector <double> operator*(const std::vector <double> &x ,const double &y);
std::vector <double> operator*(const double &x, const std::vector <double> &y);
std::vector <double> operator+(const std::vector <double> &x ,const std::vector <double> &y);
std::vector <double> operator-(const std::vector <double> &x ,const std::vector <double> &y);
double vectorNorm(const std::vector <double> &x);
double dotProduct(const std::vector <double> &x ,const std::vector <double> &y);
std::vector <double> matrixVecProduct(const CSRMatrix &A ,const std::vector <double> &x);
std::vector <double> matrixVecProduct(const MSRMatrix &A ,const std::vector <double> &x);
MSRMatrix preconditionedMatrixJacobi (CSRMatrix &A);
MSRMatrix preconditionedMatrixGaussSiedel (CSRMatrix &A);
MSRMatrix preconditionedMatrixILU (CSRMatrix &A);
#endif