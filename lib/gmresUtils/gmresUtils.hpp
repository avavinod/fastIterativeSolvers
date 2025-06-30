#ifndef _GMRESUTILS_INCLUDE_
#define _GMRESUTILS_INCLUDE_

#include "CSRMatrix/CSRMatrix.hpp"
#include <iostream>
#include <vector>
#include <tuple>

std::vector <double>  getKrylov( const CSRMatrix& A ,std::vector<std::vector <double>> &orthogonalBasisVectorMatrix,long int no_cols);
std::tuple< std::vector <double>, double> GMRES(const CSRMatrix &A ,std::vector <double> x0, const std::vector <double> &b, int m, double tolerance );

#endif