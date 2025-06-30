//Change the parameters in int main to control this code
#include <cassert>
#include "gmresUtils/gmresUtils.hpp"
#include "matrixUtils/matrixUtils.hpp"
#include "MSRMatrix/MSRMatrix.hpp"
#include "CSRMatrix/CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <libgen.h>
#include <math.h>
#include <tuple>
// #include <fstream>
// #include <string>
// #include <iomanip>
// #include <math.h>
// #include <tuple>
// #include <chrono>
// #include <typeinfo>

int testGmresUtils() 
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);
    std::string fileName = sourceDir + "/gmres_test_csr.txt";
    CSRMatrix A{fileName};
    if (!A.getSuccess())
    {
        std::cout<<A.getSuccess()<<" "<<A.getRowIndex().size()<<" "<<A.getColIndex().size()<<" "<<A.getValues().size()<<" ";
        return 1;
    }

    // This is the actual solution vector
    std::vector<double>  xActual (A.getSize(),4) ;
    xActual[0] = 1;
    xActual[A.getSize()-1] = 100;

    // Creating the RHS vector b
    std::vector<double>  b {matrixVecProduct(A,xActual)};
    // Creating the initial guess vector x0 of zeros
    std::vector<double>  x0 (A.getSize(),0) ;
    // Declaring the soluton vector xm
    std::vector<double>  xm;
    // tolerance for the solution vector
    double tol = 1e-6;

    // Doing the GMRES
    std::tuple< std::vector<double>, double> returnTemp;
    returnTemp = GMRES(A, x0, b,1000, tol);
    xm = std::get<0>(returnTemp);
    tol = std::get<1>(returnTemp);

    // Comparing solution vector xm with the actual solution vector xActual
    for (std::size_t i=0; i<xm.size(); i++)
    {
        assert(abs(xm[i] - xActual[i]) <=1e-6);
    }
    // for (std::size_t i=0; i<xm.size(); i++)
    // {
    //     std::cout<<xm[i]<<" ";
    // }
    return 0;
}

int main() {
    testGmresUtils();
    std::cout << "GMRES test passed!" << std::endl;
    return 0;
}