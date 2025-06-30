//Change the parameters in int main to control this code
#include <cassert>
#include "MSRMatrix/MSRMatrix.hpp"
#include "CSRMatrix/CSRMatrix.hpp"
#include "matrixUtils/matrixUtils.hpp"

#include <iostream>
#include <vector>
#include <libgen.h>
#include <math.h>
// #include <fstream>
// #include <string>
// #include <iomanip>
// #include <math.h>
// #include <tuple>
// #include <chrono>
// #include <typeinfo>

int testVectorUtilities() 
{

    // considers absolute path where the source file is located
    std::vector<double> x = {1,2,1};
    std::vector<double> y = {2,3,5};
    double c = -2;

    std::vector<double> b{x+y};

    assert(b[0]==3);
    assert(b[1]==5);
    assert(b[2]==6);

    b =  x - y;
    assert(b[0]==-1);
    assert(b[1]==-1);
    assert(b[2]==-4);
    
    b =  x * y;
    assert(b[0]==2);
    assert(b[1]==6);
    assert(b[2]==5);

    b =  x * c;
    assert(b[0]==-2);
    assert(b[1]==-4);
    assert(b[2]==-2);

    b =  c * x;
    assert(b[0]==-2);
    assert(b[1]==-4);
    assert(b[2]==-2);

    assert(dotProduct(x,-1*y)==-13);
    assert((vectorNorm(x)-sqrt(6))<1e-6);
    assert(b[1]==-4);
    assert(b[2]==-2);

    return 0;        

}

int testMatrixVecProductCSR()
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);

    std::vector<double> x = {1,2,1};

    // Reading the symmetric CSR matrix from file
    std::string fileName = sourceDir + "/gmres_test_csr_symm.txt";
    CSRMatrix A{fileName};
    assert(A.getSymmetry()=='s');

    //matrix vector product
    std::vector<double> b {matrixVecProduct(A, x)};
    assert(b[0]==13);
    assert(b[1]==11);
    assert(b[2]==5);
    
    // Reading the non-symmetric CSR matrix from file
    fileName = sourceDir + "/gmres_test_csr_non_symm.txt";
    A = CSRMatrix{fileName};
    assert(A.getSymmetry()=='n');

    //matrix vector product
    b = matrixVecProduct(A, x);
    assert(b[0]==13);
    assert(b[1]==8);
    assert(b[2]==0);

    return 0;        

}

int testMatrixVecProductMSR() 
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);

    std::vector<double> x = {1,2,1};

    // Reading the symmetric CSR matrix from file
    std::string fileName = sourceDir + "/gmres_test_csr_symm.txt";
    CSRMatrix A{fileName};

    // Converting the symmetric CSR matrix into MSR matrix *not implemented yet*
    // MSRMatrix B{A};
    // assert(A.getSymmetry()=='s');

    // //matrix vector product
    // std::vector<double> b {matrixVecProduct(A, x)};
    // assert(b[0]==13);
    // assert(b[1]==11);
    // assert(b[2]==5);
    
    // Reading the non-symmetric CSR matrix from file
    fileName = sourceDir + "/gmres_test_csr_non_symm.txt";
    A = CSRMatrix{fileName};
    // B = MSRMatrix{A};
    MSRMatrix B{A};
    assert(A.getSymmetry()=='n');
    
    //matrix vector product
    // b = matrixVecProduct(A, x);
    std::vector<double> b {matrixVecProduct(A, x)};
    assert(b[0]==13);
    assert(b[1]==8);
    assert(b[2]==0);

    return 0;        
}

int testPreconditionedMatrixJacobi() 
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);
    std::string fileName = sourceDir + "/gmres_test_csr_precond_jacobi.txt";
    // std::cout <<  typeid(fileName).name() << " ";
    CSRMatrix A{fileName};

    std::vector<double> x = {2,4,7};

    MSRMatrix B = preconditionedMatrixJacobi(A);

    
    std::vector<double> b {matrixVecProduct(B, x)};
    assert((b[0]-1)<1e-6);
    assert((b[1]-1)<1e-6);
    assert((b[2]-1)<1e-6);
    
    return 0;  
}

int testPreconditionedMatrixGaussSiedel () 
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);
    std::string fileName = sourceDir + "/gmres_test_csr_precond_gauss_siedel.txt";
    // std::cout <<  typeid(fileName).name() << " ";
    CSRMatrix A{fileName};

    std::vector<double> x = {1,1,1};

    MSRMatrix B = preconditionedMatrixGaussSiedel(A);
    
    std::vector<double> b {matrixVecProduct(B, x)};
    assert((b[0]-0.5)<1e-6);
    assert((b[1]-8.25)<1e-6);
    assert((b[2]-3+1/7)<1e-6);

    return 0;  
}

int testPreconditionedMatrixILU () {

    // NOT IMPLEMENTED YET
    // // considers absolute path where the source file is located
    // char sourcePath[] = __FILE__;
    // std::string sourceDir = dirname(sourcePath);
    // std::string fileName = sourceDir + "/gmres_test_csr_precond_gauss_siedel.txt";
    // // std::cout <<  typeid(fileName).name() << " ";
    // CSRMatrix A{fileName};

    // std::vector<double> x = {1,1,1};

    // MSRMatrix B = preconditionedMatrixGaussSiedel(A);
    
    // std::vector<double> b {matrixVecProduct(B, x)};
    // assert((b[0]-0.5)<1e-6);
    // assert((b[1]-8.25)<1e-6);
    // assert((b[2]-3+1/7)<1e-6);

    return 0;  
}

//     //     // std::vector<double> v1;
//     // std::vector<double> v2;
//     // std::cout <<B.getSize()<<std::endl;

//     //         for (std::size_t j=0; j<B.getValues().size(); j++)
//     //     {
//     //         std::cout<<B.getValues()[j]<<" ";
//     //     }
//     //     std::cout <<std::endl;

//     // for (std::size_t i = 0; i < A.getSize(); i++)
//     // {
//     //     // v1.assign(A.getSize(),0);
//     //     v2.assign(A.getSize(),0);

//     //     // for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
//     //     // {
//     //     //     v1[A.getColIndex()[j]] = A.getValues()[j];
//     //     // }
//     //     // for (std::size_t j=0; j<v1.size(); j++)
//     //     // {
//     //     //     std::cout<<v1[j]<<" ";
//     //     // }
//     //     // std::cout <<std::endl;

//     //     v2[i] = B.getValues()[i];
//     //     for (std::size_t j=B.getColIndex()[i]; j<B.getColIndex()[i+1]; j++)
//     //     {
//     //             v2[B.getColIndex()[j]] = B.getValues()[j];
//     //     }
//     //     for (std::size_t j=0; j<v2.size(); j++)
//     //     {
//     //         std::cout<<v2[j]<<" ";
//     //     }
//     //     std::cout <<std::endl;

//     //     // assert(v1==v2);
//     // }

//     // std::vector<double> b {matrixVecProduct(B, x)};
//     // assert((b[0]-1)<1e-6);
//     // assert((b[1]-1)<1e-6);
//     // assert((b[2]-1)<1e-6);
    
//     // fileName = sourceDir + "/gmres_test_csr_non_symm.txt";
//     // A = CSRMatrix{fileName};
//     // assert(A.getSymmetry()=='n');

//     // b = matrixVecProduct(A, x);
//     // assert(b[0]==13);
//     // assert(b[1]==8);
//     // assert(b[2]==0);

//     // if (!A.success)
//     // {
//     //     return 0;
//     // }    
//     // CRSMatrix comp1;
//     // str = 
//     // MSRMatrix B{A};
//     // std::cout << A.getSymmetry();
//     // std::cout << B.getSymmetry();
//     // std::vector<std::size_t> rowIndex = A.getRowIndex();
//     // for (std::size_t i = 0; i < rowIndex.size(); i++)
//     // {
//     //     std::cout << rowIndex[i] << " ";
//     // }

//     // std::vector<std::size_t> columnIndex = A.getColIndex();
//     // for (std::size_t i = 0; i < columnIndex.size(); i++)
//     // {
//     //     std::cout << columnIndex[i] << " ";
//     // }

//     // std::vector<double> values = A.getValues();
//     // for (std::size_t i = 0; i < values.size(); i++)
//     // {
//     //     std::cout << values[i] << " ";
//     // }
//     // assert(A.getValues()[0]==-0.168096667E+05);

//     return 0;        

// }

int main() 
{
    testVectorUtilities();
    std::cout << "Vector Utilites has passed!" << std::endl;
    testMatrixVecProductCSR();
    std::cout << "CSR Matrix product has passed!" << std::endl;
    testMatrixVecProductMSR();
    std::cout << "MSR Matrix product has passed!" << std::endl;
    testPreconditionedMatrixJacobi();
    std::cout << "Preconditioning Matrix Jacobi has passed!" << std::endl;
    testPreconditionedMatrixGaussSiedel();
    std::cout << "Preconditioning Matrix Gauss Siedel has passed!" << std::endl;    
    return 0;
}