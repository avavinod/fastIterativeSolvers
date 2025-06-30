//Change the parameters in int main to control this code
#include <cassert>
#include "MSRMatrix/MSRMatrix.hpp"
#include "CSRMatrix/CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <libgen.h>
// #include <fstream>
// #include <string>
// #include <iomanip>
// #include <math.h>
// #include <tuple>
// #include <chrono>
// #include <typeinfo>

int testMSRMatrixRead() 
{

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);
    std::string fileName = sourceDir + "/gmres_test_csr.txt";
    CSRMatrix A{fileName};

    MSRMatrix B{A};

    std::vector<double> v1;
    std::vector<double> v2;

    for (std::size_t i = 0; i < A.getSize(); i++)
    {
        v1.assign(A.getSize(),0);
        v2.assign(A.getSize(),0);

        for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
        {
            v1[A.getColIndex()[j]] = A.getValues()[j];
        }
        // for (std::size_t j=0; j<v1.size(); j++)
        // {
        //     std::cout<<v1[j]<<" ";
        // }
        // std::cout <<std::endl;

        v2[i] = B.getValues()[i];
        for (std::size_t j=B.getColIndex()[i]; j<B.getColIndex()[i+1]; j++)
        {
                v2[B.getColIndex()[j]] = B.getValues()[j];
        }
        // for (std::size_t j=0; j<v2.size(); j++)
        // {
        //     std::cout<<v2[j]<<" ";
        // }
        // std::cout <<std::endl;

        assert(v1==v2);
    }
   
    return 0;        
}

int main() {
    testMSRMatrixRead();
    std::cout << "MSR Matrix Read test passed!" << std::endl;
    return 0;
}