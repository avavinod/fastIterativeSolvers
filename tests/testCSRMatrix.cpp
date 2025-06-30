//Change the parameters in int main to control this code
#include <cassert>
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

int testCSRMatrixRead() {

    // considers absolute path where the source file is located
    char sourcePath[] = __FILE__;
    std::string sourceDir = dirname(sourcePath);
    std::string fileName = sourceDir + "/gmres_test_csr.txt";
    // std::cout <<  typeid(fileName).name() << " ";
    CSRMatrix A{fileName};

    if (!A.getSuccess())
    {
        std::cout<<A.getSuccess()<<" "<<A.getRowIndex().size()<<" "<<A.getColIndex().size()<<" "<<A.getValues().size()<<" ";
        return 1;
    }

    std::vector<std::size_t> rowIndex = A.getRowIndex();
    for (std::size_t i = 0; i < rowIndex.size(); i++)
    {
        std::cout << rowIndex[i] << " ";
    }
    std::cout << std::endl;

    std::vector<std::size_t> columnIndex = A.getColIndex();
    for (std::size_t i = 0; i < columnIndex.size(); i++)
    {
        std::cout << columnIndex[i] << " ";
    }
    std::cout << std::endl;

    std::vector<double> values = A.getValues();
    for (std::size_t i = 0; i < values.size(); i++)
    {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl;

    assert(A.getValues()[0]==-0.168096667E+05);
    return 0;        

}

int main() {
    testCSRMatrixRead();
    std::cout << "CSR Matrix Read test passed!" << std::endl;
    return 0;
}