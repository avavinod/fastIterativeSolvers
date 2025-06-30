#include "CSRMatrix.hpp"

#include <iostream>
#include <vector>
#include <fstream>

// constructor to input values from file and store it in CSR format
CSRMatrix::CSRMatrix(char const *filename)
{
    
    std::ifstream csrFile;
    csrFile.open (filename,std::ios::in);
    if (!csrFile.is_open())
    {   
        std::cout<<" Error!. CSR file cannot be opened\n";
        success = false;
        csrFile.close();
        return;
    }
    
    // reading data from file
    csrFile.seekg (0,std::ios::beg);
    csrFile>>symmetry>>size>>nonZeros;
    rowIndex.assign(size+1,0);
    colIndex.assign(nonZeros,0);
    values.assign(nonZeros,0);
    // std::size_t tempRowIndex;

    //inputting rowIndex values
    for (std::size_t i=0; i<rowIndex.size(); ++i)
    {
        csrFile >> rowIndex.at(i);
        --rowIndex.at(i);//1 is subtracted to be compatible with C++ index
    }

    //inputting colIndex,values 
    for (std::size_t i=0; i<colIndex.size(); ++i)
    {
        csrFile >> colIndex.at(i) >> values.at(i);
        --colIndex.at(i);//1 is subtracted to be compatible with C++ index
    }

    csrFile.close();
    success = true;
    // std::cout<<success<<" Success!. CSR file finished reading\n";
}

CSRMatrix::CSRMatrix(std::string const filename)
{
    *this = std::move(CSRMatrix(filename.c_str()));
}