#include <iostream>
#include <vector>
#include <fstream>
#include "MSRMatrix.hpp"
#include "CSRMatrix/CSRMatrix.hpp"

MSRMatrix::MSRMatrix( CSRMatrix const& A)
{

    symmetry = A.getSymmetry();
    size = A.getSize();
    nonZeros=A.getNonZeros();
    std::size_t currentColumn = A.getSize()+2-1;
    colIndex.assign(A.getSize()+1+A.getNonZeros(),0);
    values.assign(A.getSize()+1+A.getNonZeros(),0);

    colIndex[0] = currentColumn ;
    for (std::size_t i=0; i<A.getSize(); i++)
    {
        for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
        {
            if (i == A.getColIndex()[j])
            {
                //stores inverted values for convienience
                // values[i]=1.0/A.getValues()[j];
                values[i]=A.getValues()[j];
            }
            else
            {
                colIndex[currentColumn] = A.getColIndex()[j];
                values[currentColumn] = A.getValues()[j];
                ++currentColumn;
            }

        }
        colIndex[i+1]=currentColumn;
    }
    zerosDiag=colIndex[size]-nonZeros-1;
    colIndex.resize(nonZeros+1+zerosDiag);
    colIndex.shrink_to_fit();
    values.resize(nonZeros+1+zerosDiag);
    values.shrink_to_fit();

}