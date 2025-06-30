#include "matrixUtils.hpp"
#include "CSRMatrix/CSRMatrix.hpp"
#include "MSRMatrix/MSRMatrix.hpp"

#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
#include <math.h>

std::vector <double> operator*(const std::vector <double> &x ,const std::vector <double> &y)
{

    if ( x.size() != y.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between two vectors\n";
        exit (0);
    }

    std::vector<double> output;
    output.resize(x.size());
    for (std::size_t i=0; i<output.size(); i++)
    {
        output[i] = x[i]*y[i];
    }
    return output;
}

std::vector <double> operator*(const std::vector <double> &x ,const double &y)
{

    std::vector<double> output;
    output.resize(x.size());
    for (std::size_t i=0; i<output.size(); i++)
    {
        output[i] = x[i]*y;
    }
    return output;
}

std::vector <double> operator*(const double &x, const std::vector <double> &y)
{

    std::vector<double> output;
    output.resize(y.size());
    for (std::size_t i=0; i<output.size(); i++)
    {
        output[i] = x*y[i];
    }
    return output;
}

std::vector <double> operator+(const std::vector <double> &x ,const std::vector <double> &y) 
{
    if ( x.size() != y.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between two vectors\n";
        exit (0);
    }

    std::vector<double> output;
    output.resize(x.size());
    for (std::size_t i=0; i<output.size(); i++)
    {
        output[i] = x[i] + y[i];
    }
    return output;
}

std::vector <double> operator-(const std::vector <double> &x ,const std::vector <double> &y)
{
    if ( x.size() != y.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between two vectors\n";
        exit (0);
    }

    std::vector<double> output;
    output.resize(x.size());
    for (std::size_t i=0; i<output.size(); i++)
    {
        output[i] = x[i] - y[i];
    }
    return output;
}

double dotProduct(const std::vector <double> &x ,const std::vector <double> &y)
{

    if ( x.size() != y.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between two vectors\n";
        exit (0);
    }
    double dotProduct{0};
    for (std::size_t i=0; i<x.size(); i++)
    {
        dotProduct += x[i]*y[i];
    }
    return dotProduct;
}

double vectorNorm(const std::vector <double> &x)
{

    return sqrt(dotProduct(x,x));
}

std::vector <double> matrixVecProduct(const CSRMatrix &A ,const std::vector <double> &x)
{

    if ( A.getSize() != x.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between matrix and vector\n";
        exit (0);
    }
    std::vector<double> output;
    output.assign(A.getSize(),0);
    
    if (A.getSymmetry() == 's')
    {
        // for symmetric matrix, the bottom half of matrix is not stored.
        // so we can use the fact that A(i,j) = A(j,i) to do two multiplications
        // and addition in one row traversal i.e c(i) += A(i,j)*b(j) and c(j) += A(j,i)*b(i)
        // However this should be doen carefully for the diagonal elements
        // since then c(i) = c(j) and double addition will take place
        std::size_t jc;
        double temp;

        for (std::size_t i=0; i<A.getSize(); i++)
        {
            // c(i) += A(i,j)*b(j)
            jc = A.getRowIndex()[i];
            output[i] +=  A.getValues()[jc]*x[A.getColIndex()[jc]];

            // check if the first element in the row is diagonal
            // if not, do c(j) += A(j,i)*b(i)
            if(i!=A.getColIndex()[jc]) // 
                output[A.getColIndex()[jc]] +=  A.getValues()[jc]*x[i];
            
            for (std::size_t j=A.getRowIndex()[i]+1; j<A.getRowIndex()[i+1]; j++)
            {
                output[i] +=  A.getValues()[j]*x[A.getColIndex()[j]];
                output[A.getColIndex()[j]] +=  A.getValues()[j]*x[i];
            }
        }
    }

    else if (A.getSymmetry() == 'n')
    {
        for (std::size_t i=0; i<A.getSize(); i++)
        {
            for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
                
                output[i] += A.getValues()[j]*x[A.getColIndex()[j]];
        }
            
    }
    else
    {
        std::cout<<"Wrong symmetry character in file\n";
    }
    return output;
}

std::vector <double> matrixVecProduct(const MSRMatrix &A , const std::vector <double> &x)
{

    if ( A.getSize() != x.size() ) 
    {   
        std::cout<<" Error!. Size mismatch between matrix and vector\n";
        exit (0);
    }
    std::vector<double> output;
    output.assign(A.getSize(),0);
    
    if (A.getSymmetry() == 's')
    {
        // NOT IMPLEMENTED
        // // for symmetric matrix, the bottom half of matrix is not stored.
        // // so we can use the fact that A(i,j) = A(j,i) to do two multiplications
        // // and addition in one row traversal i.e c(i) += A(i,j)*b(j) and c(j) += A(j,i)*b(i)
        // // However this should be doen carefully for the diagonal elements
        // // since then c(i) = c(j) and double addition will take place
        // std::size_t jc;
        // double temp;

        // for (std::size_t i=0; i<A.getSize(); i++)
        // {
        //     // c(i) += A(i,j)*b(j)
        //     jc = A.getRowIndex()[i];
        //     output[i] +=  A.getValues()[jc]*x[A.getColIndex()[jc]];

        //     // check if the first element in the row is diagonal
        //     // if not, do c(j) += A(j,i)*b(i)
        //     if(i!=A.getColIndex()[jc]) // 
        //         output[A.getColIndex()[jc]] +=  A.getValues()[jc]*x[i];
            
        //     for (std::size_t j=A.getRowIndex()[i]+1; j<A.getRowIndex()[i+1]; j++)
        //     {
        //         output[i] +=  A.getValues()[j]*x[A.getColIndex()[j]];
        //         output[A.getColIndex()[j]] +=  A.getValues()[j]*x[i];
        //     }
        // }
    }

    else if (A.getSymmetry() == 'n')
    {
        for (std::size_t i=0; i<A.getSize(); i++)
        {
            output[i] += A.getValues()[i]*x[i];

            for (std::size_t j=A.getColIndex()[i]; j<A.getColIndex()[i+1]; j++)
                
                output[i] += A.getValues()[j]*x[A.getColIndex()[j]];
        }
            
    }
    else
    {
        std::cout<<"Wrong symmetry character in file\n";
    }
    return output;
}

MSRMatrix preconditionedMatrixJacobi (CSRMatrix &A)
{
    // All the diagonal elements of the CSR matrix are inverted and stored as values in MSR matrix
    // MSRMatrix M{};
    // std::size_t size{A.getSize()};
    std::vector<double> values;
    values.assign(A.getSize()+1,0);
    std::vector<std::size_t> colIndex ;
    colIndex.assign(A.getSize()+1,A.getSize()+1);

    std::size_t jc;

    if (A.getSymmetry() == 'n')
    {
        for (std::size_t i=0; i<A.getSize(); i++)
        {
            for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
            {
                if(A.getColIndex()[j]==i)
                {
                    values[i]=1.0/A.getValues()[j]; // finding inverse before hand for convenience
                    break;
                }
            }
            std::cout <<std::endl;
        }
    }
    else if (A.getSymmetry() == 's')
    {
        for (std::size_t i=0; i<A.getSize(); i++)
        {
            if(A.getColIndex()[A.getRowIndex()[i]]==i)
            {
                values[i]=1.0/A.getValues()[A.getRowIndex()[i]]; // finding inverse before hand for convenience
            }
        }
    }    
    return MSRMatrix{A.getSymmetry(),A.getSize(),A.getSize(),0,colIndex,values, true};
}

//creates the lower triangular part of the CSR matrix in MSR form for gauss siedel preconditioining
MSRMatrix preconditionedMatrixGaussSiedel (CSRMatrix &A)
{
    //  std::size_t size{A.getSize()};
    std::size_t nonZeros = 0;
    std::size_t zerosDiag = 0;
    std::size_t jm = A.getSize()+2-1;

    std::vector<std::size_t> colIndex ;
    colIndex.assign(A.getSize()+1+A.getNonZeros(),0);    
    std::vector<double> values;
    values.assign(A.getSize()+1+A.getNonZeros(),0);

    std::size_t jc;
    colIndex[0] = jm ;
    if(A.getSymmetry()=='n')
    {
        for (std::size_t i=0; i<A.getSize(); i++)
        {

            for (std::size_t j=A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
            {
                jc = A.getColIndex()[j];
                if(jc==i)
                {
                    values[i]=1.0/A.getValues()[j];// finding inverse before hand for convenience
                    ++nonZeros;
                    break;
                }
                else if(jc<i)
                {
                    colIndex[jm] = jc;
                    values[jm] = A.getValues()[j];
                    ++nonZeros;
                    ++jm;
                }
            }
            colIndex[i+1]=jm;
        }

        zerosDiag=colIndex[A.getSize()]-nonZeros-1;
        colIndex.resize(nonZeros+1+zerosDiag);
        colIndex.shrink_to_fit();
        values.resize(nonZeros+1+zerosDiag);
        values.shrink_to_fit();
    }

    // else if(A.getSymmetry()=='s') // Haven't tested
    // {
    //     vector<std::size_t> temp;
    //     temp.assign(A.getSize(),0);

    //     // First get the values of the diagonals.
    //     // Also
    //     for (std::size_t i=0; i<A.getSize(); i++)
    //     {
    //         int shift = 0;
    //         //diagonals
    //         if(A.getColIndex()[A.getRowIndex()[i]]==i)
    //         {   

    //             shift = 1;
    //             values[i]=1.0/A.getValues()[j];// finding inverse before hand for convenience
    //             ++nonZeros;
    //         }
    // 	    // If there is a diagonal, we start with the next element in the row
    //         for (std::size_t j=A.getRowIndex()[i] + shift; j<A.getRowIndex()[i+1]; j++)
    //         {   
    //             // counting how many elements belong to each row in the lower triangular
    //             // MSR matrix by using the column number from the symmetric CSR matrix.
    //             // Diagonals are excluded. 
    //             temp[A.getColIndex()[j]]++;
    //         }
    //     }
        
    //     // Now we can assign the row pointers of the MSR matrix using temp
    //     temp[0]++;
    //     for (std::size_t i=1; i<=A.getSize(); i++)
    //     {
    //         colIndex[i] = colIndex[i-1] + temp[i-1];
    //         temp[i-1] = 0;
    //     }

    //     // Finally assign the values and the columns of the non-diagonal
    //     // members of the MSR matrix.
    //     for (std::size_t i=0; i<A.getSize(); i++)
    //     {
    //         if(A.getColIndex()[A.getRowIndex()[i]]==i)
    //         {
    //             shift = 1;
    //         }
    //         for (std::size_t j=A.getRowIndex()[i] + shift; j<A.getRowIndex()[i+1]; j++)
    //         {
    //             jm = colIndex[A.getColIndex()[j]] + temp[A.getColIndex()[j]]; 
    //             colIndex[jm] = i;
    //             values[jm] = A.getValues()[j];
    //             temp[A.getColIndex()[j]]++;
    //         }
    //     }
    //     nonZeros = colIndex[A.getSize()-1] - 1;
    //     colIndex.resize(nonZeros);
    //     colIndex.shrink_to_fit();
    //     values.resize(nonZeros);
    //     values.shrink_to_fit();
    // }

    return MSRMatrix{A.getSymmetry(),A.getSize(),A.getSize(),0,colIndex,values, true};
}

//creates the LU matrix in MSR form for ILU preconditioining
MSRMatrix preconditionedMatrixILU (CSRMatrix &A)
{
    // std::size_t nonZeros = A.getNonZeros();
    char symmetry = A.getSymmetry();
    std::size_t size = A.getSize();
    std::size_t nonZeros=A.getNonZeros();
    std::size_t zerosDiag;
    // std::size_t zerosDiag = 0;
    // std::size_t jm = A.getSize()+2-1;

    std::vector<std::size_t> colIndex ;
    colIndex.assign(A.getSize()+1+A.getNonZeros(),0);    
    std::vector<double> values;
    values.assign(A.getSize()+1+A.getNonZeros(),0);
    std::vector<std::size_t> JD;
    JD.assign(A.getSize()+1+A.getNonZeros(),0);
    std::vector<std::size_t> JR;
    JR.assign(A.getSize()+1+A.getNonZeros(),0);

    // symmetry=A.getSymmetry();
    // size=A.getSize();
    // nonZeros=A.getNonZeros();

    long int jc, jk;

    long int jm = A.getSize()+2-1;
    colIndex[0] = jm ;
    for (long int i = 0; i<A.getSize(); i++)
    {
        for (long int j = A.getRowIndex()[i]; j<A.getRowIndex()[i+1]; j++)
        {
            jc = A.getColIndex()[j];
            if(jc == i)
            {
                values[i] = A.getValues()[j];
                JD[i] = jm;
                JR[jc] = i;
            }
            else
            {
                colIndex[jm] = A.getColIndex()[j];
                values[jm] = A.getValues()[j];
                JR[jc] = jm;
                ++jm;
            }

        }
        colIndex[i+1]=jm;
        for( long int j = colIndex[i]; j < JD[i]; ++j )
        {
            jc = colIndex[j];
            values[j] = values[j]/values[jc];
                
            for( long int jj = JD[jc]; jj < colIndex[jc+1]; ++jj )
            {
                jk = JR[colIndex[jj]];
                if(jk!=0)
                {
                    values[jk]-=values[j]*values[jj];
                }
            }

        }
        for( long int j = A.getRowIndex()[i]; j < A.getRowIndex()[i+1]; ++j )
        {
            JR[A.getColIndex()[j]] = 0;
        }

    }
    zerosDiag=colIndex[size]-nonZeros-1;
    colIndex.resize(nonZeros+1+zerosDiag);
    colIndex.shrink_to_fit();
    values.resize(nonZeros+1+zerosDiag);
    values.shrink_to_fit();

    return MSRMatrix{symmetry,size,nonZeros,zerosDiag,colIndex,values, true};
}
