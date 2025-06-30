#include "gmresUtils.hpp"
#include "MSRMatrix/MSRMatrix.hpp"
#include "CSRMatrix/CSRMatrix.hpp"
#include "matrixUtils/matrixUtils.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <tuple>

// Function to return the (noCols + 1)th column of Hessenberg matrix. This column is called h.
// Note that the orthogonalBasisVectorMatrix is a matrix with (noCols + 1) columns, where the first noCols columns are the orthogonal basis vectors 
// and the last column w is the new Krylov subspace basis vector.
// In this function the new Krylov subspace basis vector (noCols + 1)th is updated using a reference.
std::vector <double> getKrylov( const CSRMatrix& A, std::vector<std::vector <double>> &orthogonalBasisVectorMatrix, std::size_t noCols)
{
    // value of lastColIndex
    std::size_t lastColIndex = noCols-1;

    //declaring h at the jth column of orthogonalBasisVectorMatrix matrix
    std::vector <double> h(orthogonalBasisVectorMatrix.size()+1,0);
    double hjp1Inv;

    //finding w using matrix vector product
    std::vector <double> w = matrixVecProduct( A,orthogonalBasisVectorMatrix[lastColIndex]);
    
    for (std::size_t i = 0; i <=lastColIndex; i++)
    { 
        // dot product
        h[i] = dotProduct(orthogonalBasisVectorMatrix[i], w); 
        w = w - h[i]*orthogonalBasisVectorMatrix[i];    
    }

    //finding norm of w
    h[lastColIndex+1] = vectorNorm(w);

    // //finding v_jplus1
    orthogonalBasisVectorMatrix[lastColIndex+1] = w*(1/h[lastColIndex+1]);
    return h;
}

// Function to return the gmres solution vector x0 and the residual norm
// A is the CSRMatrix, x0 is the initial guess, b is the right hand
// m is the number of orthogonal iterations vector to be generated, and tolerance is the convergence criterion.
// If the residual norm is less than or euqal to the tolerance, then m is reduced to the actual number of orthogonal vectors generated, mActual.
// Once the orthognal basis vectors and the upper triangular Hessenberg matrix are generated, we find a solution to the linear
std::tuple< std::vector <double>, double> GMRES(const CSRMatrix& A, std::vector <double> x0, const std::vector <double>& b, int m, double tolerance )
{

    std::vector <double> r0 {b - matrixVecProduct(A,x0)}; 
    double r0Norm {vectorNorm(r0)};
    
    std::vector<std::vector <double>> orthogonalBasisVectorMatrix( (m+1), std::vector<double> (r0.size(), 0)); 
    std::vector<std::vector <double>> hessenbergMatrix(m);

    //assign the first column of orthogonalBasisVectorMatrix with r0
    orthogonalBasisVectorMatrix[0] = r0*(1/r0Norm);
        
    std::vector<double> g(m+1,0);
    g[0] = r0Norm;
    
    std::vector <double> c(m,0),s(m,0);
    double temp1,temp2;
    int mActual = m;
    double rho;

    // Creating the othognal basis vectors and the upper triangular Hessenberg matrix (which is denoted by Rm = upper triangular Hessenberg matrix)
    // The Gram-Schmidt process generates the bais vectors and the Hessenberg matrix which is not uppert triangular yet
    // The Givens rotation stored in c and s is applied to each new column of the Hessenberg matrix to make it upper triangular.
    // Rm is the upper trinagular Hessenberg matrix although it does't exist in the code as a separte matrix, it is represented by the hessenberg matrix vector of vectors.
    for( std::size_t j = 0; j < m; ++j)
    {
        
        hessenbergMatrix[j] = getKrylov(A,orthogonalBasisVectorMatrix,j+1);

        // Applying the previous Givens rotation to the current column of the Hessenberg matrix
        for (std::size_t k = 1; k <= j; ++k)
        {
            temp1 = hessenbergMatrix[j][k-1];
            temp2 = hessenbergMatrix[j][k];
            hessenbergMatrix[j][k-1] =  c[k-1]*temp1 + s[k-1]*temp2 ;
            hessenbergMatrix[j][k]   = -s[k-1]*temp1 + c[k-1]*temp2;
        }

        // Creating the Givens rotation for the current column of the Hessenberg matrix to make it upper triangular
        temp1 = 1/sqrt( hessenbergMatrix[j][j]*hessenbergMatrix[j][j] + hessenbergMatrix[j][j+1]*hessenbergMatrix[j][j+1]); 
        c[j] = hessenbergMatrix[j][j]*temp1; 
        s[j] = hessenbergMatrix[j][j+1]*temp1;

        // Applying the new Givens rotation to the current column of the Hessenberg matrix
        hessenbergMatrix[j][j] = c[j]*hessenbergMatrix[j][j] + s[j]*hessenbergMatrix[j][j+1];
        // hessenbergMatrix[j][j+1] will be 0 with the rotation so no need to set it

        // Applying the new Givens rotation to the g vector
        g[j+1] = -s[j]*g[j];
        g[j] = c[j]*g[j];

        // check if the residual norm is less than or equal to the tolerance in which case we stop the iterations
        rho = abs(g[j+1]);
        if(abs(g[j+1])<=tolerance)
        {
            mActual = j+1;
            break;
        }
        
    }

    //Back substitution to find RInvg = Rm^(-1)gm. Note that Rm is the upper triangular Hessenberg matrix.
    std::vector <double> RInvg{g};
    // This is a loop going in reverse order, so size_t which is unsigned cannot be used for i as that would lead to an infinite loop due to underflow.
    for ( long int i = mActual-1; i>=0; i=i-1) 
    {
        RInvg[i] /= hessenbergMatrix[i][i];
        for(std::size_t j = 0; j<i; ++j)
        {
            RInvg[j] -= hessenbergMatrix[i][j]*RInvg[i];
        }
    }

    //Finding Vm(Rm^(-1)gm) and adding it to x0 which gives the final solution xm although still stored in x0.
    for(std::size_t j = 0; j<mActual; ++j)    
    {
        x0 = x0 + orthogonalBasisVectorMatrix[j]*RInvg[j];
    }

    std::cout<<"No of orthogonal vectors produced = "<< mActual<<" "<<std::endl;
    std::cout<<"Tolerance = "<< rho<<" "<<std::endl;

    return {x0,abs(g[mActual])};

}