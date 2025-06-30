This repo is used to store the C++ code that I wrote as part of the fast iterative solvers course.
It is an updated version of the original code modifed with classes and more commentary.

In lib, a lot of useful functions and classes are present.
It contains a lot of matrix utilities in matrixUtils, compressed storage row format of the matrix
in CSRMatrix, modified compressed storage row format of the matrix in MSRMatrix and GMRES routines
in gmresUtils.

The testing of all these utilitiess are done in the tests folder. An example of regular GMRES 
used to solve linear system Ax=b to find x is shown in testGmresUtils.cpp
The matrix A is stored in the file gmres_test_csr.txt in CSR format.

In the future, restarting of GMRES will be added and along with preconditioning. Note that all these
features were done as part of the lectures but still needs to be adopted to utilize the new format.
