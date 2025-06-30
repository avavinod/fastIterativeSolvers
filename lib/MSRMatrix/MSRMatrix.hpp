
#ifndef _MSRMATRIX_INCLUDE_
#define _MSRMATRIX_INCLUDE_

#include <iostream>
#include <vector>
#include <string>
#include "CSRMatrix/CSRMatrix.hpp"

// stores the modified sparse row storage MSR matrix (square matrix) https://ece.uwaterloo.ca/~dwharder/aads/Algorithms/Sparse_systems/#:~:text=In%20this%20summary,%20A%20and
class MSRMatrix {
    char symmetry {'n'};
    std::size_t size{0}, nonZeros{0}, zerosDiag{0};
    std::vector<std::size_t> colIndex ;
    std::vector<double> values ;
    bool success {false};

public:
    // constructors
    MSRMatrix() = default;
    MSRMatrix(  
        char _symmetry, 
        std::size_t _size, 
        std::size_t _nonZeros, 
        std::size_t _zerosDiag,
        std::vector<std::size_t> const& _colIndex, 
        std::vector<double> const& _values,
        bool _success = false
    ) : 
        symmetry{_symmetry}, 
        size{_size}, 
        nonZeros{_nonZeros}, 
        zerosDiag{_zerosDiag}, 
        colIndex{_colIndex}, 
        values{_values}, 
        success{_success} {

    }
    // MSRMatrix(std::string const filename); // (not implemented)
    MSRMatrix( CSRMatrix const& A);

    MSRMatrix(MSRMatrix const& other) : 
        symmetry{other.symmetry}, 
        size{other.size}, 
        nonZeros{other.nonZeros}, 
        zerosDiag{other.zerosDiag}, 
        colIndex{other.colIndex}, 
        values{other.values}, 
        success{other.success} {
    }

    MSRMatrix& operator=(MSRMatrix const& other) {
        if (this == &other) return *this;
        symmetry = other.symmetry;
        size = other.size;
        nonZeros = other.nonZeros;
        zerosDiag = other.zerosDiag;
        colIndex = other.colIndex;
        values = other.values;
        success = other.success;
        return *this;
    }
    MSRMatrix(MSRMatrix&& other) : 
        symmetry{std::move(other.symmetry)}, 
        size{std::move(other.size)}, 
        nonZeros{std::move(other.nonZeros)}, 
        zerosDiag{std::move(other.zerosDiag)}, 
        colIndex{std::move(other.colIndex)}, 
        values{std::move(other.values)}, 
        success{std::move(other.success)} {
    }
    MSRMatrix& operator=(MSRMatrix&& other) {
        if (this == &other) return *this;
        symmetry = std::move(other.symmetry);
        size = std::move(other.size);
        nonZeros = std::move(other.nonZeros);
        zerosDiag = std::move(other.zerosDiag);
        colIndex = std::move(other.colIndex);
        values = std::move(other.values);
        success = std::move(other.success);
        return *this;
    }
    ~MSRMatrix() = default;

    // getters
    const char & getSymmetry()                          const { return  symmetry; }
    const std::size_t & getSize()                       const { return  size; }
    const std::size_t & getNonZeros()                   const { return  nonZeros; }
    const std::size_t & getZerosDiag()                     const { return  zerosDiag; }
    const std::vector<std::size_t> & getColIndex()      const { return  colIndex; }
    const std::vector<double> & getValues()             const { return  values; }
    const bool & getSuccess()                           const { return  success; }

    // setters
    void setValues(const std::vector<double> _values)              {values = _values; }
    // auto & getSymmetry()        const { return  symmetry; }
    // auto & getSize()        const { return  size; }
    // auto & getNonZeros()        const { return  nonZeros; }
    // auto & getRowIndex()        const { return  rowIndex; }
    // auto & getColIndex()        const { return  colIndex; }
    // auto & getValues()        const { return  values; }
    // auto & getSuccess()        const { return  success; }
        
    // T      * operator->()       { return  _data; }
    // T const* operator->() const { return  _data; }
    // T      & operator*()        { return *_data; }
    // T      & operator*()  const { return *_data; }    
};

#endif