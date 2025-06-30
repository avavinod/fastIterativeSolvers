#ifndef _CSRMATRIX_INCLUDE_
#define _CSRMATRIX_INCLUDE_

#include <iostream>
#include <vector>
#include <string>

// stores the compressed row storage CSR matrix (square)
class CSRMatrix {
    char symmetry {'n'};
    std::size_t size{0}, nonZeros{0};
    std::vector<std::size_t> rowIndex ;
    std::vector<std::size_t> colIndex ;
    std::vector<double> values ;
    bool success {false};

public:
    // constructors
    CSRMatrix() = default;
    CSRMatrix(char const *filename);
    CSRMatrix(std::string const filename);
    CSRMatrix(CSRMatrix const& other) : symmetry{other.symmetry}, size{other.size}, nonZeros{other.nonZeros}, rowIndex{other.rowIndex}, colIndex{other.colIndex}, values{other.values}, success{other.success} {}
    CSRMatrix& operator=(CSRMatrix const& other) {
        if (this == &other) return *this;
        symmetry = other.symmetry;
        size = other.size;
        nonZeros = other.nonZeros;
        rowIndex = other.rowIndex;
        colIndex = other.colIndex;
        values = other.values;
        success = other.success;
        return *this;
    }
    CSRMatrix(CSRMatrix&& other) : symmetry{std::move(other.symmetry)}, size{std::move(other.size)}, nonZeros{std::move(other.nonZeros)}, rowIndex{std::move(other.rowIndex)}, colIndex{std::move(other.colIndex)}, values{std::move(other.values)}, success{std::move(other.success)} {}
    CSRMatrix& operator=(CSRMatrix&& other) {
        if (this == &other) return *this;
        symmetry = std::move(other.symmetry);
        size = std::move(other.size);
        nonZeros = std::move(other.nonZeros);
        rowIndex = std::move(other.rowIndex);
        colIndex = std::move(other.colIndex);
        values = std::move(other.values);
        success = std::move(other.success);
        return *this;
    }
    ~CSRMatrix() = default;

    // getters
    const char & getSymmetry()                       const { return  symmetry; }
    const std::size_t & getSize()                    const { return  size; }
    const std::size_t & getNonZeros()                const { return  nonZeros; }
    const std::vector<std::size_t> & getRowIndex()   const { return  rowIndex; }
    const std::vector<std::size_t> & getColIndex()   const { return  colIndex; }
    const std::vector<double> & getValues()          const { return  values; }
    const bool & getSuccess()                        const { return  success; }


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