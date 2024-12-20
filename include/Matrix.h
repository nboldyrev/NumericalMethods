#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "Type.h"
using MyType = double;
class Matrix {
    private:
        std::vector<std::vector<MyType>> data;
        size_t rows;
        size_t cols;
        MyType eps;
    public:
        Matrix();
        Matrix(const int c,const int r, const MyType e=type<MyType>()(((MyType)0.1)));
        Matrix(const int s,const MyType pr=type<MyType>()(((MyType)0.1)));
        Matrix(const MyType pr);
        Matrix(std::initializer_list<std::vector<MyType>> list, const MyType e=type<MyType>()(((MyType)0.1)));
        Matrix(const Matrix& rhs);
        size_t getRows();
        size_t getCols();
        const size_t getRows() const;
        const size_t getCols()const;
        const MyType getPresision() const;
        void setPresision(const MyType e);

        const MyType norm(const size_t normType);


        MyType& operator()(const int col_ind, const int row_ind);
        const MyType& operator()(const int col_ind, const int row_ind) const;
        MyType& operator()(const int row_ind);
        
        Matrix getCol(const int index);
        const Matrix getCol(const int index) const ;
        Matrix setCol(const int index);
        Matrix getRow(const int index);
        Matrix setRow(const int index,const Matrix& new_row);
        std::pair<size_t, size_t> getMaxPosition(const size_t var,const int par = 0);

        Matrix operator*(const MyType scalar);
        const Matrix operator*(const MyType scalar) const;
        Matrix operator*(const Matrix &rhs);
        Matrix operator+(const Matrix& rhs);
        Matrix operator-(const Matrix &rhs);
        Matrix operator-();

        Matrix addToRow(const size_t rowLhs,const size_t rowRhs, const MyType scalar);
        Matrix addToRow(const size_t rowLhs, const Matrix rowRhs, const MyType scalar);
        Matrix append(const Matrix& col);
        Matrix swapCols(const size_t col1,const size_t col2);
        Matrix swapRows(const size_t row1, const size_t row2);
        Matrix perturb(const MyType perturbationScale);
        Matrix popCol(const size_t col);
        
        Matrix transpose();

        Matrix toUpperTriangleForm();
        Matrix toUpperTriangleForm(std::vector<std::pair<size_t, size_t>>& swaps);
        Matrix toDiagnaleForm(std::vector<std::pair<size_t, size_t>>& swaps);
        Matrix toLowerTriangleForm();

        Matrix inverce();
        Matrix getInverseMatrix();

        friend std::ostream& operator << (std::ostream &os, const Matrix& matrix);
        friend std::istream& operator >>(std::istream &os,  Matrix& matrix);// TODO:: написать чтение из файла

        
};

#endif // MATRIX_H