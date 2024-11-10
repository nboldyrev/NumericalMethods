#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
class Matrix {
    private:
        std::vector<std::vector<double>> dData;;
        std::vector<std::vector<float>> fData;
        size_t rows;
        size_t cols;
        double eps;
    public:
    //* Double
        Matrix();
        Matrix(const int c,const int r, const double e = 2.20E-16);
        Matrix(const double pr);
        Matrix(std::initializer_list<std::vector<double>> list, const double e=2.20E-16);
        Matrix(const Matrix& rhs);


        const size_t getRows() const;
        const size_t getCols()const;
        const double getPresision() const;
        void setPresision(const double e);

        const double norm(const size_t normType);
    



        double& operator()(const int col_ind, const int row_ind);
        const double& operator()(const int col_ind, const int row_ind) const;
        double& operator()(const int row_ind);
        
        Matrix setIdentity(const size_t s);
        Matrix getCol(const int index);
        const Matrix getCol(const int index) const ;
        Matrix setCol(const int index);
        Matrix getRow(const int index);
        Matrix setRow(const int index,const Matrix& new_row);
        std::pair<size_t, size_t> getMaxPosition(const size_t var,const int par = 0);

        Matrix operator*(const double scalar);
        const Matrix operator*(const double scalar) const;
        Matrix operator*(const Matrix &rhs);
        Matrix operator+(const Matrix& rhs);
        Matrix operator-(const Matrix &rhs);
        Matrix operator-();

        Matrix addToRow(const size_t rowLhs,const size_t rowRhs, const double scalar);
        Matrix addToRow(const size_t rowLhs, const Matrix rowRhs, const double scalar);
        Matrix append(const Matrix& col);
        Matrix swapCols(const size_t col1,const size_t col2);
        Matrix swapRows(const size_t row1, const size_t row2);
        Matrix perturb(const double perturbationScale);
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