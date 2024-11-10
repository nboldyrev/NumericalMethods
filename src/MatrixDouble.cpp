#include "Matrix.h"
#include <iostream>
#include <algorithm>

Matrix::Matrix():dData(),eps(2.20E-16){}
Matrix::Matrix(const int c,const int r, const double e):rows(r),cols(c), eps(e){
    dData.resize(cols,std::vector<double>(rows,0));
}

Matrix::Matrix(const double pr)
:Matrix()
{
(*this).setPresision(pr);
}
Matrix::Matrix(std::initializer_list<std::vector<double>> list, const double e) : dData(list), eps(e)
{
    cols = dData.size();
    rows = dData[0].size();
    for(size_t i = 0; i < cols; ++i){
        for(size_t j = 0; j < rows; ++j) {
            if(std::abs(dData[i][j])< eps)dData[i][j]=0;
        }
    }
}

Matrix::Matrix(const Matrix& rhs)
:dData(rhs.dData), eps(rhs.eps),cols(rhs.cols), rows(rhs.rows)
{

}

const size_t Matrix::getRows() const {
    return this->rows;
}
const size_t Matrix::getCols()const {
    return this->cols;
}

const double Matrix::getPresision() const
{
    return eps;
}

void Matrix::setPresision(const double e)
{
    eps = e;
}

const double Matrix::norm(const size_t normType)
{
    double maxNorm = 0;
    double sum = 0;
    if(normType == 0){//* Кубическая
       if(cols == 1){
           for(size_t i = 0; i < rows; ++i) {
               if(std::abs((*this)(i))>maxNorm)maxNorm = std::abs((*this)(i));
           }
       }
       else {
           for (size_t i = 0; i < rows; ++i) {
               sum = 0;
               auto col = (*this).getRow(i);
               for(size_t j = 0; i < cols; ++i) {
                   sum += std::abs(col(i));
               }
               if(sum>maxNorm)maxNorm=sum;
           }
       }
    }
    if(normType == 1) {//* Октаэдрическая
        if(cols ==1){
                for(size_t i = 0; i < rows; ++i) {
                    maxNorm+=std::abs((*this)(i));
                    if(maxNorm<=eps)maxNorm=0;
                }
        }
        else {
            for(size_t i = 0; i < cols; ++i) {
                sum = 0;
                for(size_t j = 0; j < rows; ++j) {
                    sum+=abs((*this)(i,j));
                    if(sum<=eps)sum = 0;
                }
                if(sum>=maxNorm)maxNorm=sum;
            }
        }
    }
    if(normType == 2) {//* Сферичесская
        if(cols == 1) {
            for(int i = 0; i < rows; ++i) {
                maxNorm+=std::abs(std::pow((*this)(i),2));
            }
            maxNorm = std::sqrt(maxNorm);
        }
        else {//TODO Нужны собственные числа

        }
    }
    return maxNorm;
}

double& Matrix::operator()(const int col_ind, const int row_ind) {
    return dData[col_ind][row_ind];
}
const double& Matrix::operator()(const int col_ind, const int row_ind) const {
    return dData[col_ind][row_ind];
}

double &Matrix::operator()(const int row_ind)
{
   return dData[0][row_ind];
}

Matrix Matrix::getCol(const int index) {
    const auto tmp = dData[index];
    return Matrix({tmp},eps);
}
const Matrix Matrix::getCol(const int index) const {
    const auto tmp = dData[index];
    return Matrix({tmp},eps);
}
Matrix Matrix::getRow(const int index) {
    std::vector<double> row;
    for(int i = 0; i < this->cols; ++i) {
        row.push_back(dData[i][index]);
    }
    return Matrix({row},eps);
}
Matrix Matrix::setRow(const int index,const Matrix& newRow) {
    for(int i = 0; i < this->cols; ++i) {
        (*this)(i,index)=newRow(0,i);
    }
    return (*this);
}

Matrix Matrix::setCol(const int index)
{
    return Matrix();       
}

std::pair<size_t, size_t> Matrix::getMaxPosition(const size_t var, const int par)
{
    auto max_col = var, max_row = var;
    switch (par)
    {
    case 1:
        for(size_t i = var; i < cols; --i) {
            for(size_t j = var; j < cols; --j) {
                if(std::abs((*this)(i,j))>std::abs((*this)(max_col,max_row))) {
                    max_col = i;
                    max_row = j;
                }
                if(j==0) break;
            }
            if(i==0)break;
        }
        break;
    default:
        for(size_t i = var; i < rows; ++i) {
            for(size_t j = var; j < rows; ++j) {
                if(std::abs((*this)(i,j))>std::abs((*this)(max_col,max_row))) {
                    max_col = i;
                    max_row = j;
                }
            }
        }
        break;
    }

    return std::pair<size_t, size_t>(max_col,max_row);
}
  


Matrix Matrix::operator*(const double scalar) {
    Matrix result(*(this));
    for(int i = 0; i < cols; ++i) {
        for(int j = 0; j < rows; ++j) {
            result(i,j)*=scalar;
            if(std::abs(result(i,j))<=eps)result(i,j)=0;
        }
    }
    return result;
}

const Matrix Matrix::operator*(const double scalar) const {
    Matrix result(*(this));
    for(int i = 0; i < cols; ++i) {
        for(int j = 0; j < rows; ++j) {
            result(i,j)*=scalar;
            if(std::abs(result(i,j))<=eps)result(i,j)=0;
        }
    }
    return result;
}
Matrix Matrix::operator*(const Matrix& rhs) {
    Matrix result(rhs.getCols(),rows,eps);
    if (cols != rhs.getRows()) {
        throw std::invalid_argument("Number of columns in the first matrix must equal the number of rows in the second matrix.");
    }
    for (size_t c = 0; c < rhs.getCols(); ++c) { 
        for (size_t r = 0; r < rows; ++r) { 
            for (size_t k = 0; k < cols; ++k) { 
                result(c, r) += (*this)(k,r) * rhs(c,k);
                if(std::abs(result(c,r))<=eps) result(c,r) = 0;
            }
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& rhs) {
    Matrix result(*(this));
    for(int i = 0; i < cols; ++i) {
        for(int j = 0; j < rows; ++j) {
            result(i,j)+=rhs(i,j);
        }
    }
    return result;
}
Matrix Matrix::operator-(const Matrix& rhs) {
    Matrix result(*(this));

    return result+(rhs*(-1.0));
}

Matrix Matrix::operator-()
{
    Matrix result((*this));
    return result*(-1);
}

Matrix Matrix::addToRow(const size_t rowLhs,const size_t rowRhs, const double scalar) {
    ;
    for(size_t i = 0; i < cols; ++i) {
            (*this)(i,rowLhs)+=scalar*(*this)(i,rowRhs);
            if(std::abs((*this)(i,rowLhs))<=eps)(*this)(i,rowLhs)=0;
    }
    return (*this);
}

Matrix Matrix::addToRow(const size_t rowLhs, const Matrix rowRhs, const double scalar)
{
    for(size_t i = 0; i < cols; ++i) {
            (*this)(i,rowLhs)+=scalar*rowRhs(0,i);
            if(std::abs((*this)(i,rowLhs))<=eps)(*this)(i,rowLhs)=0;
    }
    return (*this);
}

Matrix Matrix::transpose() {
    Matrix result(rows,cols,(*this).getPresision());
    for(int i = 0; i < cols; ++i){
        for(int j = 0; j < rows; ++j) {
            (result)(j,i) = dData[i][j];
        }
    }
    return (*this)=result;
}

Matrix Matrix::popCol(const size_t col)
{
    auto Col = (*this).getCol(col);
    dData.erase(dData.begin()+col);
    cols-=1;
    return Col;
}

Matrix Matrix::toUpperTriangleForm()
{
    int maxCol = 0, maxRow = 0;
    for(size_t i = 0; i < rows; ++i) {
        auto [maxCol,maxRow] = (*this).getMaxPosition(i);
        if(std::abs((*this)(maxCol,maxRow))<=eps){
            std::cerr<<"Ошибка"<<std::endl;
        }
        (*this).swapCols(i,maxCol);
        (*this).swapRows(i,maxRow);
        for(size_t j = i+1; j < rows; ++j) {
            auto c = (*this)(i,j)/(*this)(i,i);
            (*this)=(*this).addToRow(j,i,-c);
        }
    }
    return (*this);
}

Matrix Matrix::toLowerTriangleForm()
{
    for(size_t i = cols-1; i > 0; --i) {
        auto [maxCol,maxRow] = (*this).getMaxPosition(i,1);
        if(std::abs((*this)(maxCol,maxRow))<=eps){
            std::cerr<<"Ошибка"<<std::endl;
        }
        (*this).swapCols(i,maxCol);
        (*this).swapRows(i,maxRow);
        for(size_t j = i-1; j >= 0; --j) {
            auto c = (*this)(i,j)/(*this)(i,i);
            (*this)=(*this).addToRow(j,i,-c);
            if(j==0)break;
        }
    }
    return (*this);
}

Matrix Matrix::toUpperTriangleForm( std::vector<std::pair<size_t, size_t>>& swaps)
{
   int maxCol = 0, maxRow = 0;
   for(size_t i = 0; i < rows; ++i) {
        auto [maxCol,maxRow] = (*this).getMaxPosition(i);
        if(std::abs((*this)(maxCol,maxRow))<=eps){
            std::cerr<<"Ошибка"<<std::endl;
        }
        (*this).swapCols(i,maxCol);
        (*this).swapRows(i,maxRow);
        swaps.push_back(std::pair(maxCol,maxRow));
        for(size_t j = i+1; j < rows; ++j) {
            auto c = (*this)(i,j)/(*this)(i,i);
            if(std::abs(c)<=eps)continue;
            (*this)=(*this).addToRow(j,i,-c);
        }
    }
    return (*this);   
}

Matrix Matrix::toDiagnaleForm(std::vector<std::pair<size_t, size_t>> &swaps)
{
    int maxCol = 0, maxRow = 0;
    for(size_t i = 0; i < rows; ++i) {
        auto [maxCol,maxRow] = (*this).getMaxPosition(i);
        if(std::abs((*this)(maxCol,maxRow))<=eps){
            std::cerr<<"Ошибка"<<std::endl;
        }
        (*this).swapCols(i,maxCol);
        (*this).swapRows(i,maxRow);
        swaps.push_back(std::pair(maxCol,maxRow));
        for(size_t j = i+1; j < rows; ++j) {
            auto c = (*this)(i,j)/(*this)(i,i);
            if(std::abs(c)<=eps)continue;
            (*this)=(*this).addToRow(j,i,-c);
        }
        for(int j = i-1;  j >= 0; --j) {
            if(j<0)break;
            auto c = (*this)(i,j)/(*this)(i,i);
            if(std::abs(c)<=eps)continue;
            (*this)=(*this).addToRow(j,i,-c);
        }
    }
    return (*this);  
}

Matrix Matrix::inverce()
{
    std::vector<std::pair<size_t,size_t>>swaps;
    Matrix identityMatrix((int)rows);
    (*this).append(identityMatrix);
    (*this).toUpperTriangleForm(swaps);
    for(size_t i = rows-1; i > 0; --i) {
        for(size_t j = i-1; j >=0; --j) {
            if(std::abs((*this)(i,i))<=eps) {
                size_t newInd=i;
                while(std::abs((*this)(i,newInd))<=eps)newInd--;
                (*this).swapRows(i,newInd);
            }
            auto c = (*this)(i,j)/(*this)(i,i);
            (*this)=(*this).addToRow(j,i,-c);
            if(j == 0) break;
        }
        double tau =1./ (*this)(i,i);
        (*this).setRow(i,((*this).getRow(i))*tau);
    }
    (*this).setRow(0,((*this).getRow(0))*(1./(*this)(0,0)));
//    std::cout<<(*this)<<"\n\n";
    std::reverse(swaps.begin(),swaps.end());
    const auto N = swaps.size();
    for(int i = 0; i < N; ++i) {
       (*this).swapRows(N-i-1,swaps[i].first);
    }
    dData.erase(dData.begin(),dData.end()-rows);
    cols = rows; 
    return (*this);
    
}

Matrix Matrix::getInverseMatrix()
{
    auto tmp(*this);
    return tmp.inverce();
}

Matrix Matrix::append(const Matrix& col) {
    dData.insert(dData.end(),col.dData.begin(),col.dData.end());
    (*this).cols+=col.cols;
    return (*this);
}

Matrix Matrix::swapCols(const size_t col1, const size_t col2)
{
    std::swap(dData[col1],dData[col2]);
    return (*this);
}

Matrix Matrix::swapRows(const size_t row1, const size_t row2)
{
    for(int i = 0; i < cols; ++i) {
        std::swap(dData[i][row1],dData[i][row2]);
    }
    return (*this);
}

Matrix Matrix::perturb(const double perturbationScale)
{
    for(size_t i = 0; i < cols; ++i) {
        for(size_t j = 0; j < rows; ++j) {
            (*this)(i,j)+=std::pow(-1,i+j)*perturbationScale*(*this)(i,j);
        }
    }
    return (*this);
}

std::ostream& operator << (std::ostream &os, const Matrix& matrix)
{
        for(auto i = 0; i < matrix.getRows(); ++i){
            for(auto j = 0; j < matrix.getCols(); ++j) {
                os<<matrix(j,i)<<" ";
            }
            if(i<matrix.getRows()-1)std::cout<< std::endl;
        }
        return os;
}

std::istream &operator>>(std::istream &os, Matrix &matrix)
{
   int N;
   os >> N;
   matrix.dData.resize(N,std::vector<double>(N,0));
   matrix.rows=N;
   matrix.cols = N;
   for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
        os>>matrix(i,j);
    }
   }
   matrix.transpose();
   return os;
}
