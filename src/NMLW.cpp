#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <ctime>
using MyType = double;
using matrix = std::vector<std::vector<MyType>>;
using vector = std::vector<MyType>;
namespace DataComputing {

    void readData(const std::string& filename, matrix& matrix, vector& rs) { // Считываем матрицу СЛАУ и правую часть из файла
        std::ifstream file(filename);
        if (file.is_open()) {
            int N;
            file >> N;
            rs.resize(N);
            for(auto i = 0; i < N; ++i) {
                matrix.push_back(std::vector<MyType>(N));
                for(auto j = 0; j < N; ++j) {
                    file >> matrix[i][j]; //TODO: Проверка на вырожденность матрицы СЛАУ
                }
            }
            for(auto i = 0; i < N; ++i) {
                file >> rs[i];
            }
            file.close();
        } else {
            std::cerr << "Ошибка открытия файла!" << std::endl;
        }
    }

        void readData(const std::string& filename, matrix& matrix) { // Считываем матрицу СЛАУ и правую часть из файла
        std::ifstream file(filename);
        if (file.is_open()) {
            int N;
            file >> N;
            for(auto i = 0; i < N; ++i) {
                matrix.push_back(std::vector<MyType>(N));
                for(auto j = 0; j < N; ++j) {
                    file >> matrix[i][j]; //TODO: Проверка на вырожденность матрицы СЛАУ
                }
            }
            file.close();
        } else {
            std::cerr << "Ошибка открытия файла!" << std::endl;
        }
    }    

    void writeVector(const std::string& filename, const vector& X) { // Записываем вектор в файл
        std::ofstream file(filename);
        const auto N = X.size(); 
        if(file.is_open()) {
            for(auto i = 0; i < N; ++i) {
                file << X[i] << " ";
            }
            file << std::endl;
            file.close();
        } else {
            std::cerr << "Ошибка открытия файла!" << std::endl;
        }
    }
    void writeVector( const vector& X) { // Записываем вектор в файл
        const auto N = X.size(); 
       
            for(auto i = 0; i < N; ++i) {
                std::cout << X[i] << " ";
            
            std::cout << std::endl;
            }
    }
    void writeMatrix(const std::string filename, const matrix& matrix) { // Записываем матрицу в файл
        std::ofstream file (filename);
        if(file.is_open()) {
            for(auto i = 0; i < matrix.size(); ++i) {
                for(auto j = 0; j < matrix[i].size(); ++j) {
                    file << matrix[j][i] << " ";
                }
                file << std::endl;
            }
            file << std::endl;
        } else {
            std::cerr << "Ошибка открытия файла!" << std::endl;
        }
    }

    void writeMatrix(const matrix& matrix) { // Записываем матрицу в консоль
        for(auto i = 0; i < matrix.size(); ++i) {
            for(auto j = 0; j < matrix[i].size(); ++j) {
                 std::cout << matrix[i][j] << " ";
             }
             std::cout << std::endl;
         }
         std::cout << std::endl;
    }
}

namespace MatrixOperations {

    matrix matrixAugmentation(const matrix& m_matrix, const vector& rs) { // Создание расширенной матрицы с присоединением правой части
        const auto N = m_matrix.size();
        matrix augmentedMatrix;
        augmentedMatrix.reserve(N);
        vector row;
        for(auto i = 0; i < N; ++i) {
            row.reserve(N+1);
            row.insert(row.begin(), m_matrix[i].begin(), m_matrix[i].end());
            row.push_back(rs[i]);
            augmentedMatrix.push_back(move(row));
            row.clear();
        }
        return augmentedMatrix;
    }

    void swapRows(matrix& m_matrix, const int row1, const int row2) { // Меняем строки местами
        std::swap(m_matrix[row1],m_matrix[row2]);   
    }

    void swapCols(matrix& m_matrix, const int col1, const int col2) { // Меняем столбы местами
        for (auto &row : m_matrix) {
            std::swap(row[col1], row[col2]);
        }   
    }


    matrix getTransposeMatrix(const matrix& m_matrix) {
        const int hight = m_matrix.size();
        const int width = m_matrix[0].size();
        matrix transposeMatrix(width,vector(hight,0));
        for(int i = 0; i < hight; ++i) {
            for(int j = 0; j < width; ++j) {
                    transposeMatrix[j][i] = m_matrix[i][j];
            }
        }
        return transposeMatrix;
    }

    
    
   std::pair<MyType, MyType> getMaxPosition(const matrix& m_matrix, const int var) {
        const auto N = m_matrix.size();
        int max_col = var, max_row = var;
        for (auto i = var; i < N ; ++i) {
			for (auto j = i; j < N; ++j) {
				if (std::abs(m_matrix[i][j]) > std::abs(m_matrix[max_row][max_col])) {
                    max_row = i;
					max_col = j;
				}
			}
		}
		return std::make_pair(max_row, max_col);
    }
    
    void vectorSwap(const std::vector<std::pair<int,int>>& swaps, vector& x) {
        const auto N = swaps.size();
        for(int i = 0; i < N; ++i) {
            std::swap(x[swaps[i].first],x[swaps[i].second]);
        }
    }

    void matrixToUpperTriangleForm(matrix& m_matrix) { // Прямой обход метода Гаусса
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        for (auto i = 0; i < hight; i++) {
            for (auto j = i + 1; j < hight; j++) {
                auto c = m_matrix[j][i] / m_matrix[i][i];
                for (auto k = i; k < width; k++) {
                    m_matrix[j][k] -= c * m_matrix[i][k];
                }
            }
        }
    }

    void matrixToUpperTriangleFormFullChoice(matrix& m_matrix,  std::vector<std::pair<int, int>>& swaps) {
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        int max_col = 0, max_row = 0;
		for (auto i = 0; i < hight-1 ; ++i) {
			auto [max_row, max_col] = MatrixOperations::getMaxPosition(m_matrix,i);
			if (m_matrix[max_row][max_col] == 0) {
			     std::cerr << "Ошибка" << std::endl;
			}
			MatrixOperations::swapRows(m_matrix,i, max_row);
            MatrixOperations::swapCols(m_matrix,i, max_col);
            swaps.push_back(std::pair(i, max_col));
            for(auto j = i + 1; j < hight; ++j) {
                auto c = m_matrix[j][i] / m_matrix[i][i];
                for (auto k = i; k < width; k++) {
                    m_matrix[j][k] -= c * m_matrix[i][k];
                }
            }
        }
    }    

     void matrixToUpperTriangleFormFullChoice(matrix& m_matrix,  std::vector<std::pair<int, int>>& swaps, const MyType eps) {
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        int max_col = 0, max_row = 0;
		for (auto i = 0; i < hight-1 ; ++i) {
			auto [max_row, max_col] = MatrixOperations::getMaxPosition(m_matrix,i);
			if (std::abs(m_matrix[max_row][max_col])<=eps ) {
			     std::cerr << "Ошибка" << std::endl;
			}
			MatrixOperations::swapRows(m_matrix,i, max_row);
            MatrixOperations::swapCols(m_matrix,i, max_col);
            swaps.push_back(std::pair(i, max_col));
            for(auto j = i + 1; j < hight; ++j) {
                auto c = m_matrix[j][i] / m_matrix[i][i];
                if(std::abs(c)<=eps){c=0;continue;}
                for (auto k = i; k < width; k++) {
                    if(std::abs(m_matrix[i][k])<=eps){m_matrix[i][k]=0;continue;}
                    m_matrix[j][k] -= c * m_matrix[i][k];
                    if(std::abs(m_matrix[j][k])<=eps)m_matrix[j][k] = 0;
                }
            }
        }
    }

    void matrixToUpperTriangleFormFullChoice(matrix& m_matrix) {
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        int max_col = 0, max_row = 0;
		for (auto i = 0; i < hight-1 ; ++i) {
			auto [max_row, max_col] = MatrixOperations::getMaxPosition(m_matrix,i);
			if (m_matrix[max_row][max_col] == 0) {
			     std::cerr << "Ошибка" << std::endl;
			}
			MatrixOperations::swapRows(m_matrix,i, max_row);
            MatrixOperations::swapCols(m_matrix,i, max_col);
            for(auto j = i + 1; j < hight; ++j) {
                auto c = m_matrix[j][i] / m_matrix[i][i];
                for (auto k = i; k < width; k++) {
                    m_matrix[j][k] -= c * m_matrix[i][k];
                }
            }
        }
    }  

     void matrixToUpperTriangleFormFullChoice(matrix& m_matrix, const MyType eps) {
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        int max_col = 0, max_row = 0;
		for (auto i = 0; i < hight-1 ; ++i) {
			auto [max_row, max_col] = MatrixOperations::getMaxPosition(m_matrix,i);
			if (m_matrix[max_row][max_col] == 0 ||std::abs(m_matrix[max_row][max_col])<=eps ) {
			     std::cerr << "Ошибка" << std::endl;
			}
			MatrixOperations::swapRows(m_matrix,i, max_row);
            MatrixOperations::swapCols(m_matrix,i, max_col);
            for(auto j = i + 1; j < hight; ++j) {
                auto c = m_matrix[j][i] / m_matrix[i][i];
                if(std::abs(c)<=eps)c=0;
                for (auto k = i; k < width; k++) {
                    if(std::abs(m_matrix[i][k])<=eps){m_matrix[i][k]=0;continue;}
                    m_matrix[j][k] -= c * m_matrix[i][k];
                    if(std::abs(m_matrix[j][k])<=eps)m_matrix[j][k] = 0;
                }
            }
        }
    }

    matrix matrixMultScalar(const matrix& m_matrix, const MyType a) {
            const int hight = m_matrix.size();
            const int width = m_matrix[0].size();
            matrix new_matrix(std::move(m_matrix));
            for(int i = 0; i < hight; ++i) {
                for(int j = 0; j < width; ++j) {
                    new_matrix[i][j]*=a;
                }
            }
        return new_matrix;
    }

    matrix matrixMultScalar(const matrix& m_matrix, const MyType a, const MyType eps) {
            const int hight = m_matrix.size();
            const int width = m_matrix[0].size();
            matrix new_matrix(std::move(m_matrix));
            for(int i = 0; i < hight; ++i) {
                for(int j = 0; j < width; ++j) {
                    new_matrix[i][j]*=a;
                    if(std::abs(new_matrix[i][j])<=eps)new_matrix[i][j]=0;
                }
            }
        return new_matrix;
    }

    matrix matrixPlusMatrix(const matrix& a_matrix, const matrix& b_matrix) {
        if(a_matrix.size()!=b_matrix.size() || a_matrix[0].size()!=b_matrix.size()) {
            std::cerr<<"error";
        }
        matrix c_matrix(a_matrix.size(),vector(a_matrix[0].size(),0));
            const int hight = a_matrix.size();
            const int width = a_matrix[0].size();
            for(int i = 0; i < hight; ++i) {
                for(int j = 0; j < width; ++j) {
                    c_matrix[i][j]=a_matrix[i][j]+b_matrix[i][j];
                }
            }
        return c_matrix;
    }

    matrix matrixPlusMatrix(const matrix& a_matrix, const matrix& b_matrix, const MyType eps) {
        if(a_matrix.size()!=b_matrix.size() || a_matrix[0].size()!=b_matrix.size()) {
            std::cerr<<"error";
        }
        matrix c_matrix(a_matrix.size(),vector(a_matrix[0].size(),0));
            const int hight = a_matrix.size();
            const int width = a_matrix[0].size();
            for(int i = 0; i < hight; ++i) {
                for(int j = 0; j < width; ++j) {
                    c_matrix[i][j]=a_matrix[i][j]+b_matrix[i][j];
                    if(std::abs(c_matrix[i][j])<=eps)c_matrix[i][j]=0;
                }
            }
        return c_matrix;
    }
    
    matrix matrixMultMatrix(const matrix& a_matrix, const matrix& b_matrix) {
        if(a_matrix.size()!=b_matrix[0].size()) {
            std::cerr<<"error";
        }        

        matrix c_matrix(a_matrix.size(),vector(b_matrix[0].size(),0));
        for (int i = 0; i < a_matrix.size(); ++i) {
            for (int j = 0; j < b_matrix[0].size(); ++j) {
                for (int k = 0; k < b_matrix.size() ; ++k) {
                c_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
            }
        }
    }
        return c_matrix;
    }

    matrix matrixMultMatrix(const matrix& a_matrix, const matrix& b_matrix, const MyType eps) {
        if(a_matrix.size()!=b_matrix[0].size()) {
            std::cerr<<"error";
        }        

        matrix c_matrix(a_matrix.size(),vector(b_matrix[0].size(),0));
        for (int i = 0; i < a_matrix.size(); ++i) {
            for (int j = 0; j < b_matrix[0].size(); ++j) {
                for (int k = 0; k < b_matrix.size() ; ++k) {
                if(std::abs(a_matrix[i][k])<=eps || std::abs(b_matrix[k][j])<=eps) {
                    continue;
                }
                c_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
                if (std::abs(c_matrix[i][j])<=eps)c_matrix[i][j]=0;
            }

        }
    }
        return c_matrix;
    }

    vector vectorPlusVector(const vector& a_vector, const vector& b_vector) {
        if(a_vector.size()!=b_vector.size()) {
            std::cerr<<"error";
        }
        vector c_vector(std::move(b_vector));
        const int hight = b_vector.size();
        for(int i = 0; i < hight; ++i) {
                c_vector[i] = a_vector[i] + b_vector[i];
            }
        return c_vector;
    }

    vector vectorPlusVector(const vector& a_vector, const vector& b_vector, const MyType eps) {
        if(a_vector.size()!=b_vector.size()) {
            std::cerr<<"error";
        }
        vector c_vector(std::move(b_vector));
        const int hight = b_vector.size();
        for(int i = 0; i < hight; ++i) {
                c_vector[i] = a_vector[i] + b_vector[i];
                if(std::abs(c_vector[i])<=eps)c_vector[i]=0;
            }
        
        return c_vector;
    }


    vector vectorMultScalar(const vector& a_vector, const MyType val) {
        vector b_vector(std::move(a_vector));
        const int N = b_vector.size();
        for(int i = 0; i < N; ++i) {
            b_vector[i]*=val;
        }
        return b_vector;
    }

    matrix getIdentityMatrix(const int N){
        matrix identityMatrix(N,vector(N, 0));
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                if(i == j)
                    identityMatrix[i][j] = 1;
            }
        } 
        return identityMatrix;
    }

    MyType CubicNorm(const matrix& m_matrix) {
    MyType maxNorm = 0;
    MyType sumNorm = 0;
    for (size_t i = 0; i < m_matrix.size(); ++i) {
        sumNorm = 0;
        for (size_t j = 0; j < m_matrix.size(); ++j) {
            sumNorm +=std::abs(m_matrix[i][j]);
        }
        if (sumNorm > maxNorm) {
            maxNorm = sumNorm;
        }
    }
    return maxNorm;
    }

    MyType CubicNorm(const vector& v_vector) {
    MyType maxNorm = 0;
    for (size_t i = 0; i < v_vector.size(); ++i) {
        MyType sumNorm = std::abs(v_vector[i]);
        if (sumNorm > maxNorm) {
            maxNorm = sumNorm;
        }
    }
    return maxNorm;
    }

    

    MyType getMaxDiagElem(const matrix& m_matrix) {
        MyType maxElem = 0;
        for(int i = 0; i < m_matrix.size(); ++i) {
            if(maxElem < m_matrix[i][i])
                maxElem = m_matrix[i][i];
        }
        return maxElem;
    }

    vector matrixMultVector(const matrix& m_matrix, const vector& m_vec)  {
    vector resVect;
    const int N = m_matrix.size();
    resVect.resize(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            resVect[i] += m_matrix[i][j] * m_vec[j];
        }
    }
    return resVect;
    }

}
namespace LinearAlgebra {

    std::vector<MyType> upperTriangleMatrixElimination(const matrix& m_matrix) { 
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        std::vector<MyType> x(hight,0);
        for (auto i = hight-1; i >=0; i--) {
            x[i]=m_matrix[i][width-1];
            for (auto j = i + 1; j < hight; j++) {
                x[i] -= m_matrix[i][j] * x[j];
            }
            x[i] /= m_matrix[i][i];
            if(i == 0) break;
        }
        return x;
    } 

    std::vector<MyType> upperTriangleMatrixElimination(const matrix& m_matrix,const MyType eps) { 
        const auto hight = m_matrix.size();
        const auto width = m_matrix[1].size();
        std::vector<MyType> x(hight,0);
        for (auto i = hight-1; i >=0; i--) {
            x[i]=m_matrix[i][width-1];
            if(std::abs(x[i])<=eps)x[i]=0;
            for (auto j = i + 1; j < hight; j++) {
                if(std::abs(m_matrix[i][j])<=eps || std::abs(x[j])<=eps)
                    continue;
                x[i] -= m_matrix[i][j] * x[j];
                if(std::abs(x[i])<=eps)x[i]=0;
            }
            x[i] /= m_matrix[i][i];
            if(std::abs(x[i])<=eps)x[i]=0;
            if(i == 0) break;
        }
        return x;
    } 
    
    void gaussForwardElimination(matrix& augmentedMatrix) { // Прямой обход метода Гаусса
        MatrixOperations::matrixToUpperTriangleForm(augmentedMatrix);
    }

    void gaussForwardEliminationFullChoice(matrix& augmentedMatrix,  std::vector<std::pair<int, int>>& swaps) { /*Прямой обход метода Гаусса с полным выбором главного элемента*/
        MatrixOperations::matrixToUpperTriangleFormFullChoice(augmentedMatrix,swaps);
    }
    void gaussForwardEliminationFullChoice(matrix& augmentedMatrix,  std::vector<std::pair<int, int>>& swaps,const MyType eps) { /*Прямой обход метода Гаусса с полным выбором главного элемента*/
        MatrixOperations::matrixToUpperTriangleFormFullChoice(augmentedMatrix,swaps,eps);
    }

    void gaussBackwardElimination(const matrix& augmentedMatrix, vector& x) { // Обратный обход метода Гаусса
        x = upperTriangleMatrixElimination(augmentedMatrix);
    }

    void gaussBackwardElimination(const matrix& augmentedMatrix, vector& x, const MyType eps) { // Обратный обход метода Гаусса
        x = upperTriangleMatrixElimination(augmentedMatrix,eps);
    }

    void gaussSLE(const matrix& equations, const vector& rs, vector& x) { // Решение СЛАУ методом Гаусса
        matrix augmentedEquations = MatrixOperations::matrixAugmentation(equations,rs);
        gaussForwardElimination(augmentedEquations);
        gaussBackwardElimination(augmentedEquations,x);
    }

    void gaussFullChoiceSLE(const matrix& equations, const vector& rs, vector& x) {
        matrix augmentedEquations = MatrixOperations::matrixAugmentation(equations,rs);
        std::vector<std::pair<int, int>> swaps;
        gaussForwardEliminationFullChoice(augmentedEquations,swaps);    
        std::reverse(swaps.begin(), swaps.end()); 
        gaussBackwardElimination(augmentedEquations,x);
        MatrixOperations::vectorSwap(swaps,x);
    }

    void gaussFullChoiceSLE(const matrix& equations, const vector& rs, vector& x, const MyType eps) {
        matrix augmentedEquations = MatrixOperations::matrixAugmentation(equations,rs);
        std::vector<std::pair<int, int>> swaps;
        gaussForwardEliminationFullChoice(augmentedEquations,swaps,eps);    
        std::reverse(swaps.begin(), swaps.end()); 
        gaussBackwardElimination(augmentedEquations,x,eps);
        MatrixOperations::vectorSwap(swaps,x);
    }
    
    matrix getInversMatrix(const matrix& m_matrix) {
        const int N = m_matrix.size();
        matrix inverseMatrix;
        vector x;
        vector b (N,0);
        auto pos = b.begin();
        for(int i = 0; i < N; ++i) {
            b[i] = 1;
            gaussFullChoiceSLE(m_matrix,b,x);
            inverseMatrix.push_back(x);
            b[i] = 0;
        }
        return MatrixOperations::getTransposeMatrix(inverseMatrix);
    }

    matrix getInversMatrix(const matrix& m_matrix, const MyType eps) {
        const int N = m_matrix.size();
        matrix inverseMatrix;
        vector x;
        vector b (N,0);
        auto pos = b.begin();
        for(int i = 0; i < N; ++i) {
            b[i] = 1;
            gaussFullChoiceSLE(m_matrix,b,x,eps);
            inverseMatrix.push_back(x);
            b[i] = 0;
        }
        return MatrixOperations::getTransposeMatrix(inverseMatrix);
    }

    matrix forward_elimination_for_QR_method(matrix& equations) {
	MyType s = 0, c = 0, l = 0;
    const int hight = equations.size();
    const int width = equations[0].size();
    matrix T = MatrixOperations::getIdentityMatrix(hight);
	for (auto i = 0; i < hight - 1; ++i) {
		for (auto j = i + 1; j < hight; ++j) {
			l = std::sqrt(equations[i][i] * equations[i][i] + equations[j][i] * equations[j][i]);
			if (l == 0) {
				continue;
			}
			c = equations[i][i] / l;
			s = equations[j][i] / l;
            auto eqPredI((equations[i]));
            auto eqPredJ((equations[j]));
            auto tPredI((T[i]));
            auto tPredJ((T[j]));
            for(int k = 0; k < width; ++k) {
                equations[i][k] += s * eqPredJ[k] + (c-1) * eqPredI[k];
                equations[j][k] += (c-1) * eqPredJ[k] - s * eqPredI[k] ;
                T[i][k] += s * tPredJ[k] + (c-1) * tPredI[k];
                T[j][k] += (c-1) * tPredJ[k] - s * tPredI[k] ;
            }
		}
	}
    return T;
}

    vector  perturb_rhs(const vector& rhs, MyType max_perturb) {

    vector perturbed_rhs = rhs;
    std::default_random_engine generator(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_real_distribution<MyType> distribution(-max_perturb, max_perturb);

    for (auto& value : perturbed_rhs) {
        MyType perturbation = distribution(generator);
        value += perturbation;
    }

    return perturbed_rhs;
}

    MyType getTau(const matrix& m_matrix) {
        return (MyType)1./MatrixOperations::getMaxDiagElem(m_matrix);
    }  
    
    void fixedPointIteration(const matrix& equations, const vector& rs, vector& x, const MyType tolerance ) {
        MyType tau = 0.5;
        const int N = equations.size();
        matrix identityMatrix = MatrixOperations::getIdentityMatrix((int)equations.size());
        matrix C(N,vector(N,0));
        auto newA = MatrixOperations::matrixMultScalar(equations,tau);
        C = MatrixOperations::matrixMultScalar(MatrixOperations::matrixPlusMatrix(newA,
        MatrixOperations::matrixMultScalar(identityMatrix,-1)),(MyType)-1);
        auto y = MatrixOperations::vectorMultScalar(rs,tau);
        vector x_0(N,10);
        //DataComputing::writeMatrix("A.txt",C);
        x = MatrixOperations::vectorPlusVector(MatrixOperations::matrixMultVector(C,x_0),y);
        //DataComputing::writeVector("x-x_0.txt",MatrixOperations::vectorPlusVector(x, MatrixOperations::vectorMultScalar(x_0,(MyType)-1)));
        auto normC = MatrixOperations::CubicNorm(C);
        //std::cout<<normC<<std::endl;
        auto k = std::abs(((1-normC)/normC)*tolerance);
        auto Norm = MatrixOperations::CubicNorm(MatrixOperations::vectorPlusVector(x, MatrixOperations::vectorMultScalar(x_0,(MyType)-1)));
        while(MatrixOperations::CubicNorm(MatrixOperations::vectorPlusVector(x, MatrixOperations::vectorMultScalar(x_0,(MyType)-1)))>k){
            x_0 = x;
            x = MatrixOperations::vectorPlusVector(MatrixOperations::matrixMultVector(C,x),y);
            Norm = MatrixOperations::CubicNorm(MatrixOperations::vectorPlusVector(x, MatrixOperations::vectorMultScalar(x_0,(MyType)-1)));
        }
    }
}

namespace Lab1 {//* NOTE::  файлы 
    const std::string inDir = "../tests/L1/";
    const std::string outDir = "../results/L1/";

    const std::string fileIn1 = inDir+"test1.txt";
    const std::string fileIn2 = inDir+"test2.txt";
    const std::string fileIn3 = inDir+"test3.txt";
    const std::string fileIn4 = inDir+"test4.txt";
    const std::string fileIn5 = inDir+"test5.txt";
    const std::string dopFileIn1 = inDir+"dop_test1.txt";
    const std::string dopFileIn2 = inDir+"dop_test2.txt";

    const std::string fileOut1 = outDir+"result1.txt";
    const std::string fileOut2 = outDir+"result2.txt";
    const std::string fileOut3 = outDir+"result3.txt";
    const std::string fileOut4 = outDir+"result4.txt";
    const std::string fileOut5 = outDir+"result5.txt";
    const std::string dopFileOut1 = outDir+"result6.txt";
    const std::string dopFileOut2=outDir+"result7.txt";
}

namespace Lab2 {//* NOTE:: файлы
    const std::string inDir = "../tests/L2/";
    const std::string outDir = "../results/L2/";

    const std::string fileIn1 = inDir+"test1.txt";
    const std::string fileIn2 = inDir+"test2.txt";
    const std::string fileIn3 = inDir+"test3.txt";
    const std::string fileIn4 = inDir+"test4.txt";
    const std::string fileIn5 = inDir+"test5.txt";
    const std::string dopFileIn1 = inDir+"dop_test1.txt";
    const std::string dopFileIn2 = inDir+"dop_test2.txt";

    const std::string fileOut1 = outDir+"result1.txt";
    const std::string fileOut2 = outDir+"result2.txt";
    const std::string fileOut3 = outDir+"result3.txt";
    const std::string fileOut4 = outDir+"result4.txt";
    const std::string fileOut5 = outDir+"result5.txt";
    const std::string dopFileOut1 = outDir+"result6.txt";
    const std::string dopFileOut2=outDir+"result7.txt";
}

namespace Lab1 {
    void gaussNoChoiseDemo() {
        matrix equations;
        vector rightSide;
        vector x;
        DataComputing::readData(fileIn1, equations, rightSide);
        LinearAlgebra::gaussSLE(equations, rightSide, x);
        DataComputing::writeVector(fileOut1,x);
    }

    void gaussFullChoiseDemo() {
        matrix equations;
        vector rightSide;
        vector x;
        auto eps = 1.0E-20;
        DataComputing::readData(fileIn5, equations, rightSide);
        LinearAlgebra::gaussFullChoiceSLE(equations, rightSide, x,eps);
        DataComputing::writeVector(fileOut5,x);
    }

    void QRdemo(const std::string& filenameIn, const std::string& filenameOut) {
        matrix equations;
        vector rightSide;
        vector x;
        DataComputing::readData(filenameIn,equations,rightSide);
        matrix R((equations));
        matrix T = LinearAlgebra::forward_elimination_for_QR_method(R);
        std::cout<<"Q: "<<std::endl;
        DataComputing::writeMatrix(MatrixOperations::getTransposeMatrix(T));
        std::cout<<"R: "<<std::endl;
        DataComputing::writeMatrix(R);
        vector newRs = MatrixOperations::matrixMultVector(T,rightSide);
        matrix newEq = MatrixOperations::matrixAugmentation(R,newRs);
        x = LinearAlgebra::upperTriangleMatrixElimination(newEq);
        DataComputing::writeVector(filenameOut,x);
    }
    void inverseMatirxDemo() {
        matrix m_matrix;
        DataComputing::readData(fileIn5,m_matrix);
        auto eps= 1.0E-9;
        auto inv_matrix  = LinearAlgebra::getInversMatrix(m_matrix,eps);

        DataComputing::writeMatrix(inv_matrix);
        auto res = MatrixOperations::matrixMultMatrix(m_matrix,inv_matrix,eps);
        DataComputing::writeMatrix(res);

    }
    

}

namespace Lab2 {
    void fixedPointIterationsDemo() {
        matrix equations;
        vector rightSide;
        vector x;
        DataComputing::readData(fileIn1,equations, rightSide);
        LinearAlgebra::fixedPointIteration(equations,rightSide,x,0.001);
        DataComputing::writeVector(fileOut1,x);

    }
}

int a() {
    //Lab2::fixedPointIterationsDemo();
    Lab1::inverseMatirxDemo();
    Lab1::gaussFullChoiseDemo();
    return 0;
}
