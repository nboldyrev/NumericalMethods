#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "Matrix.h"
#include "algorithm"
class LinearSolver {
    protected:
        size_t normType;
        

    public:
        LinearSolver();
        LinearSolver(const size_t _normType);
        Matrix readProblem(const std::string& filename,MyType epsi=2.20E-16);
/*         Matrix SimpleIterattions(Matrix problem, Matrix rs, Matrix xStart);
        Matrix Jacobi(Matrix& problem, Matrix& rs, Matrix& xStart);
        Matrix Zeidel(Matrix problem);
        Matrix Relaxation(Matrix provlem); */
};
#endif //LINEARSOLVER_H