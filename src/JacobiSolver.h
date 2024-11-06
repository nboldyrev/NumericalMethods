#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H
#include "LinearSolver.h"
class JacobiSolver: LinearSolver{
    private:
    MyType precision;
    public:
        JacobiSolver();
        JacobiSolver(const size_t normType, const MyType presicion);
        Matrix solve(Matrix& problem, Matrix& xStart);
};
#endif 