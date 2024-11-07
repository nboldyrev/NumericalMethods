#include "SLESolver.h"
#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H
class JacobiSolver: SLESolver{
    private:
    MyType precision;
    public:
        JacobiSolver();
        JacobiSolver(const size_t normType, const MyType presicion);
        Matrix solve(Matrix& problem, Matrix& xStart);
        Matrix solve(Matrix& problem) ;
};
#endif 