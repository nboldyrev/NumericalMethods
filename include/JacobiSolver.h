#ifndef JACOBISOLVER_H
#define JACOBISOLVER_H
#include "SLESolver.h"
class JacobiSolver: public SLESolver{
    private:
    MyType precision;
    public:
        JacobiSolver();
        JacobiSolver(const size_t normType, const MyType presicion, const MyType epsilon = 2.20E-16);
        Matrix solve(Matrix& problem, Matrix& xStart);
        Matrix solve(Matrix& problem) override;
        Matrix solve(Matrix&& problem) override;
};
#endif 