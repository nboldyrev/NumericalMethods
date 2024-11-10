#ifndef RELAXATIONSOLVER_H
#include "SLESolver.h"
class RelaxtationSolver: public SLESolver {
    private:
        MyType precision;
        MyType w;
    public:
        RelaxtationSolver();
        RelaxtationSolver(const size_t normType, const MyType presicion, const MyType w, const MyType epsilon =type<MyType>()(((MyType)0.1)));
        void iteration(Matrix &x1,  Matrix &x0,const Matrix& problem, Matrix& rs);
        Matrix getC(Matrix& problem);
        Matrix solve(Matrix& problem, Matrix& xStart);
        Matrix solve(Matrix& problem) override;
        Matrix solve(Matrix&& problem) override;
};
#define RELAXATIONSOLVER_H
#endif //RELAXATIONSOLVER_H