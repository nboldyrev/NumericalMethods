#ifndef SIMPLEITERSOLVER_H
#define SIMPLEITERSOLVER_H
#include "SLESolver.h"
class SimpleIterSolver: public SLESolver{
    private:
    MyType precision;
    public:
        SimpleIterSolver();
        SimpleIterSolver(const size_t normType, const MyType presicion);
        Matrix solve(Matrix& problem, Matrix& xStart);
        Matrix solve(Matrix& problem) override;
};
#endif //SIMPLEITERSOLVER_H