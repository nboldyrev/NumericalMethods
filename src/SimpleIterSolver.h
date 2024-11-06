#ifndef SIMPLEITERSOLVER_H
#define SIMPLEITERSOLVER_H
#include "LinearSolver.h"
class SimpleIterSolver: public LinearSolver{
    private:
    MyType precision;
    public:
        SimpleIterSolver();
        SimpleIterSolver(const size_t normType, const MyType presicion);
        Matrix solve(Matrix& problem, Matrix& xStart);
};
#endif //SIMPLEITERSOLVER_H