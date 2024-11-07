#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "Matrix.h"
#include "algorithm"
#include "fstream"
class LinearSolver {
    protected:
        size_t normType;
        size_t epsilon;
    public:
        LinearSolver();
        LinearSolver(const size_t _normType, const MyType epsilon = 2.20E-16);
};
#endif //LINEARSOLVER_H