#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "Matrix.h"
#include "algorithm"
#include "fstream"
class LinearSolver {
    protected:
        size_t normType;
    public:
        LinearSolver();
        LinearSolver(const size_t _normType);
};
#endif //LINEARSOLVER_H