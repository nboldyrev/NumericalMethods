#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "Matrix.h"
#include "algorithm"
#include "fstream"
class LinearSolver {
    protected:
        size_t normType;
        MyType epsilon;
    public:
        LinearSolver();
        LinearSolver(const size_t _normType, const MyType epsilon = type<MyType>().getDefaultEps((MyType)0.1));
};
#endif //LINEARSOLVER_H