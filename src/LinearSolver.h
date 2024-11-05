#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H
#include "Matrix.h"
class LinearSolver {
    protected:
        MyType epsilon;
        MyType presicion;
        size_t iterrationLimit;
        size_t normType;


    public:
        LinearSolver(const MyType _epsilon = 2.20E-16,const MyType _presicion=0, const size_t _iterrationLimit = 14881488, const size_t _normType = 0 );

};
#endif LINEARSOLVER_H