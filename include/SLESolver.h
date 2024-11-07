#ifndef SLESOLVER_H
#define SLESOLVER_H
#include "LinearSolver.h"
class SLESolver: public LinearSolver
{

public:
    SLESolver();
    SLESolver(const size_t normType,const MyType _epsilon=2.20E-16);
    Matrix readLSE(const std::string& filename);   
    const MyType conditionNumber(Matrix& problem) const;  
    virtual Matrix solve(Matrix& solve)=0;
};
#endif //SLESOLVER_H
