#ifndef SLESOLVER_H
#define SLESOLVER_H
#include "LinearSolver.h"
class SLESolver: public LinearSolver
{
public:
    SLESolver();
    SLESolver(const size_t normType);
    Matrix readLSE(const std::string& filename,MyType epsi=2.20E-16);    
    const MyType conditionNumber(Matrix& problem) const;  
    virtual Matrix solve(Matrix& solve)=0;
};
#endif //SLESOLVER_H
