#ifndef SLESOLVER_H
#define SLESOLVER_H
#include "LinearSolver.h"
class SLESolver: public LinearSolver
{

public:
    SLESolver();
    SLESolver(const size_t normType,const MyType _epsilon=2.20E-16);
    Matrix readLSE(const std::string& filename); 
    Matrix solve(const std::string filename);  
    const MyType conditionNumber(Matrix& problem) const;  
    virtual Matrix solve(Matrix& problem)=0;
    virtual Matrix solve(Matrix&& problem)=0;
    Matrix fsolve(const std::string filename);
};
#endif //SLESOLVER_H
