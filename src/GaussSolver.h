#ifndef GAUSSSOLVER_H
#define GAUSSSOLVER_H
#include "LinearSolver.h"
class GaussSolver:public LinearSolver {

    public:
        GaussSolver();
        GaussSolver(const size_t _normType);
        Matrix gaussForwardElim(Matrix& problem,std::vector<std::pair<size_t,size_t>>&swaps);
        Matrix gaussBackwardElim(Matrix& augProblem);
        Matrix gaussModifiedElim(Matrix& problem, Matrix& rs,std::vector<std::pair<size_t,size_t>>&swaps);
        Matrix solve(Matrix& problem);
        Matrix solveModified(Matrix&A);
};
#endif //GAUSSSOLVER_H
