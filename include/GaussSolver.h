#ifndef GAUSSSOLVER_H
#define GAUSSSOLVER_H
#include "SLESolver.h"
class GaussSolver:public SLESolver {

    public:
        GaussSolver();
        GaussSolver(const size_t _normType, const MyType _epsilon=type<MyType>()(((MyType)0.1)));
        Matrix gaussForwardElim(Matrix& problem,std::vector<std::pair<size_t,size_t>>&swaps);
        Matrix gaussBackwardElim(Matrix& augProblem);
        Matrix gaussForwardElimMod(Matrix& problem,std::vector<std::pair<size_t,size_t>>&swaps);
        Matrix solve(Matrix& problem) override;
        Matrix solve(Matrix&& problem) override;
        Matrix solveMod(Matrix&A);
};
#endif //GAUSSSOLVER_H
