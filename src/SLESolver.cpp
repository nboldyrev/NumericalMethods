#include "SLESolver.h"
SLESolver::SLESolver():LinearSolver()
{
}

SLESolver::SLESolver(const size_t normType):LinearSolver(normType)
{
}

Matrix SLESolver::readLSE(const std::string &filename, MyType epsi)
{
    std::ifstream file(filename);
    int N;
    file >> N;
    Matrix problem(N+1,N);
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N+1; ++j){
            file>>problem(j,i);
        }
    }
    return problem;
}
