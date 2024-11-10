#include "SLESolver.h"
SLESolver::SLESolver():LinearSolver()
{
}

SLESolver::SLESolver(const size_t normType, const MyType _epsilon):
LinearSolver(normType,_epsilon)
{
}

Matrix SLESolver::readLSE(const std::string &filename)
{
    std::ifstream file(filename);
    int N;
    file >> N;
    Matrix problem(N+1,N);
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N+1; ++j){
            file>>problem(j,i);
            if(std::abs(problem(j,i))<epsilon)problem(j,i)=0;
        }
    }
    return problem;
}

Matrix SLESolver::fsolve(const std::string filename)
{

    return (*this).solve(((*this).readLSE(filename)));
}
