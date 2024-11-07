#include "GaussSolver.h"
 GaussSolver::GaussSolver():SLESolver()
{
}

GaussSolver::GaussSolver(const size_t _normType, const MyType _epsilon)
:SLESolver(_normType,_epsilon)
{
}

Matrix GaussSolver::gaussForwardElim(Matrix &problem, std::vector<std::pair<size_t, size_t>> &swaps)
{
    problem.setPresision(epsilon);
    problem.toUpperTriangleForm(swaps);
    return problem;
}

Matrix GaussSolver::gaussBackwardElim(Matrix &augProblem)
{
    augProblem.setPresision(epsilon);
    const auto hight = augProblem.getRows();
    const auto width = augProblem.getCols();
    Matrix result(1,(int)hight);
    for(int i = width-1; i>=0; --i) {
        result(i) = augProblem(width-1,i);
        for(int j = width-2; j>i;--j) {
                result(i)-=augProblem(j,i)*result(j);
        }
       result(i)/=augProblem(i,i);
    }
    return result;
}


Matrix GaussSolver::solve(Matrix &A) 
{
    A.setPresision(epsilon);
    std::vector<std::pair<size_t,size_t>> swaps;
    (*this).gaussForwardElim(A,swaps);
    auto result=(*this).gaussBackwardElim(A);
    const auto N = swaps.size(); 
    for(int i = 0; i < N; ++i) {
       result.swapRows(i,swaps[i].first);
    } 
    std::cout<<result<<"\n\n";
   return Matrix();
}

Matrix GaussSolver::solve(Matrix &&problem)
{
    return (*this).solve(problem);
}
