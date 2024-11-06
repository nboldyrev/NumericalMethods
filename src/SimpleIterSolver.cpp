#include "SimpleIterSolver.h"

SimpleIterSolver::SimpleIterSolver():LinearSolver(),precision(0.01)
{
}

Matrix SimpleIterSolver::solve(Matrix &problem, Matrix &xStart)
{
    Matrix rs=problem.popCol(problem.getCols()-1);
    Matrix E(problem.getRows());
    MyType tau=0.01;
    Matrix C(E-problem*tau);
    Matrix y=rs*tau;

      if(C.norm(normType)<1) {
        Matrix x= C*xStart+y;
        while((x-xStart).norm(normType)>
        ((1-C.norm(normType))* precision)
        /C.norm(normType)) {
            xStart = x;
            x=C*xStart+y;
        }
        return x;
      } 
    else{
        //TODO tau
    }
    return Matrix();
}