#include "LinearSolver.h"

#include <iostream>
#include <fstream>

LinearSolver::LinearSolver(): normType(0){
    
}

LinearSolver::LinearSolver(const size_t _normType):
normType(_normType)
{
}


/* Matrix LinearSolver::SimpleIterattions(Matrix problem, Matrix rs, Matrix xStart)
{
    Matrix E(problem.getRows());
     MyType tau = 0.01;
    Matrix C = (E-problem*tau);
    Matrix y = rs* tau;
    if(C.norm(normType)<1) {
        Matrix x= C*xStart+y;
        while((x-xStart).norm(normType)>
        ((1-C.norm(normType))* 0.01)/C.norm(normType)) {
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




Matrix LinearSolver::Jacobi(Matrix& problem, Matrix& rs, Matrix& xStart)
{

    Matrix C(problem.getCols(),problem.getRows());
/*     Matrix L(problem.getCols(),problem.getRows());
    Matrix U(problem.getCols(),problem.getRows());
    Matrix D(problem.getCols(),problem.getRows()); 
    Matrix y(1,problem.getRows());


    for(int i = 0; i < C.getCols(); ++i) {
        for(int j = 0; j < C.getRows(); ++j) {
/*              if(i==j)D(i,j)=problem(i,j);
            if(i<j)L(i,j)=problem(i,j);
            if(i>j)U(i,j)=problem(i,j);  
           if(i!=j) {
                C(j,i)=(problem(j,i)/problem(i,i))*(-1);
                y(i)=rs(i)/problem(i,i);
            } 
        }
    }
  //  std::cout<<L<<"\n\n"<<D<<"\n\n"<<U<<"\n\n";
    std::cout<<C<<"\n\n"<<y<<"\n\n";
        if(C.norm(normType)<1) {
        Matrix x= C*xStart+y;
        while((x-xStart).norm(normType)>
        ((1-C.norm(normType))* 0.01)
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
 */