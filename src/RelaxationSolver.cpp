#include "RelaxationSolver.h"

RelaxtationSolver::RelaxtationSolver()
:SLESolver(), 
precision(0.01),
w(0)
{
}

RelaxtationSolver::RelaxtationSolver(const size_t normType, const MyType _precision, const MyType _w, const MyType epsilon)
:SLESolver(normType, epsilon),
precision(_precision),
w(_w)
{
}

void RelaxtationSolver::iteration(Matrix &x1, Matrix& x0,const Matrix& problem, Matrix& rs)
{
    for(int i = 0; i < x1.getRows(); ++i) {
    MyType sum1=0;
    for(int j = 0; j < i; ++j) {
        sum1+=((problem(j,i)/problem(i,i))*x1(j));
    }
    sum1*=-w;
    MyType sum2 =(1-w)*x0(i);
    MyType sum3=0;
    for(int j = i+1; j < problem.getRows(); ++j){
        sum3+=((problem(j,i)/problem(i,i))*x0(j));
    }
    sum3*=-w;
    MyType sum= sum1+sum2+sum3+w*(rs(i)/problem(i,i));
    x1(i) = sum;
}
}

Matrix RelaxtationSolver::getC(Matrix &problem)
{   
    Matrix C(problem.getCols(),(int)problem.getRows(),epsilon);
    Matrix E(problem.getRows(),(MyType)epsilon);
    Matrix L(problem.getCols(),(int)problem.getRows());
    Matrix U(problem.getCols(),(int)problem.getRows());
    Matrix D(problem.getCols(),(int)problem.getRows()); 
     for(int i = 0; i < C.getCols(); ++i) {
        for(int j = 0; j < C.getRows(); ++j) {
            if(i==j)D(i,j)=problem(i,j);
            if(i<j)L(i,j)=problem(i,j);
            if(i>j)U(i,j)=problem(i,j); 
        }
    }
    Matrix Di=D.getInverseMatrix();
    /* std::cout<<L<<"\n\n"<<D<<"\n\n"<<U<<"\n\n"; */
    C = (E+(Di*w)*L).getInverseMatrix()*((E*(1-w))-(Di*w)*U);
    return C;
}

Matrix RelaxtationSolver::solve(Matrix &problem, Matrix &xStart)
{
    problem.setPresision(epsilon);
    Matrix rs=problem.popCol(problem.getCols()-1);
    Matrix result(1,(int)problem.getRows(),epsilon);
    this->iteration(result,xStart,problem,rs);
    const auto Norm = this->getC(problem).norm(normType);
    size_t itCount=0;
    while((result-xStart).norm(normType)> ((1-Norm)* precision)/Norm ) {
        xStart = result;
        this->iteration(result,xStart,problem,rs);
        itCount+=1;
    } 
    std::cout<<itCount<<" iterations\n";
    return result;
}

Matrix RelaxtationSolver::solve(Matrix &problem)
{
    problem.setPresision(epsilon);
    Matrix startX(1,(int)problem.getRows(),epsilon);
    return (*this).solve(problem, startX);
}

Matrix RelaxtationSolver::solve(Matrix &&problem)
{
    return (*this).solve(problem);
}
