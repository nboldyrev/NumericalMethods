#include "QRSolver.h"
#include "GaussSolver.h"
std::pair<Matrix, Matrix> QRSolver::calcQR(Matrix &A)
{
	MyType s = 0, c = 0, l = 0;
    if(A.getCols()-A.getRows()!=0)A.popCol(A.getCols()-1);
    const int hight = A.getRows();
    const int width = A.getCols();
    Matrix T(hight);
	for (auto i = 0; i < hight - 1; ++i) {
		for (auto j = i + 1; j < hight; ++j) {
			l = std::sqrt(A(i,i) * A(i,i) + A(i,j) * A(i,j));
			if (l == 0) {
				continue;
			}
			c = A(i,i) / l;
			s = A(i,j) / l;
            auto iPred=A.getRow(i);
            auto jPred = A.getRow(j);
            auto iPredT=T.getRow(i);
            auto jPredT = T.getRow(j);
            A.addToRow(i,iPred,c-1);
            A.addToRow(i,jPred,s);
            A.addToRow(j,iPred,-s);
            A.addToRow(j,jPred,c-1);
            T.addToRow(i,iPredT,c-1);
            T.addToRow(i,jPredT,s);
            T.addToRow(j,iPredT,-s);
            T.addToRow(j,jPredT,c-1);
		}
	}
    Q=A;
    R=T.transpose();
    return std::pair<Matrix, Matrix>(Q,R);
}

std::pair<Matrix, Matrix> QRSolver::getQR()
{
    return std::pair<Matrix, Matrix>(Q,R);
}

std::pair<Matrix, Matrix> QRSolver::getQR(Matrix &A)
{
    (*this).calcQR(A);
    return std::pair<Matrix, Matrix>(Q,R);
}

Matrix QRSolver::solve(Matrix &rhs)
{
    if(rhs.getCols()==1){
        auto newRs=(R.transpose())*rhs;
        R.transpose();
        Q.append(newRs);
        GaussSolver a;
        auto x=a.gaussBackwardElim(Q);
        return x;
    }
    else {
        auto rs = rhs.popCol(rhs.getCols()-1);
        (*this).calcQR(rhs);
        return (*this).solve(rs);
    }

}

QRSolver::QRSolver() : SLESolver()
{
}

QRSolver::QRSolver(const size_t _normType):
 SLESolver(_normType)
{

}

QRSolver::QRSolver(Matrix& A):SLESolver()
{

     (*this).calcQR(A);
}
