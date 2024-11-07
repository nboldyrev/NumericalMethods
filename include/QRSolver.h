#ifndef QRSOLVER_H
#define QRSOLVER_H
#include "SLESolver.h"
class QRSolver :public SLESolver {
    private:
        Matrix Q;
        Matrix R;
    public:
        QRSolver();
        QRSolver(const size_t _normType);
        QRSolver( Matrix& A);
        std::pair<Matrix,Matrix> calcQR(Matrix& A);
        std::pair<Matrix,Matrix> getQR();
        std::pair<Matrix,Matrix> getQR(Matrix& A);
        Matrix solve(Matrix& rs) override;
};
#endif //QRSOLVER_H