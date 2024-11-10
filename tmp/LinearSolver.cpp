#include "LinearSolver.h"

#include <iostream>
#include <fstream>

LinearSolver::LinearSolver(): normType(0), epsilon(type<MyType>()(((MyType)0.1))){
    
}

LinearSolver::LinearSolver(const size_t _normType, const MyType _epsilon):
normType(_normType),
epsilon(_epsilon)
{
}
