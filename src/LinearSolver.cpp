#include "LinearSolver.h"

LinearSolver::LinearSolver(const MyType _epsilon, const MyType _presicion, const size_t _iterrationLimit, const size_t _normType):
epsilon(_epsilon),
presicion(_presicion),
iterrationLimit(_iterrationLimit),
normType(_normType)
{}

