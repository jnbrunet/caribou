#include <SofaCaribou/Solver/LUSolver.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::solver {

static int SparseLUSolverClass = sofa::core::RegisterObject("Caribou Sparse LU linear solver")
    .add< LUSolver<Eigen::SparseLU<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::AMDOrdering<int>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< LUSolver<Eigen::PardisoLU<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>> >()
#endif
;

} // namespace SofaCaribou::solver