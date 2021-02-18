#pragma once

#include <SofaCaribou/Solver/LLTSolver.h>
#include <SofaCaribou/Solver/EigenSparseSolver.inl>

#include<Eigen/SparseCholesky>

#include <algorithm>
#include <cctype>
#include <string>

#ifdef CARIBOU_WITH_MKL

// Bug introduced in Eigen 3.3.8, fixed in bfdd4a9
#ifndef EIGEN_USING_STD
#define EIGEN_USING_STD(a) EIGEN_USING_STD_MATH(a)
#endif

#include <Eigen/PardisoSupport>
#endif

namespace SofaCaribou::solver {

namespace {
template <typename T>
struct solver_traits {};

template<typename MatrixType, int UpLo, typename Ordering>
struct solver_traits<Eigen::SimplicialLLT < MatrixType, UpLo, Ordering>> {
static auto BackendName() -> std::string { return "Eigen"; }
};

#ifdef CARIBOU_WITH_MKL
template<typename MatrixType, int UpLo>
struct solver_traits <Eigen::PardisoLLT< MatrixType, UpLo >> {
    static auto BackendName() -> std::string {return "Pardiso";}
};
#endif
}

template<class EigenSolver>
std::string LLTSolver<EigenSolver>::BackendName() {
    return solver_traits<EigenSolver>::BackendName();
}

template<typename EigenSolver>
LLTSolver<EigenSolver>::LLTSolver()
: d_backend(initData(&d_backend
, "backend"
, R"(
    Solver backend used.

    Available backends are:
    Eigen:   Eigen LLT solver (SimplicialLLT) [default].
    Pardiso: Pardiso LLT solver.
  )", true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
{
    d_backend.setValue(sofa::helper::OptionsGroup(std::vector < std::string > {
            "Eigen", "Pardiso"
    }));


    // Put the backend name in lower case
    std::string backend_str = solver_traits<EigenSolver>::BackendName();
    std::transform(backend_str.begin(), backend_str.end(), backend_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    sofa::helper::WriteAccessor <Data<sofa::helper::OptionsGroup >> backend = d_backend;
    if (backend_str == "pardiso") { // Case insensitive
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }

    // Explicitly state that the matrix is symmetric (would not be possible to do an LLT decomposition otherwise)
    this->set_symmetric(true);
}


} // namespace SofaCaribou::solver
