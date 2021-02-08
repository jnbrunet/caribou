#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/OptionsGroup.h>

#include <SofaCaribou/Topology/FictitiousGrid.h>

#include <Eigen/Sparse>

#include <Caribou/Geometry/Hexahedron.h>
#include <Caribou/Geometry/RectangularHexahedron.h>

namespace SofaCaribou::forcefield {

using namespace sofa::core;
using namespace sofa::core::objectmodel;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;
using sofa::defaulttype::Vec3Types;

class FictitiousGridElasticForce : public ForceField<Vec3Types>
{
public:
    SOFA_CLASS(FictitiousGridElasticForce, SOFA_TEMPLATE(ForceField, Vec3Types));

    // Type definitions
    using DataTypes = Vec3Types;
    using Inherit  = ForceField<DataTypes>;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using Real     = typename Coord::value_type;

    using FictitiousGrid = SofaCaribou::topology::FictitiousGrid<DataTypes>;

    using Hexahedron = caribou::geometry::Hexahedron<caribou::Linear>;
    using RectangularHexahedron = caribou::geometry::RectangularHexahedron<caribou::Linear>;
    static constexpr INTEGER_TYPE NumberOfNodes = caribou::geometry::traits<Hexahedron>::NumberOfNodesAtCompileTime;


    template<int nRows, int nColumns, int Options=Eigen::RowMajor>
    using Matrix = Eigen::Matrix<Real, nRows, nColumns, Options>;

    template<int nRows, int nColumns>
    using Map = Eigen::Map<const Matrix<nRows, nColumns, Eigen::RowMajor>>;

    template<int nRows, int Options=0>
    using Vector = Eigen::Matrix<Real, nRows, 1, Options>;

    template<int nRows>
    using MapVector = Eigen::Map<const Vector<nRows, Eigen::ColMajor>>;

    using Mat33   = Matrix<3, 3, Eigen::RowMajor>;
    using Vec3   = Vector<3>;
    using Mat2424 = Matrix<24, 24, Eigen::RowMajor>;
    using Vec24   = Vector<24>;

    template <typename ObjectType>
    using Link = SingleLink<FictitiousGridElasticForce, ObjectType, BaseLink::FLAG_STRONGLINK>;

    // Data structures

    struct GaussNode {
        Real weight = 0;
        Matrix<NumberOfNodes, 3> dN_dx = Matrix<NumberOfNodes, 3, Eigen::RowMajor>::Zero();
        Mat33 F = Mat33::Identity();
    };

    /// Integration method used to integrate the stiffness matrix.
    enum class IntegrationMethod : unsigned int {
        /// Regular 8 points gauss integration
            Regular = 0,

        /// Hexas are recursively subdivided into cuboid subcells and the later
        /// are used to compute the inside volume of the regular hexa's gauss points.
        /// ** Requires a sparse grid topology **
            SubdividedVolume = 1,

        /// Hexas are recursively subdivided into cuboid subcells and the later
        /// are used to add new gauss points. Gauss points outside of the boundary are ignored.
        /// ** Requires a sparse grid topology **
            SubdividedGauss = 2
    };

    // Public methods

    FictitiousGridElasticForce();

    void init() override;
    void reinit() override;

    void addForce(
        const MechanicalParams* mparams,
        Data<VecDeriv>& d_f,
        const Data<VecCoord>& d_x,
        const Data<VecDeriv>& d_v) override;

    void addDForce(
        const MechanicalParams* /*mparams*/,
        Data<VecDeriv>& /*d_df*/,
        const Data<VecDeriv>& /*d_dx*/) override;

    void draw(const sofa::core::visual::VisualParams* vparams) override;

    SReal getPotentialEnergy(
        const MechanicalParams* /* mparams */,
        const Data<VecCoord>& /* d_x */) const override
    {return 0;}

    void addKToMatrix(sofa::defaulttype::BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override;

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

    template <typename T>
    inline
    Hexahedron hexahedron(std::size_t hexa_id, const T & x) const
    {
        auto * grid = d_grid_container.get();
        const auto &node_indices = grid->get_node_indices_of(hexa_id);

        Matrix<8, 3> m;
        for (std::size_t j = 0; j < 8; ++j) {
            const auto &node_id = node_indices[j];
            m.row(j) = MapVector<3>(&x[node_id][0]);
        }

        return Hexahedron(m);
    }

    inline
    IntegrationMethod integration_method() const
    {
        const auto m = static_cast<IntegrationMethod> (d_integration_method.getValue().getSelectedId());
        switch (m) {
            case IntegrationMethod::Regular:
                return IntegrationMethod::Regular;
            case IntegrationMethod::SubdividedVolume:
                return IntegrationMethod::SubdividedVolume;
            case IntegrationMethod::SubdividedGauss:
                return IntegrationMethod::SubdividedGauss;
            default:
                return IntegrationMethod::Regular;
        }
    }

    inline
    std::string integration_method_as_string() const
    {
        return d_integration_method.getValue().getSelectedItem();
    }

    inline
    void set_integration_method(const IntegrationMethod & m) {
        auto index = static_cast<unsigned int > (m);
        sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> d = d_integration_method;
        d->setSelectedItem(index);
    }

    const std::vector<GaussNode> & gauss_nodes_of(std::size_t hexahedron_id) const {
        return p_quadrature_nodes[hexahedron_id];
    }

    const Matrix<24, 24> & stiffness_matrix_of(std::size_t hexahedron_id) const {
        return p_stiffness_matrices[hexahedron_id];
    }

    /** Get the complete tangent stiffness matrix */
    const Eigen::SparseMatrix<Real> & K();

    /** Get the eigen values of the tangent stiffness matrix */
    const Vector<Eigen::Dynamic> & eigenvalues();

    /** Get the condition number of the tangent stiffness matrix */
    Real cond();

private:
    /** (Re)Compute the tangent stiffness matrix */
    virtual void compute_K();

protected:
    Data< Real > d_youngModulus;
    Data< Real > d_poissonRatio;
    Data< bool > d_linear_strain;
    Data< bool > d_corotated;
    Data< sofa::helper::OptionsGroup > d_integration_method;
    Link<FictitiousGrid> d_grid_container;

private:
    bool recompute_compute_tangent_stiffness = false;
    std::vector<Matrix<24, 24>> p_stiffness_matrices;
    std::vector<std::vector<GaussNode>> p_quadrature_nodes;
    std::vector<Mat33> p_initial_rotation;
    std::vector<Mat33> p_current_rotation;
    Eigen::SparseMatrix<Real> p_K;
    Vector<Eigen::Dynamic> p_eigenvalues;
    bool K_is_up_to_date;
    bool eigenvalues_are_up_to_date;

};

} // namespace SofaCaribou::forcefield
