#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Topology/CaribouTopology.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/Mass.h>
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
namespace sofa { using Index = unsigned int; }
#endif

namespace SofaCaribou::mass {

// Traits to get the Sofa vector type from the dimension
template <std::size_t Dim> struct SofaVecType {};
template <> struct SofaVecType<1> { using Type = sofa::defaulttype::Vec1Types; };
template <> struct SofaVecType<2> { using Type = sofa::defaulttype::Vec2Types; };
template <> struct SofaVecType<3> { using Type = sofa::defaulttype::Vec3Types; };

/**
 * Implementation of a consistent Mass matrix.
 *
 * The assembly of this mass matrix takes the form of
 *
 * \f{eqnarray*}{
 *     \mat{M}_{IK}  = \int_{\Omega_e} \rho_0 N_I N_K d\Omega \mat{I}
 * \f}
 *
 * where \f$K\f$ and \f$K\f$ are a pair of indices of the element \f$e\f$ nodes.
 *
 * @tparam Element The element type of this mass. It must inherits from caribou::geometry::Element.
 */
template <typename Element>
class CaribouMass : public sofa::core::behavior::Mass<typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(CaribouMass, Element), SOFA_TEMPLATE(sofa::core::behavior::Mass, typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type));

    // Type aliases
    using DataTypes = typename SofaVecType<caribou::geometry::traits<Element>::Dimension>::Type;
    using Coord    = typename DataTypes::Coord;
    using VecCoord = typename DataTypes::VecCoord;
    using DataVecCoord = typename sofa::core::objectmodel::Data<VecCoord>;
    using Deriv    = typename DataTypes::Deriv;
    using VecDeriv = typename DataTypes::VecDeriv;
    using DataVecDeriv = typename sofa::core::objectmodel::Data<VecDeriv>;
    using Real     = typename DataTypes::Real;

    template<int nRows, int nColumns>
    using Matrix = Eigen::Matrix<Real, nRows, nColumns, Eigen::RowMajor>;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<CaribouMass<Element>, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

    // Constant properties
    static constexpr INTEGER_TYPE Dimension = caribou::geometry::traits<Element>::Dimension;
    static constexpr INTEGER_TYPE NumberOfNodesPerElement = caribou::geometry::traits<Element>::NumberOfNodesAtCompileTime;
    static constexpr INTEGER_TYPE NumberOfGaussNodesPerElement = caribou::geometry::traits<Element>::NumberOfGaussNodesAtCompileTime;

    // Data structures
    struct GaussNode {
        Real weight;                          ///< Weight of this integration point
        Real jacobian_determinant;            ///< Determinant of the transformation matrix Jacobian, ie dx = J du
        Matrix<NumberOfNodesPerElement, 1> N; ///< shape value of each nodes evaluated at this integration point
    };

    // The container of Gauss points (for each elements) is an array if the number of integration
    // points per element is known at compile time, or a dynamic vector otherwise.
    using GaussContainer = typename std::conditional<
            NumberOfGaussNodesPerElement != caribou::Dynamic,
            std::array<GaussNode, static_cast<std::size_t>(NumberOfGaussNodesPerElement)>,
            std::vector<GaussNode>
    >::type;

    // Public methods
    CARIBOU_API
    CaribouMass();

    CARIBOU_API
    void init() override;

    template <typename Derived>
    CARIBOU_API
    static auto canCreate(Derived * o, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool;

    template <typename DerivedElement>
    static auto templateName(const CaribouMass<DerivedElement>* = nullptr) -> std::string {
        return SofaCaribou::topology::CaribouTopology<DerivedElement>::templateName(nullptr);
    }

    /** Get the number of elements contained in this field **/
    [[nodiscard]] inline
    auto number_of_elements() const noexcept -> std::size_t {
        return (p_topology ? p_topology->number_of_elements() : 0);
    }

    [[nodiscard]] inline
    auto topology() const noexcept -> typename SofaCaribou::topology::CaribouTopology<Element>::SPtr {
        return p_topology;
    }

    CARIBOU_API
    void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& /*f*/, const DataVecCoord& /*x*/, const DataVecDeriv& /*v*/) override {}

    /**
     * @return True if the mass matrix is lumped, false otherwise. The lumping is done by placing integration point
     * directly onto the nodes of each elements.
     */
    CARIBOU_API
    bool isDiagonal() const override {return d_lumped.getValue();}

    /**
     * This will force an assembly of the consistent mass matrix. It is called automatically during the
     * scene graph initialization phase. It is up to the user the call back this method when the topology
     * changes.
     */
    CARIBOU_API
    void update_mass_matrix();

    /** Get the set of Gauss integration nodes of an element */
    CARIBOU_API
    inline auto gauss_nodes_of(std::size_t element_id) const -> const auto & {
        return p_elements_quadrature_nodes[element_id];
    }

private:
    // These private methods are implemented but can be overridden

    /** Get the set of Gauss integration nodes of the given element */
    virtual auto get_gauss_nodes(const std::size_t & element_id, const Element & element) const -> GaussContainer;

    /** Compute and store the shape functions and their derivatives for every integration points */
    virtual void initialize_elements();

    // Data members
    /// This link is specifically set to point towards a very general BaseObject since it can be either a
    /// BaseMeshTopology (SOFA topology container), or a CaribouTopology.
    Link<sofa::core::objectmodel::BaseObject> d_topology_container;

    /// Whether or not the mass matrix should be lumped (diagonal) using
    /// integration points located at element node positions.
    sofa::core::objectmodel::Data<bool> d_lumped;

    /// Mass density of the material.
    sofa::core::objectmodel::Data<Real> d_density;

    // Private variables
    /// Pointer to a CaribouTopology. This pointer will be null if a CaribouTopology
    /// is found within the scene graph and linked using the d_topology_container data
    /// parameter. Otherwise, if a compatible SOFA's topology (see SofaCaribou::topology::CaribouTopology::mesh_is_compatible())
    /// is found and linked, an internal CaribouTopology component will be created
    /// and its pointer will be stored here.
    typename SofaCaribou::topology::CaribouTopology<Element>::SPtr p_topology;

    /// Consistent mass matrix (used when d_lumped == false)
    Eigen::SparseMatrix<Real> p_M;

    /// Diagonal mass matrix (used when d_lumped == true)
    Eigen::DiagonalMatrix<Real, Eigen::Dynamic> p_Mdiag;

    /// Integration points of each elements
    std::vector<GaussContainer> p_elements_quadrature_nodes;

};

} // namespace SofaCaribou::mass
