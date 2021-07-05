#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <SofaCaribou/Topology/CaribouTopology.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/MechanicalParams.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::mass {

template<typename Element>
CaribouMass<Element>::CaribouMass()
: d_topology_container(initLink(
        "topology",
        "Topology container containing the elements on which this mass will be computed."))
, d_lumped(initData(
        &d_lumped,
        false,
        "lumped",
        "Whether or not the mass matrix should be lumped (diagonal) "
        "using integration points located at element node positions."))
, d_density(initData(
        &d_density,
        Real(1),
        "density",
        "Mass density of the material."))
{}

template<typename Element>
void CaribouMass<Element>::init() {
    using sofa::core::topology::BaseMeshTopology;
    using sofa::core::objectmodel::BaseContext;
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    Inherit1::init();

    auto *context = this->getContext();

    if (!this->mstate) {
        msg_warning() << "No mechanical object found in the current context node. The data parameter "
                      << "'" << this->mstate.getName() << "' can be use to set the path to a mechanical "
                      << "object having a template of '" << DataTypes::Name() << "'";
        return;
    }

    // If not topology is specified, try to find one automatically in the current context
    if (not d_topology_container.get()) {
        // No topology specified. Try to find one suitable.
        auto caribou_containers = context->template getObjects<CaribouTopology>(BaseContext::Local);
        auto sofa_containers = context->template getObjects<BaseMeshTopology>(BaseContext::Local);
        std::vector<BaseMeshTopology *> sofa_compatible_containers;
        for (auto container : sofa_containers) {
            if (CaribouTopology::mesh_is_compatible(container)) {
                sofa_compatible_containers.push_back(container);
            }
        }
        if (caribou_containers.empty() and sofa_compatible_containers.empty()) {
            msg_warning() << "Could not find a topology container in the current context. "
                          << "Please add a compatible one in the current context or set the "
                          << "container's path using the '" << d_topology_container.getName()
                          << "' data parameter.";
        } else {
            if (caribou_containers.size() + sofa_compatible_containers.size() > 1) {
                msg_warning() << "Multiple topologies were found in the context node. "
                              << "Please specify which one contains the elements on "
                              << "which this force field will be applied "
                              << "by explicitly setting the container's path in the  '"
                              << d_topology_container.getName() << "' data parameter.";
            } else {
                // Prefer caribou's containers first
                if (not caribou_containers.empty()) {
                    d_topology_container.set(caribou_containers[0]);
                } else {
                    d_topology_container.set(sofa_compatible_containers[0]);
                }
            }

            msg_info() << "Automatically found the topology '" << d_topology_container.get()->getPathName()
                       << "'.";
        }
    }

    // Create a caribou internal Domain over the topology
    if (d_topology_container.get()) {
        auto sofa_topology = dynamic_cast<BaseMeshTopology *>(d_topology_container.get());
        auto caribou_topology = dynamic_cast<CaribouTopology *>(d_topology_container.get());
        if (sofa_topology) {
            // Initialize a new caribou topology from the SOFA topology
            p_topology = sofa::core::objectmodel::New<CaribouTopology>();
            p_topology->findData("indices")->setParent(CaribouTopology::get_indices_data_from(sofa_topology));
#if (defined(SOFA_VERSION) && SOFA_VERSION < 201299)
            using CaribouTopologyMechanicalLink = typename SofaCaribou::topology::CaribouTopology<Element>::template Link<sofa::core::State<DataTypes>>;
            auto state_link = dynamic_cast<CaribouTopologyMechanicalLink*>(p_topology->findLink("state"));
            state_link->set(this->getMState());
#else
            p_topology->findLink("state")->set(this->getMState());
#endif
            p_topology->init();
        } else {
            // A Caribou topology already exists in the scene
            p_topology = caribou_topology;
        }

        if (number_of_elements() == 0) {
            msg_warning() << "No element found in the topology '" << d_topology_container.get()->getPathName() << "'";
        }
    }

    initialize_elements();
}

template<typename Element>
template<typename Derived>
auto CaribouMass<Element>::canCreate(Derived *o, sofa::core::objectmodel::BaseContext *context,
                                     sofa::core::objectmodel::BaseObjectDescription *arg) -> bool {
    using namespace sofa::core::objectmodel;
    using CaribouTopology = SofaCaribou::topology::CaribouTopology<Element>;

    std::string requested_element_type = arg->getAttribute( "template", "");
    std::string this_element_type = Derived::templateName(o);

    // to lower
    std::string requested_element_type_lower = requested_element_type, this_element_type_lower = this_element_type;
    std::transform(requested_element_type.begin(), requested_element_type.end(), requested_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(this_element_type.begin(), this_element_type.end(), this_element_type_lower.begin(), [](unsigned char c){ return std::tolower(c); });

    if (requested_element_type_lower == this_element_type_lower) {
        return Inherit1::canCreate(o, context, arg);
    }

    if (not requested_element_type.empty()) {
        arg->logError("Requested element type is not '"+this_element_type+"'.");
        return false;
    }

    std::string topology_path = arg->getAttribute("topology", "");
    if (not topology_path.empty()) {
        topology_path = topology_path.substr(1); // removes the "@"
        // Make sure the specified topology has elements of type Element
        auto topology = context->get<BaseObject>(topology_path);
        auto caribou_topology = dynamic_cast<CaribouTopology *>(topology);
        auto sofa_topology = dynamic_cast<sofa::core::topology::BaseMeshTopology *>(topology);
        if (not caribou_topology and (not sofa_topology or not CaribouTopology::mesh_is_compatible(sofa_topology))) {
            arg->logError("Cannot deduce the element type from the specified mesh topology '" + topology_path + "'.");
            return false;
        }
    } else {
        // Try to find a compatible topology in the current context
        BaseObject * topology = nullptr;
        auto objects = context->getObjects<BaseObject>(BaseContext::SearchDirection::Local);
        for (auto * object : objects) {
            auto caribou_topology = dynamic_cast<const CaribouTopology *>(object);
            auto sofa_topology = dynamic_cast<const sofa::core::topology::BaseMeshTopology *>(object);
            if (caribou_topology  or (sofa_topology and CaribouTopology::mesh_is_compatible(sofa_topology))) {
                topology = object;
                break;
            }
        }
        if (not topology) {
            arg->logError("Cannot find a topology in the current context from which the template '"+this_element_type+"' can be deduced.");
            return false;
        }

        if (Inherit1::canCreate(o, context, arg)) {
            arg->setAttribute("topology", "@" + topology->getPathName());
            return true;
        } else {
            return false;
        }
    }

    return Inherit1::canCreate(o, context, arg);
}

template<typename Element>
void CaribouMass<Element>::update_mass_matrix() {
    using namespace sofa::core::objectmodel;

    const auto density = d_density.getValue();
    if (density < std::numeric_limits<Real>::epsilon()) {
        return;
    }

    const auto state = this->mstate.get();
    if (not state) {
        return;
    }

    const auto lumped = d_lumped.getValue();

    static const auto Id = Matrix<3, 3>::Identity();
    const auto nb_elements = this->number_of_elements();

    const sofa::helper::ReadAccessor<Data<VecCoord>> X = this->mstate->readRestPositions();
    const auto nDofs = X.size() * 3;

    /// Triplets are used to store matrix entries before the call to 'compress'.
    /// Duplicates entries are summed up. This is only used when the matrix isn't
    /// lumped to a diagonal matrix
    std::vector<Eigen::Triplet<Real>> triplets;

    if (lumped) {
        // Only the lumped (diagonal) mass matrix is initialized
        p_Mdiag.setZero(nDofs);
        p_M.resize(0, 0);
    } else {
        // Only the sparse mass matrix is initialized
        p_Mdiag.resize(0);
        p_M.resize(nDofs, nDofs);
        p_M.setZero();
        triplets.reserve(nDofs*3);
    }

    sofa::helper::AdvancedTimer::stepBegin("CaribouMass::update_mass_matrix");
    for (int element_id = 0; element_id < static_cast<int>(nb_elements); ++element_id) {
        // Fetch the node indices of the element
        auto node_indices = this->topology()->domain()->element_indices(element_id);

        if (lumped) { // This condition should be optimized out of the loop by the compiler
            // Only the lumped (diagonal) mass matrix is computed here. To do so, we use the node location
            // for positioning the integration points. This way, the shape value of a node i evaluated at any
            // node j != i  will be 0.
        } else {
            for (const auto & gauss_node : gauss_nodes_of(element_id)) {
                // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
                const auto & detJ = gauss_node.jacobian_determinant;

                // Gauss quadrature node weight
                const auto & w = gauss_node.weight;

                // Shape functions at gauss node
                const auto & N = gauss_node.N;

                // Computation of the consistent mass sub-matrix M_ij
                for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
                    const auto & Ni = N[i];
                    for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {

                    }
                }
            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("CaribouMass::update_mass_matrix");
}

template <typename Element>
auto CaribouMass<Element>::get_gauss_nodes(const std::size_t & /*element_id*/, const Element & element) const -> GaussContainer {
    GaussContainer gauss_nodes {};
    if constexpr (NumberOfGaussNodesPerElement == caribou::Dynamic) {
        gauss_nodes.resize(element.number_of_gauss_nodes());
    }

    const auto nb_of_gauss_nodes = gauss_nodes.size();
    for (std::size_t gauss_node_id = 0; gauss_node_id < nb_of_gauss_nodes; ++gauss_node_id) {
        const auto & g = element.gauss_node(gauss_node_id);

        const auto J = element.jacobian(g.position);
        const auto detJ = std::abs(J.determinant());

        // Derivatives of the shape functions at the gauss node with respect to global coordinates x,y and z
        const Matrix<NumberOfNodesPerElement, 1> N = element.L(g.position);


        GaussNode & gauss_node = gauss_nodes[gauss_node_id];
        gauss_node.weight               = g.weight;
        gauss_node.jacobian_determinant = detJ;
        gauss_node.N                    = N;
    }

    return gauss_nodes;
}

template <typename Element>
void CaribouMass<Element>::initialize_elements()
{
    using namespace sofa::core::objectmodel;

    sofa::helper::AdvancedTimer::stepBegin("CaribouMass::initialize_elements");

    if (!this->mstate)
        return;

    // Resize the container of elements'quadrature nodes
    const auto nb_elements = this->number_of_elements();
    if (p_elements_quadrature_nodes.size() != nb_elements) {
        p_elements_quadrature_nodes.resize(nb_elements);
    }

    // Translate the Sofa's mechanical state vector to Eigen vector type
    sofa::helper::ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
    const Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), sofa_x0.size(), Dimension);

    // Loop on each element and compute the shape functions and their derivatives for every of their integration points
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

        // Get an Element instance from the Domain
        const auto initial_element = this->topology()->element(element_id);

        // Fill in the Gauss integration nodes for this element
        p_elements_quadrature_nodes[element_id] = get_gauss_nodes(element_id, initial_element);
    }

    // Compute the volume
    Real v = 0.;
    for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {
        for (const auto & gauss_node : gauss_nodes_of(element_id)) {
            // Jacobian of the gauss node's transformation mapping from the elementary space to the world space
            const auto detJ = gauss_node.jacobian_determinant;

            // Gauss quadrature node weight
            const auto w = gauss_node.weight;

            v += detJ*w;
        }
    }

    msg_info() << "Total mass of the geometry is " << v*d_density.getValue();

    sofa::helper::AdvancedTimer::stepEnd("CaribouMass::initialize_elements");
}

} // namespace SofaCaribou::mass