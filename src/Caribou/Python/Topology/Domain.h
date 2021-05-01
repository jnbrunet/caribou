#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>

#include <Caribou/Python/Caribou.h>
#include <Caribou/Python/Topology/BarycentricContainer.h>

namespace caribou::topology::bindings {
template <typename Element, typename NodeIndex>
void declare_domain(pybind11::module & m, const std::string & name) {
    using D = Domain<Element, NodeIndex>;
    pybind11::class_<D, std::unique_ptr<D, py::nodelete>> c(m, name.c_str());

    c.def("canonical_dimension", &D::canonical_dimension);
    c.def("number_of_nodes_per_elements", &D::number_of_nodes_per_elements);
    c.def("number_of_elements", &D::number_of_elements);
    c.def("element", [](const D & domain, const UNSIGNED_INTEGER_TYPE & element_id) {
        return domain.element(element_id);
    }, py::arg("element_id"));

    c.def("element_indices", &D::element_indices, py::arg("index"));
    c.def("mesh", &D::mesh);

    // Mesh's add_domain binding for Domain<Element, NodeIndex> type
    m.def("add_domain", [](BaseMesh & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime> & node_indices) {
        auto d = new Domain<Element, NodeIndex>(nullptr, node_indices);
        mesh.add_domain(d, domain_name);
        return d;
    }, py::arg("mesh"), py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"), py::return_value_policy::reference);
    m.def("add_domain", [](BaseMesh & mesh, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, geometry::traits<Element>::NumberOfNodesAtCompileTime> & node_indices) {
        auto d = new Domain<Element, NodeIndex>(nullptr, node_indices);
        mesh.add_domain(d);
        return d;
    }, py::arg("mesh"), py::arg("element_type"), py::arg("node_indices"), py::return_value_policy::reference);

    m.def("add_domain", [](BaseMesh & mesh, const std::string & domain_name, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic> & node_indices) {
        if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and node_indices.cols() != geometry::traits<Element>::NumberOfNodesAtCompileTime) {
            std::ostringstream ss;
            ss << "You tried to create a domain containing elements of " << geometry::traits<Element>::NumberOfNodesAtCompileTime
               << " nodes, but elements having " << node_indices.cols() << " nodes were found instead. The node indices "
               << "matrix should have NxM indices where N is the number of elements, and M is the number of nodes per element.";
            throw py::key_error(ss.str());
        }

        auto d = new Domain<Element, NodeIndex>(nullptr, node_indices);
        mesh.add_domain(d, domain_name);
        return d;
    }, py::arg("mesh"), py::arg("domain_name"), py::arg("element_type"), py::arg("node_indices"), py::return_value_policy::reference);

    m.def("add_domain", [](BaseMesh & mesh, const Element &, const Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic> & node_indices) {
        if (geometry::traits<Element>::NumberOfNodesAtCompileTime != caribou::Dynamic and node_indices.cols() != geometry::traits<Element>::NumberOfNodesAtCompileTime) {
            std::ostringstream ss;
            ss << "You tried to create a domain containing elements of " << geometry::traits<Element>::NumberOfNodesAtCompileTime
               << " nodes, but elements having " << node_indices.cols() << " nodes were found instead. The node indices "
               << "matrix should have NxM indices where N is the number of elements, and M is the number of nodes per element.";
            throw py::key_error(ss.str());
        }

        auto d = new Domain<Element, NodeIndex>(nullptr, node_indices);
        mesh.add_domain(d);
        return d;
    }, py::arg("mesh"), py::arg("element_type"), py::arg("node_indices"), py::return_value_policy::reference);

    declare_barycentric_container(c);
}

} // namespace caribou::topology::bindings