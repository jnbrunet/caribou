#pragma once
#include <SofaCaribou/Forcefield/CaribouForcefield.h>
#include <SofaCaribou/Topology/CaribouTopology[Hexahedron].h>
#include <Caribou/Geometry/Hexahedron.h>

namespace SofaCaribou::forcefield {

// Hexahedron linear specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>::triangulate_face(const caribou::geometry::Hexahedron < caribou::Linear> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron < caribou::Linear>>;

// Hexahedron quadratic specialization
template <> auto CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::templateName(const CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>> *) -> std::string;
template <> void CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>::triangulate_face(const caribou::geometry::Hexahedron < caribou::Quadratic> & e, const std::size_t & face_id, std::vector<sofa::type::Vector3> & triangles_nodes);
extern template class CaribouForcefield<caribou::geometry::Hexahedron < caribou::Quadratic>>;

}