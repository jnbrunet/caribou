#pragma once

#include <SofaCaribou/Mass/CaribouMass.h>
#include <Caribou/Geometry/Tetrahedron.h>

namespace SofaCaribou::mass {

// Tetrahedron linear specialization
extern template class CaribouMass<caribou::geometry::Tetrahedron<caribou::Linear>>;

// Tetrahedron quadratic specialization
extern template class CaribouMass<caribou::geometry::Tetrahedron<caribou::Quadratic>>;

} // namespace SofaCaribou::mass