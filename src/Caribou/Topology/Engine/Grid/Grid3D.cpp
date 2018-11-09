#include <Caribou/Topology/Engine/Grid/Grid.h>

namespace caribou
{

namespace topology
{

namespace engine
{


template <>
Grid3D::Index
Grid3D::cell_index(const VecInt & grid_coordinates) const
{
    const VecInt::ValueType & i = grid_coordinates[0];
    const VecInt::ValueType & j = grid_coordinates[1];
    const VecInt::ValueType & k = grid_coordinates[2];

    const auto & n = number_of_subdivision();
    const auto & nx = n[0];
    const auto & ny = n[1];
    const auto & nz = n[2];

    if (i > nx-1 || j > ny-1 || k > nz-1) {
        throw std::out_of_range(
                "Trying to access a cell at an invalid grid coordinate (" + std::to_string(i) + ", " + std::to_string(j) + "," + std::to_string(k) + ")"
        );
    }

    // Index is generated by looking at the cells as a flatten array
    return k*ny*nx + j*nx + i;
};

template <>
Grid3D::VecInt
Grid3D::grid_coordinates(const Index & cell_index) const
{
    const auto & n = number_of_subdivision();
    const auto & nx = n[0];
    const auto & ny = n[1];

    const Index k = cell_index / (nx*ny);
    const Index j = (cell_index - (k*nx*ny)) / nx;
    const Index i = cell_index - ((k*nx*ny) + (j*nx));

    return {i, j, k};
};

template <>
std::array<typename Grid3D::Index, Grid3D::NumberOfNodes>
Grid3D::nodes(const VecInt & grid_coordinates) const
{
    const VecInt::ValueType & i = grid_coordinates[0];
    const VecInt::ValueType & j = grid_coordinates[1];
    const VecInt::ValueType & k = grid_coordinates[2];

    const auto & n = number_of_subdivision();
    const auto & nx = n[0]+1; // Number of nodes in the x direction
    const auto & ny = n[1]+1; // Number of nodes in the y direction
    const auto & nz = n[2]+1; // Number of nodes in the z direction

    if (i > nx || j > ny || k > nz) {
        throw std::out_of_range(
                "Trying to access a cell at an invalid grid coordinate (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")"
        );
    }

    return {
            (k+0)*ny*nx + (j+0)*nx + (i+0),
            (k+0)*ny*nx + (j+0)*nx + (i+1),
            (k+0)*ny*nx + (j+1)*nx + (i+1),
            (k+0)*ny*nx + (j+1)*nx + (i+0),
            (k+1)*ny*nx + (j+0)*nx + (i+0),
            (k+1)*ny*nx + (j+0)*nx + (i+1),
            (k+1)*ny*nx + (j+1)*nx + (i+1),
            (k+1)*ny*nx + (j+1)*nx + (i+0)
    };
};

template <>
Grid3D::VecFloat
Grid3D::position(const Index & node_id) const
{
    const auto & n = number_of_subdivision();
    const auto & nx = n[0]+1; // Number of nodes in the x direction
    const auto & ny = n[1]+1; // Number of nodes in the y direction

    const Index k = node_id / (nx*ny); // Node indice in the z direction
    const Index j = (node_id - (k*nx*ny)) / nx; // Node indice in the y direction
    const Index i = node_id - ((k*nx*ny) + (j*nx)); // Node indice in the x direction

    const auto & h = cell_size();
    const auto & hx = h[0]; // Width of a cell
    const auto & hy = h[1]; // Height of a cell
    const auto & hz = h[2]; // Dept of a cell

    // Relative position of the node within the grid
    const VecFloat p = {
            i*hx,
            j*hy,
            k*hz
    };

    return anchor + p; // World position
};

template <>
Grid3D::Grid(VecFloat anchor, VecInt subdivisions, VecFloat dimensions)
        : anchor(anchor), nSubdivisions(subdivisions), dimensions(dimensions)
{
    const auto & nx = subdivisions[0];
    const auto & ny = subdivisions[1];
    const auto & nz = subdivisions[2];

    // Here we loop on the z, y and x order to proper align in memory cells so that when loping in the inverse
    // order (x, y and z), the cells are closer to each other in the memory
    cells.resize(nz*ny*nx);
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                Index index = cell_index({i, j, k});
                cells[index].reset(new CellType (this, index));
            }
        }
    }
};

} // namespace engine

} // namespace topology

} // namespace caribou