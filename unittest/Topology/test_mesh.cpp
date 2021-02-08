#ifndef CARIBOU_TOPOLOGY_TEST_MESH_H
#define CARIBOU_TOPOLOGY_TEST_MESH_H

#include <gtest/gtest.h>
#include "topology_test.h"
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Tetrahedron.h>

TEST(Mesh, Mesh) {
    using namespace caribou;
    using namespace caribou::geometry;
    using caribou::topology::Mesh;


    {
        // Construct by copy (from std::vector)
        std::vector<Mesh<3>::WorldCoordinates> initial_positions = {{0,0,0}, {1,1,1}};
        Mesh<3> mesh (initial_positions);
        auto positions = mesh.positions({1, 0});
        EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
        EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
    }

    {
        // Construct by copy (from Eigen::Matrix)
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3> initial_positions;
        initial_positions << 0, 0, 0,
                             1, 1, 1;
        Mesh<3> mesh (initial_positions);
        auto positions = mesh.positions({1, 0});
        EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
        EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
    }

    {
        // Construct by reference (from Eigen::Matrix)
        Eigen::Matrix<FLOATING_POINT_TYPE, 2, 3> initial_positions;
        initial_positions << 0, 0, 0,
            1, 1, 1;
        Mesh mesh (&initial_positions);
        auto positions = mesh.positions({1, 0});
        EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
        EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
        EXPECT_MATRIX_EQUAL(mesh.position(0), initial_positions.row(0));
        EXPECT_MATRIX_EQUAL(mesh.position(1), initial_positions.row(1));
        EXPECT_EQ(&(mesh.position(0).coeff(0)), &(initial_positions.row(0).coeff(0)));
    }

    {
        // Construct by reference (from Eigen::Map)
        std::vector<float> initial_positions = {
            0, 1, 2, 3, 4, 5, 6,
            7, 8, 9, // Ignored
            10, 11, 12, 13, 14, 15, 16,
            17, 18, 19, // Ignored
            20, 21, 22, 23, 24, 25, 26
        };
        Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor>, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 3>> map(initial_positions.data(), {10, 3});
        Mesh mesh (map);
        using WorldCoordinates = Eigen::Matrix<float, 3, 1>;
        auto positions = mesh.positions({1, 0});
        EXPECT_MATRIX_EQUAL(mesh.position(1), positions.row(0));
        EXPECT_MATRIX_EQUAL(mesh.position(0), positions.row(1));
        EXPECT_MATRIX_EQUAL(mesh.position(0), WorldCoordinates(0, 3, 6).transpose());
        EXPECT_MATRIX_EQUAL(mesh.position(1), WorldCoordinates(10, 13, 16).transpose());
        EXPECT_MATRIX_EQUAL(mesh.position(2), WorldCoordinates(20, 23, 26).transpose());
        EXPECT_EQ(&(mesh.position(0).coeff(0)), initial_positions.data());
    }
}

TEST(Mesh, Segment) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Segment<_3D>>;

    Mesh mesh({{0, 0, 0}, {1, 1, 1}});

    Domain::ElementsIndices domain_indices(1, 2); // One element of two nodes
    domain_indices << 0, 1;

    Domain * domain = mesh.add_domain<Segment<3>>("segments", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_EQ(indices, Domain::ElementIndices ({0, 1}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Segment<3>>("segments", domain_indices));

    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Segment<3> segment(positions);
        EXPECT_MATRIX_EQUAL(segment.center(), Segment<3>::WorldCoordinates(0.5, 0.5, 0.5));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

TEST(Mesh, Triangle) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Triangle<_3D>>;

    Mesh mesh ({{50,50,33}, {60, 50, 21}, {55, 55, -4}});

    Domain::ElementsIndices domain_indices(1, 3); // One element of three nodes
    domain_indices << 2, 1, 0;

    Domain * domain = mesh.add_domain<Triangle<_3D, Linear>>("triangles", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_MATRIX_EQUAL(indices, Domain::ElementIndices({2, 1, 0}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Triangle<3>>("triangles", domain_indices));


    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Triangle<3> triangle (positions);
        EXPECT_MATRIX_EQUAL(triangle.node(0), Triangle<_3D>::WorldCoordinates(55, 55, -4));
    }

    {
        auto triangle = domain->element(0);
        EXPECT_MATRIX_EQUAL(triangle.node(0), Triangle<_3D>::WorldCoordinates(55, 55, -4));
        EXPECT_MATRIX_EQUAL(triangle.node(1), Triangle<_3D>::WorldCoordinates(60, 50, 21));
        EXPECT_MATRIX_EQUAL(triangle.node(2), Triangle<_3D>::WorldCoordinates(50,50,33));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

TEST(Mesh, Tetrahedron) {
    using namespace caribou;
    using namespace caribou::topology;
    using namespace caribou::geometry;
    using Mesh = Mesh<_3D>;
    using Domain = Mesh::Domain<Tetrahedron<Linear>>;

    Mesh mesh ({{50,50,33}, {50,50,33}, {55, 55, -4}, {55, 55, -40}});

    Domain::ElementsIndices domain_indices(1, 4); // One element of four nodes
    domain_indices << 3, 2, 1, 0;

    Domain * domain = mesh.add_domain<Tetrahedron<Linear>>("tetras", domain_indices);

    EXPECT_EQ(domain->number_of_elements(), 1);

    {
        const auto &indices = domain->element_indices(0);
        EXPECT_MATRIX_EQUAL(indices, Domain::ElementIndices({3, 2, 1, 0}));
    }

    EXPECT_ANY_THROW(mesh.add_domain<Tetrahedron<Linear>>("tetras", domain_indices));


    {
        auto positions = mesh.positions(domain->element_indices(0));
        const Tetrahedron<Linear> tetra (positions);
        EXPECT_MATRIX_EQUAL(tetra.node(0), Tetrahedron<Linear>::WorldCoordinates(55, 55, -40));
    }

    {
        auto tetra = domain->element(0);
        EXPECT_MATRIX_EQUAL(tetra.node(1), Tetrahedron<Linear>::WorldCoordinates(55, 55, -4));
    }

    mesh.remove(domain);
    EXPECT_EQ(mesh.number_of_domains(), 0);
}

#endif //CARIBOU_TOPOLOGY_TEST_MESH_H