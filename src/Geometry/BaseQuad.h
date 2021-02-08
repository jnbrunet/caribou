#pragma once

#include <Caribou/config.h>
#include <Caribou/constants.h>
#include <Caribou/Geometry/Element.h>
#include <Eigen/Core>

namespace caribou::geometry {

template<typename Derived>
struct BaseQuad : public Element<Derived> {
    // Types
    using Base = Element<Derived>;

    using LocalCoordinates = typename Base::LocalCoordinates;
    using WorldCoordinates = typename Base::WorldCoordinates;

    using GaussNode = typename Base::GaussNode;

    template <UNSIGNED_INTEGER_TYPE Dim>
    using Vector = typename Base::template Vector<Dim>;

    template <UNSIGNED_INTEGER_TYPE Rows, UNSIGNED_INTEGER_TYPE Cols>
    using Matrix = typename Base::template Matrix<Rows, Cols>;

    // Constants
    static constexpr auto CanonicalDimension = Base::CanonicalDimension;
    static constexpr auto Dimension = Base::Dimension;
    static constexpr auto NumberOfNodesAtCompileTime = Base::NumberOfNodesAtCompileTime;
    static constexpr auto NumberOfGaussNodesAtCompileTime = Base::NumberOfGaussNodesAtCompileTime;

    static_assert(Dimension == 2 or Dimension == 3, "Quads can only be of dimension 2 or 3.");

    /** Default empty constructor */
    BaseQuad() = default;

    /** Constructor from an Eigen matrix containing the positions of the quad's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseQuad(Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from an Eigen matrix containing the positions of the quad's nodes */
    template<typename EigenType, REQUIRES(EigenType::RowsAtCompileTime == NumberOfNodesAtCompileTime)>
    explicit BaseQuad(const Eigen::EigenBase<EigenType> & nodes) :p_nodes(nodes.derived().template cast<typename Base::Scalar>()) {}

    /** Constructor from a serie of nodes. */
    template <
        typename ...Nodes,
        REQUIRES(NumberOfNodesAtCompileTime == sizeof...(Nodes)+1)
    >
    explicit BaseQuad(const WorldCoordinates & first_node, Nodes&&...remaining_nodes)
    {
        construct_from_nodes<0>(first_node, std::forward<Nodes>(remaining_nodes)...);
    }

    // Public methods common to all quad types

    /**
     * Get the list of node indices of the edges.
     * \sa Element::boundary_elements_node_indices
     */
    inline auto edges() const {
        return self().get_boundary_elements_nodes();
    }

    /**
     * Extract the orthogonal frame of the element by computing the cross product of the unit vectors
     * from the center position to its projection on opposite faces.
     *
     * \sa BaseQuad::frame(const LocalCoordinates & local_point)
     *
     * \warning If the quad isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted quads.
     */
    inline auto frame() const -> Matrix<Dimension, Dimension>
    {
        return frame(LocalCoordinates::Zero());
    }

    /**
     * Extract the frame positioned at the given position (in local coordinates) on the quad by computing the cross
     * product of the unit vectors from the given position its projection on opposite edges.
     *
     * This function will return a matrix of the form:
     *
     *    2D              3D
     *
     * | ux vx |     | ux vx wx |
     * | uy vy |     | uy vy wy |
     *               | uz vz wz |
     *
     * Where (ux, uy[, uz]), (vx, vy[, vz]) and (wx, wy[, wz]) are orthogonal unitary vectors representing
     * the u, v [and w] frames in the current quad. If the quad is rectangular and not rotated, this matrix is the
     * Identity matrix. If it is rectangular but rotated, rotating the quad by the transposed of this frame should
     * align the u,v [,w] axis to the x,y[,z] world frame (identity matrix).
     *
     * \warning If the quad isn't rectangular, the frame extracted by this function will be a rough
     *          approximation that could be far from the real solution, especially for strongly deformed or
     *          inverted quads.
     */
    inline
    auto
    frame(const LocalCoordinates & local_point) const -> Matrix<Dimension, Dimension>
    {
        // Position of the point inside the quad where the frame should be computed
        const auto p = this->world_coordinates( local_point );

        // Project of the point on the edge facing the u axis
        const auto projected_on_u = this->world_coordinates({1,local_point[1]});

        // Project of the point on the edge facing the v axis
        const auto projected_on_v = this->world_coordinates({local_point[0], 1});

        // Vector from the point to its projection on the edge facing the u axis
        const auto point_to_u = projected_on_u - p;

        // Vector from the point to its projection on the edge facing the v axis
        const auto point_to_v = projected_on_v - p;

        // The u-axis of the computed frame
        const auto u = point_to_u.normalized();

        // The v-axis of the computed frame
        auto v = point_to_v.normalized();

        Matrix<Dimension, Dimension> m {};
        if constexpr (Dimension == 3) {
            // The w-axis of the computed frame
            const WorldCoordinates w = u.cross(v).normalized();
            m << u, v, w;
        } else {
            m << u, v;
        }

        return m;
    }

private:
    // Implementations
    friend struct Element<Derived>;
    [[nodiscard]]
    inline auto get_number_of_nodes() const {return NumberOfNodesAtCompileTime;}
    inline auto get_number_of_gauss_nodes() const {return NumberOfGaussNodesAtCompileTime;}
    inline auto get_node(const UNSIGNED_INTEGER_TYPE & index) const {return WorldCoordinates(p_nodes.row(index));};
    inline auto get_nodes() const -> const auto & {return p_nodes;};
    inline auto get_center() const {return Base::world_coordinates(LocalCoordinates({0, 0}));};
    inline auto get_number_of_boundary_elements() const -> UNSIGNED_INTEGER_TYPE {return 4;};
    inline auto get_contains_local(const LocalCoordinates & xi, const FLOATING_POINT_TYPE & eps) const -> bool {
        const auto & u = xi[0];
        const auto & v = xi[1];
        return IN_CLOSED_INTERVAL(-1-eps, u, 1+eps) and
               IN_CLOSED_INTERVAL(-1-eps, v, 1+eps);
    }

    auto self() -> Derived& { return *static_cast<Derived*>(this); }
    auto self() const -> const Derived& { return *static_cast<const Derived*>(this); }

    template <size_t index, typename ...Nodes, REQUIRES(sizeof...(Nodes) >= 1)>
    inline
    void construct_from_nodes(const WorldCoordinates & first_node, Nodes&&...remaining_nodes) {
        p_nodes.row(index) = first_node;
        construct_from_nodes<index+1>(std::forward<Nodes>(remaining_nodes)...);
    }

    template <size_t index>
    inline
    void construct_from_nodes(const WorldCoordinates & last_node) {
        p_nodes.row(index) = last_node;
    }
protected:
    Matrix<NumberOfNodesAtCompileTime, Dimension> p_nodes;
};

}