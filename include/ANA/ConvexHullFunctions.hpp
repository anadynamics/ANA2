#ifndef ANA_CONVEXHULLFUNCTIONS_H
#define ANA_CONVEXHULLFUNCTIONS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/PDB.hpp>
#include <ANA/PrimitivesUtils.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH);

// Cicles through every CH triangle and its normals to find which vertices are
// IN/OUT.
auto is_cell_inside(Tetrahedron const &cell, ConvexHull const &CH)
    -> std::tuple<int, std::vector<int>, std::vector<int>>;

// Get dot product between triangle normal and diff vector between test_point
// and p.
inline double which_side(
    Point const &test_point, Point const &p, Vector const &normal) {

    Vector const diff_vtor = normalize(test_point - p);
    return dot_product(diff_vtor, normal);
}

// Check if test_point is on the same side v is pointing, with respect to p.
// CH triangles normals are pointing outwards.
// If the dot product gives a value between zero_bot and zero_bot_3, check again
// with the other 2 triangle points. If they don't add up to 3*zero_bot_3, then
// it's out. This is to prevent numerical errors from vertices too close to CH
// surface.
auto is_vtx_inside(Point const &test_point, Triangle const &triangle,
    Vector const &normal_0, Vector const &normal_1, Vector const &normal_2)
    -> bool;

// Returns the intersection point between the segment that joins the 2 input
// points and the convex hull.
Point get_intersection_point(
    Point const &r, Point const &q, ConvexHull const &CH);

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Tetrahedron const &cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Tetrahedron const cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Tetrahedron const cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<Point, Point, Point,
    Point, Point, Point, double, double, double>;

} // namespace ANA

#endif // _H