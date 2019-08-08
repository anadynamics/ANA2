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

// Check if p_test is on the same side v is pointing, with respect to p.
inline bool is_vtx_inside(
    Point const &test_point, Point const &p, Vector const v) {

    Vector test_vtor = normalize(test_point - p);
    double const test_dot_pdt = dot_product(test_vtor, v);
    return test_dot_pdt < zero_bot;
}

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