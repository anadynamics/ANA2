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

auto is_vtx_inside(
    Point const &test_point, Point const &p_0, Vector const &normal_0) -> bool;

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

// Returns the intersection point between the segment that joins the 2 input
// points and the convex hull.
auto get_intersection_point(
    Point const &p_in, Point const &p_out, ConvexHull const &CH) -> Point;

// Returns the intersection point between the segment that joins the 2 input
// points and the convex hull.
auto try_get_intersection_point(
    Point const &p_in, Point const &p_out, ConvexHull const &CH) -> Point;

} // namespace ANA

#endif // _H