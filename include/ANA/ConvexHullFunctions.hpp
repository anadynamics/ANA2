#ifndef ANA_CONVEXHULLFUNCTIONS_H
#define ANA_CONVEXHULLFUNCTIONS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CCH_into_cavity(CCavity &hueco, CConvexHull const &CH);

// Returns the intersection point between the segment and the convex hull.
inline CPoint cget_intersection_point(Segment const &s, CConvexHull const &CH) {

    for (auto const &t : CH._triangles) {
        Object const inter_obj = CGAL::intersection(s, t);
        if (CPoint const *inter_point = CGAL::object_cast<CPoint>(&inter_obj)) {

            return *inter_point;
        }
    }
    throw("Fatal error, cget_intersection_point() could not find an "
          "intersection.");
}

// Returns the intersection point between the vector and the convex hull.
inline Point get_intersection_point(
    Point const &p_in, Point const &p_out, ConvexHull const &CH) {

    for (size_t i = 0; i < CH._normals.size(); ++i) {
        // Check if this triangle is on the same halfspace as this point.
        double const pdt = dot_product(normalize(p_out - p_in), CH._normals[i]);
        // if (pdt < 0) {
        //     std::abs(a - b) < a * .001
        // }
    }
    //  for (auto const &t :) {
    //     Object const inter_obj = CGAL::intersection(s, t);
    //     if (CPoint const *inter_point =
    //     CGAL::object_cast<CPoint>(&inter_obj)) {

    //         return *inter_point;
    //     }
    // }
    throw("Fatal error, could not find an intersection.");
}

// Get intersection points between the cell and the included area.
auto cget_vertices_3_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<CPoint, CPoint, CPoint, CPoint, double>;

// Get intersection points between the cell and the included area.
auto cget_vertices_2_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<CPoint, CPoint, CPoint,
    CPoint, CPoint, CPoint, double, double>;

// Get intersection points between the cell and the included area.
auto cget_vertices_1_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<CPoint, CPoint, CPoint,
    CPoint, CPoint, CPoint, double, double, double>;

// Check if p_test is on the same side v is pointing, with respect to p.
inline bool is_vtx_inside(
    Point const &test_point, Point const &p, Vector const v) {

    Vector const test_vtor = normalize(test_point - p);
    double const pdt = dot_product(test_vtor, v);
    return pdt < 0;
}

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH);

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Tetrahedron const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Tetrahedron const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Tetrahedron const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double,
        double>;

} // namespace ANA

#endif // _H