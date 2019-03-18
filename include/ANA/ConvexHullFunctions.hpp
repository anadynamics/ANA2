#ifndef ANA_CONVEXHULLFUNCTIONS_H
#define ANA_CONVEXHULLFUNCTIONS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CCH_into_cavity(Cavity &hueco, ConvexHull const &CH);

// Returns the intersection point between the segment and the convex hull.
inline CPoint get_intersection_point(Segment const &s, ConvexHull const &CH) {

    for (auto const &t : CH._triangles) {
        Object const inter_obj = CGAL::intersection(s, t);
        if (CPoint const *inter_point = CGAL::object_cast<CPoint>(&inter_obj)) {
            return *inter_point;
        }
    }
    throw std::runtime_error("Fatal error, get_intersection_point() could not "
                             "find an intersection.");
}

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<CPoint, CPoint, CPoint, CPoint, double>;

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<CPoint, CPoint, CPoint, CPoint, CPoint, CPoint, double,
        double>;

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Finite_cells_iterator const cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<CPoint, CPoint, CPoint, CPoint, CPoint, CPoint, double,
        double, double>;

// Check if p_test is on the same side v is pointing, with respect to p.
inline bool is_vtx_inside(
    CPoint const &test_point, CPoint const &p, CVector const v) {

    CVector test_vtor = test_point - p;
    test_vtor =
        test_vtor / std::sqrt(CGAL::to_double(test_vtor.squared_length()));
    double const test_dot_pdt = CGAL::to_double(test_vtor * v);
    return test_dot_pdt > 0;
}

} // namespace ANA

#endif // _H