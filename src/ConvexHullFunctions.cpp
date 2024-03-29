#include <ANA/ConvexHullFunctions.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH) {

    for (size_t k = 0; k < hueco._all_cells.size(); ++k) {
        Tetrahedron const &cell = hueco._all_cells[k];
        TetraInfo const &info = hueco._all_info[k];

        // Find which vertices are in and which are out of the CH. Also get the
        // total count of the out vertices.
        auto [vertices_out_cnt, vertices_in, vertices_out] =
            is_cell_inside(cell, CH);

        if (vertices_out_cnt == 0) {
            // cell is entirely contained in the included area.
            hueco.add_inner_cell(cell, info);
        } else if (vertices_out_cnt == 4) {
            // cell is outside the included area.
            continue;
        } else {
            // Cell has 1-3 vertices inside the included area.
            // Get the points of the intersections between the convex hull
            // and the tetrahedron's segments that go from the inner
            // point(s) to the outer point(s). The polyhedron delimited by these
            // points will be added to the cavity.
            hueco.add_outer_cell(cell, info);

            switch (vertices_out_cnt) {
            case 3: {
                auto const [p_in_0, ip1, ip2, ip3, radius_0] =
                    get_vertices_3_out(
                        cell, info, CH, vertices_in, vertices_out);

                hueco.add_border_tetra(p_in_0, ip1, ip2, ip3, radius_0);
                break;
            }
            case 2: {
                auto const [p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0,
                    radius_1] = get_vertices_2_out(cell, info, CH, vertices_in,
                    vertices_out);

                hueco.add_border_penta_2_2(
                    p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0, radius_1);
                break;
            }
            case 1: {
                auto const [p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, radius_0,
                    radius_1, radius_2] = get_vertices_1_out(cell, info, CH,
                    vertices_in, vertices_out);

                hueco.add_border_penta_3_1(p_in_0, p_in_1, p_in_2, ip1, ip2,
                    ip3, radius_0, radius_1, radius_2);
                break;
            }
            }
        }
    }

    return;
}

// Cicles through every CH triangle and its normals to find which vertices
// are IN/OUT.
auto is_cell_inside(Tetrahedron const &cell, ConvexHull const &CH)
    -> std::tuple<int, std::vector<int>, std::vector<int>> {

    std::vector<int> vertices_out, vertices_in;
    // Only inside points will give a >0 dot product against all normals.
    for (std::size_t i = 0; i < 4; ++i) {
        bool vtx_is_inside = true;
        Point const &test_point = cell[i];

        for (std::size_t j = 0; j < CH._triangles.size(); ++j) {
            vtx_is_inside =
                is_vtx_inside(test_point, CH._triangles[j][0], CH._normal_0[j]);
            if (!vtx_is_inside) {
                vertices_out.push_back(i);
                break;
            }
        }
        if (vtx_is_inside) {
            vertices_in.push_back(i);
        }
    }
    return {vertices_out.size(), vertices_in, vertices_out};
}

// Check if test_point is on the same side v is pointing, with respect to p.
// CH triangles normals are pointing outwards.
// If the vertex is close to a CH triangle, it will be moved along the
// triangle's normal
auto is_vtx_inside(
    Point const &test_point, Point const &p_0, Vector const &normal_0) -> bool {

    Vector const diff_vtor = normalize(test_point - p_0);
    double const test_dot_pdt0 = dot_product(diff_vtor, normal_0);

    bool const atroden = test_dot_pdt0 < 0.;

    return atroden;
}

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Tetrahedron const &cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, double> {

    Point const &p_in_0 = cell[vertices_in[0]];
    Point const &p_out_1 = cell[vertices_out[0]];
    Point const &p_out_2 = cell[vertices_out[1]];
    Point const &p_out_3 = cell[vertices_out[2]];

    Point const ip1 = get_intersection_point(p_in_0, p_out_1, CH);
    Point const ip2 = get_intersection_point(p_in_0, p_out_2, CH);
    Point const ip3 = get_intersection_point(p_in_0, p_out_3, CH);

    return {p_in_0, ip1, ip2, ip3, info._radius[vertices_in[0]]};
}

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Tetrahedron const cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double> {

    Point const &p_in_0 = cell[vertices_in[0]];
    Point const &p_in_1 = cell[vertices_in[1]];
    Point const &p_out_2 = cell[vertices_out[0]];
    Point const &p_out_3 = cell[vertices_out[1]];

    Point const ip1 = get_intersection_point(p_in_0, p_out_2, CH);
    Point const ip2 = get_intersection_point(p_in_0, p_out_3, CH);
    Point const ip3 = get_intersection_point(p_in_1, p_out_2, CH);
    Point const ip4 = get_intersection_point(p_in_1, p_out_3, CH);

    return {p_in_0, p_in_1, ip1, ip2, ip3, ip4, info._radius[vertices_in[0]],
        info._radius[vertices_in[1]]};
}

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Tetrahedron const cell, TetraInfo const &info,
    ConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<Point, Point, Point,
    Point, Point, Point, double, double, double> {

    Point const &p_in_0 = cell[vertices_in[0]];
    Point const &p_in_1 = cell[vertices_in[1]];
    Point const &p_in_2 = cell[vertices_in[2]];
    Point const &p_out_3 = cell[vertices_out[0]];

    Point ip1 = get_intersection_point(p_in_0, p_out_3, CH);
    Point ip2 = get_intersection_point(p_in_1, p_out_3, CH);
    Point ip3 = get_intersection_point(p_in_2, p_out_3, CH);

    return {p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, info._radius[vertices_in[0]],
        info._radius[vertices_in[1]], info._radius[vertices_in[2]]};
}

// Returns the intersection point between the segment that joins the 2 input
// points and the convex hull.
Point get_intersection_point(
    Point const &p_in, Point const &p_out, ConvexHull const &CH) {

    try {
        auto const intersection_point =
            try_get_intersection_point(p_in, p_out, CH);
        return intersection_point;

    } catch (const std::logic_error &le) {
        // Move the point that's outside the CH along each triangle normal until
        // it works.
        for (size_t i = 0; i < CH._triangles.size(); i++) {
            try {
                auto const intersection_point = try_get_intersection_point(
                    p_in, p_out + CH._normal_0[i], CH);
                return intersection_point;
            } catch (...) { }
        }
        // It didn't work.
        throw std::runtime_error(
            "Fatal error, get_intersection_point() could not find an "
            "intersection even after perturbing the out point.");
    }
}

// Returns the intersection point between the segment that joins the 2 input
// points and the convex hull.
Point try_get_intersection_point(
    Point const &p_in, Point const &p_out, ConvexHull const &CH) {

    Vector const s = p_out - p_in;
    // Check they are not the same.
    if ((s[0] < zero_top and s[0] > zero_bot) &&
        (s[1] < zero_top and s[1] > zero_bot) &&
        (s[2] < zero_top and s[2] > zero_bot)) {
        std::cerr << "p_in: " << p_in << " -- p_out: " << p_out << '\n';
        throw std::runtime_error(
            "Fatal error, try_get_intersection_point() could not "
            "find an intersection. IN/OUT points are the same.");
    }

    for (size_t t = 0; t < CH._triangles.size(); ++t) {

        Vector const d = p_out - CH._triangles[t][0];
        // First, check if the ray crosses the triangle's plane.
        double const dota = dot_product(normalize(s), CH._normal_0[t]);
        if (dota < zero_top) {
            continue;
        }

        double const det = determinant(s, CH._v01[t], CH._v02[t]);
        double const det2 = determinant(s, d, CH._v02[t]);
        double const det3 = determinant(s, CH._v01[t], d);

        double const l2 = det2 / det;
        double const l3 = det3 / det;
        double const l1 = 1 - l2 - l3;

        if (l1 > neg_delta and l2 > neg_delta and l3 > neg_delta) {
            double const sum = l1 + l2 + l3;
            if (sum > uno_bot and sum < uno_top) {
                return (CH._triangles[t][0] * l1 + CH._triangles[t][1] * l2 +
                    CH._triangles[t][2] * l3);
            }
        }
    }

    // SHOULD OUTPUT THIS WHEN IN DEBUG RELEASE
    // std::cerr << "p_in: " << p_in << " -- p_out: " << p_out << '\n';
    // std::cerr << "s: " << s << '\n';
    throw std::logic_error("");
}

} // namespace ANA
