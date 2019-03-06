#include <ANA/ConvexHullFunctions.hpp>

namespace ANA {

// Discard voids outside the convex hull.
void carve_CCH_into_cavity(CCavity &hueco, CConvexHull const &CH) {

    for (auto const &cell : hueco._all_cells) {
        std::vector<int> vertices_out, vertices_in;
        // Only inside points will give a >0 dot product against all normals.
        for (std::size_t i = 0; i < 4; ++i) {
            bool vtx_is_inside = true;
            CPoint const test_point(cell->vertex(i)->point());

            for (std::size_t j = 0; j < CH._triangles.size(); j++) {
                CVector test_vtor = test_point - CH._triangles[j][1];
                test_vtor = test_vtor /
                    std::sqrt(CGAL::to_double(test_vtor.squared_length()));
                double const test_dot_pdt =
                    CGAL::to_double(test_vtor * CH._normals[j]);
                if (test_dot_pdt < 0) {
                    vtx_is_inside = false;
                    vertices_out.push_back(i);
                    break;
                }
            }
            if (vtx_is_inside) {
                vertices_in.push_back(i);
            }
        }
        auto const vertices_out_cnt = vertices_out.size();
        if (vertices_out_cnt == 0) {
            // cell is entirely contained in the included area.
            hueco._volume = hueco._volume + volume(cell);
            hueco._included_cells.push_back(cell);
        } else if (vertices_out_cnt == 4) {
            // cell is outside the included area.
            continue;
        } else {
            // Cell has 1-3 vertices inside the included area.
            // Get the points of the intersections between the convex hull
            // and the tetrahedron's segments that go from the inner
            // point(s) to the outer point(s). The polyhedron delimited by these
            // points will be added to the cavity.
            switch (vertices_out_cnt) {
            case 3: {
                auto const [p_in_0, ip1, ip2, ip3, radius_0] =
                    cget_vertices_3_out(cell, CH, vertices_in, vertices_out);

                hueco.cadd_border_tetra(p_in_0, ip1, ip2, ip3, radius_0);
                break;
            }
            case 2: {
                auto const [p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0,
                    radius_1] =
                    cget_vertices_2_out(cell, CH, vertices_in, vertices_out);

                hueco.cadd_border_penta(
                    p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0, radius_1);
                break;
            }
            case 1: {
                auto const [p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, radius_0,
                    radius_1, radius_2] =
                    cget_vertices_1_out(cell, CH, vertices_in, vertices_out);

                hueco.cadd_border_penta(p_in_0, p_in_1, p_in_2, ip1, ip2, ip3,
                    radius_0, radius_1, radius_2);
                break;
            }
            }
        }
    }

    return;
}

// Get intersection points between the cell and the included area.
auto cget_vertices_3_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out)
    -> std::tuple<CPoint, CPoint, CPoint, CPoint, double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;
    CPoint p_in_0(vtx_in_0->point());

    CPoint const p_out_1(cell->vertex(vertices_out[0])->point());
    CPoint const p_out_2(cell->vertex(vertices_out[1])->point());
    CPoint const p_out_3(cell->vertex(vertices_out[2])->point());

    CPoint const ip1 = cget_intersection_point(Segment(p_in_0, p_out_1), CH);
    CPoint const ip2 = cget_intersection_point(Segment(p_in_0, p_out_2), CH);
    CPoint const ip3 = cget_intersection_point(Segment(p_in_0, p_out_3), CH);

    return {p_in_0, ip1, ip2, ip3, radius_0};
}

// Get intersection points between the cell and the included area.
auto cget_vertices_2_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<CPoint, CPoint, CPoint,
    CPoint, CPoint, CPoint, double, double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;
    CPoint p_in_0(vtx_in_0->point());

    auto const vtx_in_1{cell->vertex(vertices_in[1])};
    double const radius_1 = vtx_in_1->info()._radius;
    CPoint p_in_1(vtx_in_1->point());

    CPoint const p_out_2(cell->vertex(vertices_out[0])->point());
    CPoint const p_out_3(cell->vertex(vertices_out[1])->point());

    CPoint const ip1 = cget_intersection_point(Segment(p_in_0, p_out_2), CH);
    CPoint const ip2 = cget_intersection_point(Segment(p_in_0, p_out_3), CH);
    CPoint const ip3 = cget_intersection_point(Segment(p_in_1, p_out_2), CH);
    CPoint const ip4 = cget_intersection_point(Segment(p_in_1, p_out_3), CH);

    return {p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0, radius_1};
}

// Get intersection points between the cell and the included area.
auto cget_vertices_1_out(Finite_cells_iterator const cell,
    CConvexHull const &CH, std::vector<int> const &vertices_in,
    std::vector<int> const &vertices_out) -> std::tuple<CPoint, CPoint, CPoint,
    CPoint, CPoint, CPoint, double, double, double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;
    CPoint p_in_0(vtx_in_0->point());

    auto const vtx_in_1{cell->vertex(vertices_in[1])};
    double const radius_1 = vtx_in_1->info()._radius;
    CPoint p_in_1(vtx_in_1->point());

    auto const vtx_in_2{cell->vertex(vertices_in[2])};
    double const radius_2 = vtx_in_2->info()._radius;
    CPoint p_in_2(vtx_in_2->point());

    CPoint const p_out_3(cell->vertex(vertices_out[0])->point());

    CPoint const ip1 = cget_intersection_point(Segment(p_in_0, p_out_3), CH);
    CPoint const ip2 = cget_intersection_point(Segment(p_in_1, p_out_3), CH);
    CPoint const ip3 = cget_intersection_point(Segment(p_in_2, p_out_3), CH);

    return {
        p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, radius_0, radius_1, radius_2};
}

// Discard voids outside the convex hull.
void carve_CH_into_cavity(Cavity &hueco, ConvexHull const &CH) {

    for (Tetrahedron const &cell : hueco._all_cells) {
        std::vector<int> vertices_out, vertices_in;
        // Only inside points will give a >0 dot product against all normals.
        for (std::size_t i = 0; i < 4; ++i) {
            bool vtx_is_inside = true;
            Point const test_point(cell[i]);
            for (std::size_t j = 0; j < CH._triangles.size(); j++) {

                vtx_is_inside = is_vtx_inside(
                    test_point, CH._triangles[j][1], CH._normals[j]);

                if (!vtx_is_inside) {
                    vertices_out.push_back(i);
                    break;
                }
            }
            if (vtx_is_inside) {
                vertices_in.push_back(i);
            }
        }
        auto const vertices_out_cnt = vertices_out.size();
        if (vertices_out_cnt == 0) {
            // cell is entirely contained in the included area.
            hueco._volume = hueco._volume + volume(cell);
            hueco._included_cells.push_back(cell);
        } else if (vertices_out_cnt == 4) {
            // cell is outside the included area.
            continue;
        } else {
            // Cell has 1-3 vertices inside the included area.
            // Get the points of the intersections between the convex hull
            // and the tetrahedron's segments that go from the inner
            // point(s) to the outer point(s). The polyhedron delimited by
            // these points will be added to the cavity.
            switch (vertices_out_cnt) {
            case 3: {
                auto const [p_in_0, ip1, ip2, ip3, radius_0] =
                    get_vertices_3_out(cell, CH, vertices_in, vertices_out);

                hueco.add_border_tetra(p_in_0, ip1, ip2, ip3, radius_0);
                break;
            }
            case 2: {
                auto const [p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0,
                    radius_1] =
                    get_vertices_2_out(cell, CH, vertices_in, vertices_out);

                hueco.add_border_penta(
                    p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0, radius_1);
                break;
            }
            case 1: {
                auto const [p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, radius_0,
                    radius_1, radius_2] =
                    get_vertices_1_out(cell, CH, vertices_in, vertices_out);

                hueco.add_border_penta(p_in_0, p_in_1, p_in_2, ip1, ip2, ip3,
                    radius_0, radius_1, radius_2);
                break;
            }
            }
        }
    }

    return;
}

// Get intersection points between the cell and the included area.
auto get_vertices_3_out(Tetrahedron const &cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;

    Point const ip1 =
        get_intersection_point(cell[vertices_in[0]], cell[vertices_out[0]] CH);
    Point const ip2 =
        get_intersection_point(cell[vertices_in[1]], cell[vertices_out[1]] CH);
    Point const ip3 =
        get_intersection_point(cell[vertices_in[2]], cell[vertices_out[2]] CH);

    return {p_in_0, ip1, ip2, ip3, radius_0};
}

// Get intersection points between the cell and the included area.
auto get_vertices_2_out(Tetrahedron const &cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;
    Point p_in_0(vtx_in_0->point());

    auto const vtx_in_1{cell->vertex(vertices_in[1])};
    double const radius_1 = vtx_in_1->info()._radius;
    Point p_in_1(vtx_in_1->point());

    Point const p_out_2(cell->vertex(vertices_out[0])->point());
    Point const p_out_3(cell->vertex(vertices_out[1])->point());

    Point const ip1 = get_intersection_point(p_in_0, p_out_2, CH);
    Point const ip2 = get_intersection_point(p_in_0, p_out_3, CH);
    Point const ip3 = get_intersection_point(p_in_1, p_out_2, CH);
    Point const ip4 = get_intersection_point(p_in_1, p_out_3, CH);

    return {p_in_0, p_in_1, ip1, ip2, ip3, ip4, radius_0, radius_1};
}

// Get intersection points between the cell and the included area.
auto get_vertices_1_out(Tetrahedron const &cell, ConvexHull const &CH,
    std::vector<int> const &vertices_in, std::vector<int> const &vertices_out)
    -> std::tuple<Point, Point, Point, Point, Point, Point, double, double,
        double> {

    auto const vtx_in_0{cell->vertex(vertices_in[0])};
    double const radius_0 = vtx_in_0->info()._radius;
    Point p_in_0(vtx_in_0->point());

    auto const vtx_in_1{cell->vertex(vertices_in[1])};
    double const radius_1 = vtx_in_1->info()._radius;
    Point p_in_1(vtx_in_1->point());

    auto const vtx_in_2{cell->vertex(vertices_in[2])};
    double const radius_2 = vtx_in_2->info()._radius;
    Point p_in_2(vtx_in_2->point());

    Point const p_out_3(cell->vertex(vertices_out[0])->point());

    Point const ip1 = get_intersection_point(p_in_0, p_out_3, CH);
    Point const ip2 = get_intersection_point(p_in_1, p_out_3, CH);
    Point const ip3 = get_intersection_point(p_in_2, p_out_3, CH);

    return {
        p_in_0, p_in_1, p_in_2, ip1, ip2, ip3, radius_0, radius_1, radius_2};
}

} // namespace ANA
