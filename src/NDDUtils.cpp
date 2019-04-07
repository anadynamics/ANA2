#include <ANA/NDDUtils.hpp>

namespace ANA::NDD {

// NDD Specific function for PDB input
void ndd_read_PDB_get_cells(std::string const &filename,
    NDD_IVector const &in_void_cells_indices, NDD_Vector &output_cells) {

    // Read molecule
    chemfiles::Trajectory in_traj(filename);
    chemfiles::Frame in_frame = in_traj.read();
    auto in_xyz = in_frame.positions();
    chemfiles::Topology in_top = in_frame.topology();
    VertexInfo vi1;
    NDD_Element temp_cell_coords;

    for (NDD_IElement const &idx : in_void_cells_indices) {
        // Iterate over each cell of the current pocket
        for (std::size_t i = 0; i < 4; ++i) {
            // Iterate over each vertex of the current cell
            CPoint const p0 =
                CPoint(in_xyz[idx[i]][0], in_xyz[idx[i]][1], in_xyz[idx[i]][2]);
            temp_cell_coords[i] =
                std::make_pair(p0, in_top[idx[i]].vdw_radius().value_or(1.5));
        }
        output_cells.push_back(temp_cell_coords);
    }

    return;
}

// NDD Specific function for PDB input. Hi precision method
void ndd_read_PDB_get_cells(std::string const &filename,
    NDD_IVector const &in_void_cells_indices,
    const std::vector<int> &include_CH_atoms, NDD_Vector &output_cells,
    Triang_Vector &CH_triangs) {

    // Read molecule
    chemfiles::Trajectory in_traj(filename);
    chemfiles::Frame in_frame = in_traj.read();
    auto in_xyz = in_frame.positions();
    chemfiles::Topology in_top = in_frame.topology();
    VertexInfo vi1;
    NDD_Element temp_cell_coords;
    std::vector<CPoint> incl_area_points;

    for (NDD_IElement const &idx : in_void_cells_indices) {
        // Iterate over each cell of the current pocket
        for (std::size_t i = 0; i < 4; ++i) {
            // Iterate over each vertex of the current cell
            CPoint const p0 =
                CPoint(in_xyz[idx[i]][0], in_xyz[idx[i]][1], in_xyz[idx[i]][2]);
            temp_cell_coords[i] =
                std::make_pair(p0, in_top[idx[i]].vdw_radius().value_or(1.5));
        }
        output_cells.push_back(temp_cell_coords);
    }

    // Actualize convex hull of the include area.
    for (auto const &each : include_CH_atoms) {
        incl_area_points.push_back(
            CPoint(in_xyz[each][0], in_xyz[each][1], in_xyz[each][2]));
    }

    Polyhedron CH;
    CGAL::convex_hull_3(incl_area_points.begin(), incl_area_points.end(), CH);
    P_Facet_const_iterator f_end = CH.facets_end();
    for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
         ++f_ite) {
        // Fix around the weirdest CGAL bug.
        P_Halfedge_around_facet_const_circulator he_ite = f_ite->facet_begin();
        auto const he_ite_0 = he_ite++;
        auto const he_ite_1 = he_ite++;
        auto const he_ite_2 = he_ite;

        CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
            he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
    }

    return;
}

// On-site NDD.
void ndd(Molecule const &protein, Cavity const &hueco, ConvexHull const &CH,
    IncludedAreaOptions const &IA_opts, NDDOptions const &NDD_opts,
    std::string const &pdb_filename) {

    std::vector<double> output_volumes;
    std::vector<int> all_indices;

    // NDD_IVector const cells_indices = get_vertices(cavity_void_cells);

    Modes const modos(NDD_opts._modes_ndd_filename);

    // if (CH._included_resis.size() != 0) {
    //     for (size_t j = 0; j < modos._j; ++j) {
    //         Molecule prote_ndd = protein;
    //         for (const auto nres : CH._included_resis) {
    //             // nres is 1-indexed.
    //             int const nres_3 = nres * 3 - 1;
    //             prote_ndd._data[nres].first +
    //                 CVector(modos._evectors[j][nres_3],
    //                     modos._evectors[j][nres_3 + 1],
    //                     modos._evectors[j][nres_3 + 2]);
    //         }
    //         ANA::ConvexHull const CH_ndd =
    //             create_convex_hull(prote_ndd, IA_opts);

    //         Cavity hueco_ndd;
    //         hueco_ndd._all_cells.reserve(
    //             hueco._inner_cells.size() + hueco._outer_cells.size());
    //         for (const auto &cell : hueco._inner_cells) {
    //             auto const p0_info = cell->vertex(0)->info();
    //             auto const p1_info = cell->vertex(1)->info();
    //             auto const p2_info = cell->vertex(2)->info();
    //             auto const p3_info = cell->vertex(3)->info();

    //             // _resn is 1-indexed.
    //             int const nres_3 = p0_info._resn * 3 - 1;
    //             cell->vertex(0)->point() +
    //                 CVector(modos._evectors[j][nres_3],
    //                     modos._evectors[j][nres_3 + 1],
    //                     modos._evectors[j][nres_3 + 2]);
    //         }
    //     }
    // } else if (CH._included_atoms.size() != 0) {
    // } else {
    // }

    // ANA::NDD::ndd_write_out_file(output_volumes, out_file);

    return;
}

// Perform Non Delaunay Dynamics.
void ndd_nondelaunay_dynamics_old(NA_Vector const &cavity_void_cells,
    std::string const &pdb_list, bool const precision,
    const std::vector<int> include_CH_atoms, std::string const &out_file) {
    std::vector<double> output_volumes;
    std::vector<int> all_indices;

    NDD_IVector cells_indices;
    std::ifstream input_pdbs_filename(pdb_list);
    std::string pdb_filename;

    if (input_pdbs_filename.is_open()) {
        // Get the atoms involved in each pocket.
        // new_cells_coordinates and cells_indices have similar structure.
        // They are vectors of vectors of arrays of 4 elements. Each array
        // refers
        // to a cell. Each vector of arrays refer to a pocket and all of
        // these pockets (vectors of arrays) are compiled into a vector.

        auto cells_indices = get_vertices(cavity_void_cells);

        if (precision == 1) {
            while (std::getline(input_pdbs_filename, pdb_filename)) {
                NDD_Vector init_cells_coordinates, new_cells_coordinates,
                    cavity_void_coords, cavity_intersecting_coords;
                Triang_Vector CH_triangs;
                std::vector<std::array<bool, 4>> intersecting_bool;
                std::vector<int> intersecting_total;
                Poly_Vector border_poly;

                // Get cell coordinates and the new include area triangles
                ANA::NDD::ndd_read_PDB_get_cells(pdb_filename, cells_indices,
                    include_CH_atoms, new_cells_coordinates, CH_triangs);
                // Separate fully inside cells and the intersecting ones.
                // Get the volume of this last group.
                ndd_discard_CH_0(new_cells_coordinates, CH_triangs,
                    cavity_void_coords, cavity_intersecting_coords,
                    intersecting_bool, intersecting_total);
                double poly_vol =
                    ndd_discard_CH_1(cavity_intersecting_coords, CH_triangs,
                        intersecting_bool, intersecting_total, border_poly);
                // Store result
                output_volumes.push_back(
                    poly_vol + ndd_get_void_volume(cavity_void_coords));
            }
        } else {
            while (std::getline(input_pdbs_filename, pdb_filename)) {
                NDD_Vector init_cells_coordinates, new_cells_coordinates,
                    cavity_void_coords;
                // Get cell coordinates
                ANA::NDD::ndd_read_PDB_get_cells(
                    pdb_filename, cells_indices, new_cells_coordinates);
                // Store result
                output_volumes.push_back(
                    ndd_get_void_volume(new_cells_coordinates));
            }
        }

        ANA::NDD::ndd_write_out_file(output_volumes, out_file);
    } else
        std::cerr << "Unable to open " << pdb_list << " for NDD" << '\n';

    return;
}

// Get the indices of the atoms involved in the given cells
NDD_IVector get_vertices(NA_Vector const &cavity_void_cells) {

    NDD_IVector cells_indices;
    cells_indices.reserve(cavity_void_cells.size() * 4);

    for (Finite_cells_iterator const ac_ite : cavity_void_cells) {
        NDD_IElement temp {ac_ite->vertex(0)->info()._index,
            ac_ite->vertex(1)->info()._index, ac_ite->vertex(2)->info()._index,
            ac_ite->vertex(3)->info()._index};
        cells_indices.push_back(std::move(temp));
    }

    return cells_indices;
}

// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(CPoint const &p_0, CPoint const &p_1,
    CPoint const &p_2, CPoint const &p_3, double const radius) {
    // get 1st point of the mini tetrahedron
    CVector vec_1 = p_1 - p_0;
    vec_1 =
        vec_1 / (std::sqrt(CGAL::to_double(vec_1.squared_length()))) * radius;
    CPoint point_1 = p_0 + vec_1;
    // get 2nd point of the mini tetrahedron
    CVector vec_2 = p_2 - p_0;
    vec_2 =
        vec_2 / (std::sqrt(CGAL::to_double(vec_2.squared_length()))) * radius;
    CPoint point_2 = p_0 + vec_2;
    // get 3rd point of the mini tetrahedron
    CVector vec_3 = p_3 - p_0;
    vec_3 =
        vec_3 / (std::sqrt(CGAL::to_double(vec_3.squared_length()))) * radius;
    CPoint point_3 = p_0 + vec_3;

    // Now, get the volume of the sphere slice
    // get the normal vector
    CVector plane_vec_1 = p_2 - p_1;
    CVector plane_vec_2 = p_3 - p_2;
    CVector plane_normal = CGAL::cross_product(plane_vec_1, plane_vec_2);
    // normalize the normal vector
    plane_normal = plane_normal /
        (std::sqrt(CGAL::to_double(plane_normal.squared_length())));
    // get the distance between the delaunay vertex and the plane formed by
    // p_1, p_2 and p_3
    double dist_to_plane = CGAL::to_double(plane_normal * vec_1);
    double h = radius - dist_to_plane;
    double spherical_cap_volume =
        1 / 3 * M_PI * std::pow(h, 2) * (3 * radius - h);
    double mini_tetrahedron_vol =
        CGAL::to_double(CGAL::volume(p_0, point_1, point_2, point_3));

    // Add the 2 volumes that represent the space occupied by the atom with
    // coordinates p0
    double volume =
        std::abs(mini_tetrahedron_vol) + std::abs(spherical_cap_volume);
    if (isnan(volume)) {
        volume = 0;
        //    std::cerr << "Warning: sphere_sector_vol gave NaN value. Void
        //    calculation "
        //                 "keeps going."
        //              << '\n';
    }

    return volume;
}

// Calc volume of the input cells. Reedited for array container.
double ndd_get_void_volume(NDD_Vector const &cavity_void_cells) {

    double current_cell_vol, volume = 0;

    for (const NDD_Element &cell : cavity_void_cells) {
        // Iterate over each cell of and get the total volume of the
        // tetrahedron
        current_cell_vol = CGAL::to_double(CGAL::volume(
            cell[0].first, cell[1].first, cell[2].first, cell[3].first));
        // Substract the volume filled by the 4 atoms in the vtces
        current_cell_vol = ndd_refine_cell_volume(current_cell_vol, cell);
        volume = volume + current_cell_vol;
    }

    return volume;
}
// Discard cells without a vertex inside the specified convex hull. Hi
// precision. NDD version
void ndd_discard_CH_0(NDD_Vector const &in_coords,
    Triang_Vector const &CH_triangs, NDD_Vector &out_coords,
    NDD_Vector &out_intersecting_coords,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &intersecting_total) {

    std::vector<CVector> CH_normals;
    std::vector<CPoint> CH_vtces;
    // CTriangle normals point inwards. Only inside points will give a
    // positive dot product against all normals
    for (auto const &triangle : CH_triangs) {
        CVector v1 = triangle.vertex(1) - triangle.vertex(0);
        CVector v2 = triangle.vertex(2) - triangle.vertex(1);
        CVector normal = CGAL::cross_product(v2, v1);
        normal = normal / std::sqrt(CGAL::to_double(normal.squared_length()));
        CH_normals.push_back(normal);
        CH_vtces.push_back(triangle.vertex(1));
    }

    // Now discard outside cells
    for (auto const &ndd_array : in_coords) {
        std::array<bool, 4> vtx_inside_bool = {false, false, false, false};
        int total = 0;
        // Get nbr of vertices that lie outside the include area.
        for (std::size_t i = 0; i <= 3; ++i) {
            CPoint test_point(ndd_array[i].first);

            for (std::size_t j = 0; j < CH_vtces.size(); j++) {
                CVector test_vtor = test_point - CH_vtces[j];
                test_vtor = test_vtor /
                    std::sqrt(CGAL::to_double(test_vtor.squared_length()));
                double test_dot_pdt =
                    CGAL::to_double(test_vtor * CH_normals[j]);

                if (test_dot_pdt < 0) {
                    vtx_inside_bool[i] = true;
                    ++total;
                    break;
                }
            }
        }
        if (total == 0) {
            // cell is entirely contained in the included area
            out_coords.push_back(ndd_array);
        } else if (total == 4) {
            // cell is outside the included area
            continue;
        } else {
            // cell instersects the included area
            out_intersecting_coords.push_back(ndd_array);
            intersecting_bool.push_back(vtx_inside_bool);
            intersecting_total.push_back(total);
        }
    }

    return;
}

// Discard parts of cells outside the specified triangulation using
// intersecitons. NDD version
double ndd_discard_CH_1(NDD_Vector const &in_intersecting_coords,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<int> &intersecting_total, Poly_Vector &border_poly) {
    // Use this vector to obtain the indices of the vtces used to form the
    // intersecting segments
    std::vector<std::size_t> v = {0, 1, 2, 3};
    double volume = 0;

    for (std::size_t each = 0; each < in_intersecting_coords.size(); ++each) {

        switch (intersecting_total[each]) {
        case 3: {
            for (std::size_t i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    // Got the handles of the included_area cell(s) that
                    // contain vertices of the intersecting cell 1 vertex
                    // inside included area
                    std::vector<CPoint> i_points;

                    // Get the points of the vtces of the intersecting cell
                    // p0 is the point inside the included area
                    CPoint p0(in_intersecting_coords[each][i].first);
                    double const VdW_radius_0 =
                        double(in_intersecting_coords[each][i].second);
                    CPoint p1(
                        in_intersecting_coords[each][get_i_not_equal(v, i, 1)]
                            .first);
                    CPoint p2(
                        in_intersecting_coords[each][get_i_not_equal(v, i, 2)]
                            .first);
                    CPoint p3(
                        in_intersecting_coords[each][get_i_not_equal(v, i, 3)]
                            .first);

                    // Get the intersecting segments
                    std::vector<Segment> segments = {
                        Segment(p0, p1), Segment(p0, p2), Segment(p0, p3)};

                    // Get the intersections
                    for (auto const &s : segments) {
                        for (auto const &t : CH_triangs) {
                            Object inter_obj = CGAL::intersection(s, t);
                            if (CPoint const *inter_point =
                                    CGAL::object_cast<CPoint>(&inter_obj)) {
                                i_points.push_back(*inter_point);
                                break;
                            }
                        }
                    }

                    // Check if included area and intersecting cell are
                    // sharing vertices, creating degeneracies.
                    bool degeneracies_bool = false;
                    for (std::size_t i = 0; i < i_points.size() - 1; i++) {
                        for (std::size_t k = i + 1; k < i_points.size(); k++) {
                            if (i_points[i] == i_points[k]) {
                                i_points[k] = i_points[k] +
                                    CVector(0.01 * i, 0.01 * k, 0.01);
                                degeneracies_bool = true;
                            }
                        }
                    }
                    if (degeneracies_bool == true) {
                        p0 = p0 - CVector(0.01, 0.01, 0.01);
                    }

                    // Construct the polyhedron and get its volume.
                    Polyhedron CH;
                    CH.make_tetrahedron(
                        p0, i_points[0], i_points[1], i_points[2]);
                    border_poly.push_back(CH);
                    // Add the volume of this polyhedron.
                    volume = volume +
                        std::abs(CGAL::to_double(CGAL::volume(
                            p0, i_points[0], i_points[1], i_points[2])));

                    // Substract the volume of the sphere sector.
                    volume = volume -
                        sphere_sector_vol(p0, i_points[0], i_points[1],
                            i_points[2], VdW_radius_0);
                }
            }
            break;
        }
        case 2: {
            for (std::size_t i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    for (std::size_t j = i + 1; j < 4; j++) {
                        if (!intersecting_bool[each][j]) {
                            // 2 vertices inside included area
                            std::vector<CPoint> i_points;
                            // Get the points of the vtces of the
                            // intersecting cell
                            CPoint p0(in_intersecting_coords[each][i].first);
                            CPoint p1(in_intersecting_coords[each][j].first);
                            // p0 p1 are the points inside the included area
                            double const VdW_radius_0 =
                                double(in_intersecting_coords[each][i].second);
                            double const VdW_radius_1 =
                                double(in_intersecting_coords[each][j].second);
                            std::vector<std::size_t> query_vec = {i, j};
                            CPoint p2(
                                in_intersecting_coords[each][get_i_not_equal(v,
                                                                 query_vec, 1)]
                                    .first);
                            CPoint p3(
                                in_intersecting_coords[each][get_i_not_equal(v,
                                                                 query_vec, 2)]
                                    .first);

                            // Get the intersecting segments
                            std::vector<Segment> segments = {Segment(p0, p2),
                                Segment(p0, p3), Segment(p1, p2),
                                Segment(p1, p3)};

                            // Get the intersections
                            for (auto const &s : segments) {
                                for (auto const &t : CH_triangs) {
                                    Object inter_obj = CGAL::intersection(s, t);
                                    if (CPoint const *inter_point =
                                            CGAL::object_cast<CPoint>(
                                                &inter_obj)) {
                                        i_points.push_back(*inter_point);
                                        break;
                                    }
                                }
                            }

                            // Check if included area and intersecting cell
                            // are sharing vertices, creating degeneracies.
                            bool degeneracies_bool = false;
                            for (std::size_t i = 0; i < i_points.size() - 1;
                                 i++) {
                                for (std::size_t k = i + 1; k < i_points.size();
                                     k++) {
                                    if (i_points[k] == i_points[k]) {
                                        i_points[k] = i_points[k] +
                                            CVector(0.01, 0.01, 0.01);
                                        degeneracies_bool = true;
                                    }
                                }
                            }
                            if (degeneracies_bool == true) {
                                p0 = p0 - CVector(0.01, 0.01, 0.01);
                                p1 = p1 - CVector(0.01, -0.01, 0.01);
                            }

                            // Construct the polyhedron and get its volume
                            Polyhedron CH;
                            CH.make_tetrahedron(
                                p0, i_points[0], i_points[1], p1);
                            CH.make_tetrahedron(
                                p1, i_points[0], i_points[1], i_points[2]);
                            CH.make_tetrahedron(
                                p1, i_points[1], i_points[2], i_points[3]);
                            border_poly.push_back(CH);
                            // Add the volume of this polyhedron
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(
                                    p0, i_points[0], i_points[1], p1)));
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(p1,
                                    i_points[0], i_points[1], i_points[2])));
                            volume = volume +
                                std::abs(CGAL::to_double(CGAL::volume(p1,
                                    i_points[1], i_points[2], i_points[3])));

                            // Substract the volume of the sphere sector
                            volume = volume -
                                sphere_sector_vol(p0, i_points[0], i_points[1],
                                    p1, VdW_radius_0);
                            volume = volume -
                                sphere_sector_vol(p1, i_points[0], i_points[1],
                                    i_points[2], VdW_radius_1);
                            volume = volume -
                                sphere_sector_vol(p1, i_points[1], i_points[2],
                                    i_points[3], VdW_radius_1);
                        }
                    }
                }
            }
            break;
        }
        case 1: {
            for (std::size_t i = 0; i < 4; i++) {
                if (!intersecting_bool[each][i]) {
                    for (std::size_t j = i + 1; j < 4; j++) {
                        if (!intersecting_bool[each][j]) {
                            for (std::size_t k = j + 1; k < 4; k++) {
                                if (!intersecting_bool[each][k]) {
                                    // 3 vertices inside the included area
                                    std::vector<CPoint> i_points;
                                    // Get the points of the vtces of the
                                    // intersecting cell
                                    CPoint p0(
                                        in_intersecting_coords[each][i].first);
                                    CPoint p1(
                                        in_intersecting_coords[each][j].first);
                                    CPoint p2(
                                        in_intersecting_coords[each][k].first);
                                    // p0 p1 p3 are the points inside the
                                    // included area
                                    double const VdW_radius_0 = double(
                                        in_intersecting_coords[each][i].second);
                                    double const VdW_radius_1 = double(
                                        in_intersecting_coords[each][j].second);
                                    double const VdW_radius_2 = double(
                                        in_intersecting_coords[each][k].second);
                                    std::vector<std::size_t> query_vec = {
                                        i, j, k};
                                    CPoint p3(in_intersecting_coords
                                                  [each][get_i_not_equal(
                                                             v, query_vec, 1)]
                                                      .first);

                                    // Get the intersecting segments
                                    std::vector<Segment> segments = {
                                        Segment(p0, p3), Segment(p1, p3),
                                        Segment(p2, p3)};

                                    // Get the intersections
                                    for (auto const &s : segments) {
                                        for (auto const &t : CH_triangs) {
                                            Object inter_obj =
                                                CGAL::intersection(s, t);
                                            if (CPoint const *inter_point =
                                                    CGAL::object_cast<CPoint>(
                                                        &inter_obj)) {
                                                i_points.push_back(
                                                    *inter_point);
                                                break;
                                            }
                                        }
                                    }

                                    // Check if included area and
                                    // intersecting cell are sharing
                                    // vertices, creating degeneracies.
                                    bool degeneracies_bool = false;
                                    for (std::size_t i = 0;
                                         i < i_points.size() - 1; i++) {
                                        for (std::size_t k = i + 1;
                                             k < i_points.size(); k++) {
                                            if (i_points[k] == i_points[k]) {
                                                i_points[k] = i_points[k] +
                                                    CVector(0.01, 0.01, 0.01);
                                                degeneracies_bool = true;
                                            }
                                        }
                                    }
                                    if (degeneracies_bool == true) {
                                        p0 = p0 - CVector(0.01, 0.01, 0.01);
                                        p1 = p1 - CVector(0.01, -0.01, 0.01);
                                        p2 = p2 - CVector(0.01, -0.01, -0.01);
                                    }

                                    // Construct the polyhedron and get its
                                    // volume
                                    Polyhedron CH;
                                    CH.make_tetrahedron(
                                        p0, p1, p2, i_points[0]);
                                    CH.make_tetrahedron(
                                        p0, p2, i_points[0], i_points[1]);
                                    CH.make_tetrahedron(p2, i_points[0],
                                        i_points[1], i_points[2]);
                                    border_poly.push_back(CH);
                                    // Add the volume of this polyhedron
                                    volume = volume +
                                        std::abs(CGAL::to_double(CGAL::volume(
                                            p0, p1, p2, i_points[0])));
                                    volume = volume +
                                        std::abs(CGAL::to_double(CGAL::volume(
                                            p0, p2, i_points[0], i_points[1])));
                                    volume = volume +
                                        std::abs(CGAL::to_double(
                                            CGAL::volume(p2, i_points[0],
                                                i_points[1], i_points[2])));

                                    // Substract the volume of the sphere
                                    // sector 1st tetra
                                    volume = volume -
                                        sphere_sector_vol(p0, p1, p2,
                                            i_points[0], VdW_radius_0);
                                    volume = volume -
                                        sphere_sector_vol(p1, p2, p0,
                                            i_points[0], VdW_radius_1);
                                    volume = volume -
                                        sphere_sector_vol(p2, p0, p1,
                                            i_points[0], VdW_radius_2);
                                    // 2nd tetra
                                    volume = volume -
                                        sphere_sector_vol(p0, p2, i_points[0],
                                            i_points[1], VdW_radius_0);
                                    volume = volume -
                                        sphere_sector_vol(p2, p0, i_points[0],
                                            i_points[1], VdW_radius_2);
                                    // 3rd tetra
                                    volume = volume -
                                        sphere_sector_vol(p2, i_points[0],
                                            i_points[1], i_points[2],
                                            VdW_radius_2);
                                }
                            }
                        }
                    }
                }
            }
            break;
        }
        }
    }

    return volume;
}

} // namespace ANA::NDD