#include <ANA/Write.hpp>

namespace ANA {

std::ofstream out_vol_stream;

// Draw pockets in .PDB format. Static version
void draw_raw_PDB(NA_Vector const &list_of_pockets, const Poly_Vector &polys,
    std::string const &out_filename) {

    int atom_cnt = 4;
    // Write PDB header.
    chemfiles::Topology ana_void_top;
    chemfiles::Frame ana_void_frame;
    auto out_traj = chemfiles::Trajectory(out_filename, 'w');

    for (Finite_cells_iterator const cell : list_of_pockets) {
        // Add the 4 vertices to the topology and their coords to the frame.
        ana_void_top.add_atom(chemfiles::Atom("N", "N0"));
        ana_void_frame.add_atom(chemfiles::Atom("N"),
            chemfiles::Vector3D(CGAL::to_double(cell->vertex(0)->point().x()),
                CGAL::to_double(cell->vertex(0)->point().y()),
                CGAL::to_double(cell->vertex(0)->point().z())));

        ana_void_top.add_atom(chemfiles::Atom("N", "N1"));
        ana_void_frame.add_atom(chemfiles::Atom("N"),
            chemfiles::Vector3D(CGAL::to_double(cell->vertex(1)->point().x()),
                CGAL::to_double(cell->vertex(1)->point().y()),
                CGAL::to_double(cell->vertex(1)->point().z())));

        ana_void_top.add_atom(chemfiles::Atom("N", "N2"));
        ana_void_frame.add_atom(chemfiles::Atom("N"),
            chemfiles::Vector3D(CGAL::to_double(cell->vertex(2)->point().x()),
                CGAL::to_double(cell->vertex(2)->point().y()),
                CGAL::to_double(cell->vertex(2)->point().z())));

        ana_void_top.add_atom(chemfiles::Atom("N", "N3"));
        ana_void_frame.add_atom(chemfiles::Atom("N"),
            chemfiles::Vector3D(CGAL::to_double(cell->vertex(3)->point().x()),
                CGAL::to_double(cell->vertex(3)->point().y()),
                CGAL::to_double(cell->vertex(3)->point().z())));

        // Connect the vertices
        connect_cell_residue(ana_void_top, atom_cnt);
        ana_void_top.add_residue(make_cell_residue(atom_cnt * 0.25, atom_cnt));
        atom_cnt += 4;
    }

    // Now, the polyhedron of the intersecting cells, if there are any.
    int res_cnt = atom_cnt * 0.25;
    atom_cnt -= 4;
    for (auto const &polyhedron : polys) {
        P_Facet_const_iterator fc_end = polyhedron.facets_end();
        int first_atom = atom_cnt;

        for (P_Facet_const_iterator fc_ite = polyhedron.facets_begin();
             fc_ite != fc_end; ++fc_ite) {
            P_Halfedge_around_facet_const_circulator he_circ =
                fc_ite->facet_begin();
            do {
                // Add the 3 vertices of each polyhedron.
                ana_void_top.add_atom(chemfiles::Atom("N", "NB"));
                ana_void_frame.add_atom(chemfiles::Atom("N"),
                    chemfiles::Vector3D(
                        CGAL::to_double(he_circ->vertex()->point().x()),
                        CGAL::to_double(he_circ->vertex()->point().y()),
                        CGAL::to_double(he_circ->vertex()->point().z())));
            } while (++he_circ != fc_ite->facet_begin());

            // Connect the vertices.
            connect_facet(ana_void_top, atom_cnt);
            atom_cnt += 3;
        }
        // Show the polyhedron as a residue
        ana_void_top.add_residue(
            make_polyhedron_residue(res_cnt, first_atom, atom_cnt));
        ++res_cnt;
    }

    // Applying the topology to the frame will set the atom names and bonds
    ana_void_frame.set_topology(ana_void_top);
    // Write
    out_traj.write(ana_void_frame);
    return;
}

// // Draw pockets in .PDB format. MD version
// void draw_raw_PDB(NA_Vector const &pocket, const Poly_Vector &polys,
//     std::string out_filename) {

//     int atom_cnt = 4;
//     // Write PDB header
//     chemfiles::Topology ana_void_top;
//     chemfiles::Frame ana_void_frame;
//     auto out_traj = chemfiles::Trajectory(out_filename, 'a');

//     for (Finite_cells_iterator const cell : pocket) {
//         // Add the 4 vertices to the topology and their cords to the frame
//         ana_void_top.add_atom(chemfiles::Atom("N", "N0"));
//         ana_void_frame.add_atom(chemfiles::Atom("N"),
//             chemfiles::Vector3D(CGAL::to_double(cell->vertex(0)->point().x()),
//                 CGAL::to_double(cell->vertex(0)->point().y()),
//                 CGAL::to_double(cell->vertex(0)->point().z())));

//         ana_void_top.add_atom(chemfiles::Atom("N", "N1"));
//         ana_void_frame.add_atom(chemfiles::Atom("N"),
//             chemfiles::Vector3D(CGAL::to_double(cell->vertex(1)->point().x()),
//                 CGAL::to_double(cell->vertex(1)->point().y()),
//                 CGAL::to_double(cell->vertex(1)->point().z())));

//         ana_void_top.add_atom(chemfiles::Atom("N", "N2"));
//         ana_void_frame.add_atom(chemfiles::Atom("N"),
//             chemfiles::Vector3D(CGAL::to_double(cell->vertex(2)->point().x()),
//                 CGAL::to_double(cell->vertex(2)->point().y()),
//                 CGAL::to_double(cell->vertex(2)->point().z())));

//         ana_void_top.add_atom(chemfiles::Atom("N", "N3"));
//         ana_void_frame.add_atom(chemfiles::Atom("N"),
//             chemfiles::Vector3D(CGAL::to_double(cell->vertex(3)->point().x()),
//                 CGAL::to_double(cell->vertex(3)->point().y()),
//                 CGAL::to_double(cell->vertex(3)->point().z())));

//         // Connect the vertices
//         connect_cell_residue(ana_void_top, atom_cnt);
//         ana_void_top.add_residue(make_cell_residue(atom_cnt * 0.25,
//         atom_cnt)); atom_cnt += 4;
//     }

//     // Now, the polyhedron of the intersecting cells.
//     int res_cnt = atom_cnt * 0.25;
//     atom_cnt -= 4;
//     for (auto const &polyhedron : polys) {
//         P_Facet_const_iterator fc_end = polyhedron.facets_end();
//         int first_atom = atom_cnt;

//         for (P_Facet_const_iterator fc_ite = polyhedron.facets_begin();
//              fc_ite != fc_end; ++fc_ite) {
//             P_Halfedge_around_facet_const_circulator he_circ =
//                 fc_ite->facet_begin();
//             do {
//                 // Add the 3 vertices of each polyhedron.
//                 ana_void_top.add_atom(chemfiles::Atom("N", "NB"));
//                 ana_void_frame.add_atom(chemfiles::Atom("N"),
//                     chemfiles::Vector3D(
//                         CGAL::to_double(he_circ->vertex()->point().x()),
//                         CGAL::to_double(he_circ->vertex()->point().y()),
//                         CGAL::to_double(he_circ->vertex()->point().z())));
//             } while (++he_circ != fc_ite->facet_begin());

//             // Connect the vertices.
//             connect_facet(ana_void_top, atom_cnt);
//             atom_cnt += 3;
//         }
//         // Show the polyhedron as a residue
//         ana_void_top.add_residue(
//             make_polyhedron_residue(res_cnt, first_atom, atom_cnt));
//         ++res_cnt;
//     }

//     // Applying the topology to the frame will set the atom names and bonds
//     ana_void_frame.set_topology(ana_void_top);
//     // Write
//     out_traj.write(ana_void_frame);
//     return;
// }

// Draw pockets in .PDB format. MD version, whole trajectory.
void draw_raw_PDB(const NDD_Matrix &list_of_pockets,
    const Poly_Matrix &list_of_polys, std::string &out_filename,
    int const max_atom_cnt) {

    out_filename.append(".pdb");
    auto out_traj = chemfiles::Trajectory(out_filename, 'w');

    for (size_t i = 0; i < list_of_pockets.size(); i++) {
        int atom_cnt = 4;
        // Write PDB header.
        chemfiles::Topology ana_void_top;
        chemfiles::Frame ana_void_frame;

        for (auto const &cell : list_of_pockets[i]) {
            // Add the 4 vertices to the topology and their cords to the frame
            ana_void_top.add_atom(chemfiles::Atom("N", "N0"));

            ana_void_frame.add_atom(chemfiles::Atom("N"),
                chemfiles::Vector3D(CGAL::to_double(cell[0].first.x()),
                    CGAL::to_double(cell[0].first.y()),
                    CGAL::to_double(cell[0].first.z())));

            ana_void_top.add_atom(chemfiles::Atom("N", "N1"));
            ana_void_frame.add_atom(chemfiles::Atom("N"),
                chemfiles::Vector3D(CGAL::to_double(cell[1].first.x()),
                    CGAL::to_double(cell[1].first.y()),
                    CGAL::to_double(cell[1].first.z())));

            ana_void_top.add_atom(chemfiles::Atom("N", "N2"));
            ana_void_frame.add_atom(chemfiles::Atom("N"),
                chemfiles::Vector3D(CGAL::to_double(cell[2].first.x()),
                    CGAL::to_double(cell[2].first.y()),
                    CGAL::to_double(cell[2].first.z())));

            ana_void_top.add_atom(chemfiles::Atom("N", "N3"));
            ana_void_frame.add_atom(chemfiles::Atom("N"),
                chemfiles::Vector3D(CGAL::to_double(cell[3].first.x()),
                    CGAL::to_double(cell[3].first.y()),
                    CGAL::to_double(cell[3].first.z())));

            // Connect the vertices
            connect_cell_residue(ana_void_top, atom_cnt);
            ana_void_top.add_residue(
                make_cell_residue(atom_cnt * 0.25, atom_cnt));
            atom_cnt += 4;
        }

        int res_cnt = atom_cnt * 0.25;

        if (list_of_polys.size() != 0) {
            atom_cnt -= 4;
            // Now, the polyhedron of the intersecting cells, if there are any.
            for (auto const &polyhedron : list_of_polys[i]) {
                P_Vertex_const_iterator v_end = polyhedron.vertices_end();
                P_Vertex_const_iterator v_beg = polyhedron.vertices_begin();
                int first_atom = atom_cnt;

                for (P_Vertex_const_iterator v_ite = v_beg; v_ite != v_end;
                     ++v_ite) {
                    // Add the 4/12 vertices of each polyhedron.
                    ana_void_top.add_atom(chemfiles::Atom("N", "NB"));
                    ana_void_frame.add_atom(chemfiles::Atom("N"),
                        chemfiles::Vector3D(CGAL::to_double(v_ite->point().x()),
                            CGAL::to_double(v_ite->point().y()),
                            CGAL::to_double(v_ite->point().z())));
                    atom_cnt += 1;
                }
                // Connect the vertices.
                if ((atom_cnt - first_atom) == 4) {
                    // These polyhedron is composed of 1 tetrahedron.
                    for (int j = first_atom; j < atom_cnt - 1; j++) {
                        for (int k = j + 1; k < atom_cnt; k++) {
                            ana_void_top.add_bond(j, k);
                        }
                    }
                } else if ((atom_cnt - first_atom) == 12) {
                    // These polyhedron is composed of 3 tetrahedrons.
                    int const low_treshold = first_atom + 4,
                              hi_treshold = first_atom + 8;

                    for (int j = first_atom; j < low_treshold - 1; j++) {
                        for (int k = j + 1; k < low_treshold; k++) {
                            ana_void_top.add_bond(j, k);
                        }
                    }
                    for (int j = low_treshold; j < hi_treshold - 1; j++) {
                        for (int k = j + 1; k < hi_treshold; k++) {
                            ana_void_top.add_bond(j, k);
                        }
                    }
                    for (int j = hi_treshold; j < atom_cnt - 1; j++) {
                        for (int k = j + 1; k < atom_cnt; k++) {
                            ana_void_top.add_bond(j, k);
                        }
                    }
                } else {
                    std::cerr << "whaaaa??" << '\n';
                }

                // Show the polyhedron as a residue
                ana_void_top.add_residue(
                    make_polyhedron_residue(res_cnt, first_atom, atom_cnt));
                ++res_cnt;
            }
        }
        // In order to make it VMD-friendly, every model needs to have the same
        // nbr of atoms.
        for (int j = atom_cnt; j < max_atom_cnt; j++) {
            ana_void_top.add_atom(chemfiles::Atom("N", "NE"));
            ana_void_frame.add_atom(
                chemfiles::Atom("N"), chemfiles::Vector3D(0., 0., 0.));
        }
        // Show the extra atoms as a residue.
        ana_void_top.add_residue(
            make_polyhedron_residue(res_cnt, atom_cnt, max_atom_cnt));
        ++res_cnt;

        // Applying the topology to the frame will set the atom names and bonds.
        ana_void_frame.set_topology(ana_void_top);

        // Write
        out_traj.write(ana_void_frame);
    }

    return;
}

// Draw pockets as dots in a PDB. MD version. Whole trajectory.
void draw_grid_pdb(const NDD_Matrix &list_of_pockets,
    const Poly_Matrix &list_of_polys,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<std::vector<int>> &list_intersecting_total,
    int const sphere_count, int const precision, std::string &out_filename) {

    // Write PDB header.
    int frame_nbr = list_of_pockets.size(), max_atom_cnt = 0;
    std::vector<int> atom_cnt_list(frame_nbr), res_cnt_list(frame_nbr);

    std::vector<chemfiles::Frame> list_ana_void_frame(frame_nbr);
    std::vector<chemfiles::Topology> list_ana_void_top(frame_nbr);
    out_filename.append(".pdb");
    auto out_traj = chemfiles::Trajectory(out_filename, 'w');

    for (int i = 0; i < frame_nbr; i++) {
        int res_cnt = 0;
        chemfiles::Topology ana_void_top;
        chemfiles::Frame ana_void_frame;
        atom_cnt_list[i] = make_grid_pdb(list_of_pockets[i], ana_void_top,
            ana_void_frame, sphere_count, res_cnt);

        if (precision == 1) {
            atom_cnt_list[i] = make_grid_pdb_polyhedrons(ana_void_top,
                ana_void_frame, in_vtces_radii, list_intersecting_total[i],
                list_of_polys[i], sphere_count, atom_cnt_list[i], res_cnt);
        }
        res_cnt_list[i] = res_cnt;

        // Add a bond so pymol does not assign a topology.
        ana_void_top.add_bond(1, 2);

        // Store the frame and topologies, for later addition of extra atoms.
        list_ana_void_frame[i] = std::move(ana_void_frame);
        list_ana_void_top[i] = std::move(ana_void_top);
        // Keep track of the atom count
        if (atom_cnt_list[i] > max_atom_cnt) {
            max_atom_cnt = atom_cnt_list[i];
        }
    }

    // Write.
    for (int i = 0; i < frame_nbr; i++) {
        // In order to make it VMD-friendly, every model needs to have the same
        // nbr of atoms. Adding extra atoms.
        for (int j = atom_cnt_list[i]; j < max_atom_cnt; j++) {
            list_ana_void_top[i].add_atom(chemfiles::Atom("N", "NE"));
            list_ana_void_frame[i].add_atom(
                chemfiles::Atom("N"), chemfiles::Vector3D(0., 0., 0.));
        }
        // Show the extra atoms as a residue.
        list_ana_void_top[i].add_residue(make_polyhedron_residue(
            res_cnt_list[i], atom_cnt_list[i], max_atom_cnt));

        // Applying the topology to the frame will set the atom names and bonds.
        list_ana_void_frame[i].set_topology(list_ana_void_top[i]);
        out_traj.write(list_ana_void_frame[i]);
    }
    return;
}

// Draw cells as dots in a PDB. Using "NDD_Vector" data structure.
int make_grid_pdb(NDD_Vector const &cells_to_draw,
    chemfiles::Topology &ana_void_top, chemfiles::Frame &ana_void_frame,
    int const sphere_count, int &res_cnt) {

    int n = sphere_count, atom_cnt_old, atom_cnt = 0;
    int const corner = (int)(sphere_count / 3);

    for (auto const &cell : cells_to_draw) {
        atom_cnt_old = atom_cnt;
        ++res_cnt;
        CPoint p0 = cell[0].first;
        CPoint p1 = cell[1].first;
        CPoint p2 = cell[2].first;
        CPoint p3 = cell[3].first;

        double const VdW_0 = cell[0].second;
        double const VdW_1 = cell[1].second;
        double const VdW_2 = cell[2].second;
        double const VdW_3 = cell[3].second;

        std::array<double, 3> p0_ = {CGAL::to_double(p0.x()),
            CGAL::to_double(p0.y()), CGAL::to_double(p0.z())};
        std::array<double, 3> p1_ = {CGAL::to_double(p1.x()),
            CGAL::to_double(p1.y()), CGAL::to_double(p1.z())};
        std::array<double, 3> p2_ = {CGAL::to_double(p2.x()),
            CGAL::to_double(p2.y()), CGAL::to_double(p2.z())};
        std::array<double, 3> p3_ = {CGAL::to_double(p3.x()),
            CGAL::to_double(p3.y()), CGAL::to_double(p3.z())};

        for (int i = 0; i < sphere_count; i++) {
            for (int j = 0; j <= sphere_count - i; j++) {
                for (int k = 0; k <= sphere_count - i - j; k++) {
                    std::array<double, 3> point;
                    int l = n - i - j - k;

                    point[0] =
                        (i * p0_[0] + j * p1_[0] + k * p2_[0] + l * p3_[0]) / n;
                    point[1] =
                        (i * p0_[1] + j * p1_[1] + k * p2_[1] + l * p3_[1]) / n;
                    point[2] =
                        (i * p0_[2] + j * p1_[2] + k * p2_[2] + l * p3_[2]) / n;

                    if (i <= corner || j <= corner || k <= corner) {
                        // Check if any sphere is in conflict with a vtx
                        double const dist0 =
                            std::sqrt(pow((point[0] - p0_[0]), 2.0) +
                                pow((point[1] - p0_[1]), 2.0) +
                                pow((point[2] - p0_[2]), 2.0));
                        if (dist0 < VdW_0) {
                            continue;
                        }
                        double const dist1 =
                            std::sqrt(pow((point[0] - p1_[0]), 2.0) +
                                pow((point[1] - p1_[1]), 2.0) +
                                pow((point[2] - p1_[2]), 2.0));
                        if (dist1 < VdW_1) {
                            continue;
                        }
                        double const dist2 =
                            std::sqrt(pow((point[0] - p2_[0]), 2.0) +
                                pow((point[1] - p2_[1]), 2.0) +
                                pow((point[2] - p2_[2]), 2.0));
                        if (dist2 < VdW_2) {
                            continue;
                        }
                        double const dist3 =
                            std::sqrt(pow((point[0] - p3_[0]), 2.0) +
                                pow((point[1] - p3_[1]), 2.0) +
                                pow((point[2] - p3_[2]), 2.0));
                        if (dist3 < VdW_3) {
                            continue;
                        }
                    }

                    ana_void_top.add_atom(chemfiles::Atom("H", "H"));
                    ana_void_frame.add_atom(chemfiles::Atom("H"),
                        chemfiles::Vector3D(point[0], point[1], point[2]));
                    ++atom_cnt;
                }
            }
        }
        // Collect all this points as a residue.
        ana_void_top.add_residue(
            make_grid_residue(atom_cnt_old, atom_cnt, res_cnt));
    }

    return atom_cnt;
}

// Draw a vector of polyhedrons as dots in a PDB.
int make_grid_pdb_polyhedrons(chemfiles::Topology &ana_void_top,
    chemfiles::Frame &ana_void_frame,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<int> &intersecting_total, const Poly_Vector &CH_vec,
    double const sphere_count, int atom_cnt, int &res_cnt) {

    int n = sphere_count, atom_cnt_old;
    int const corner = (int)(sphere_count / 3);

    for (std::size_t ii = 0; ii < CH_vec.size(); ++ii) {

        atom_cnt_old = atom_cnt;
        ++res_cnt;
        std::vector<CPoint> points;
        P_Vertex_const_iterator v_end = CH_vec[ii].vertices_end();
        for (P_Vertex_const_iterator v_ite = CH_vec[ii].vertices_begin();
             v_ite != v_end; ++v_ite) {
            points.push_back(v_ite->point());
        }

        switch (intersecting_total[ii]) {
        case 3: {
            // This polyhedron is a tetrahedron from the intersection of a cell
            // with
            // 1
            // vtx inside the included area
            std::array<double, 3> p0_ = {CGAL::to_double(points[0].x()),
                CGAL::to_double(points[0].y()), CGAL::to_double(points[0].z())};
            std::array<double, 3> p1_ = {CGAL::to_double(points[1].x()),
                CGAL::to_double(points[1].y()), CGAL::to_double(points[1].z())};
            std::array<double, 3> p2_ = {CGAL::to_double(points[2].x()),
                CGAL::to_double(points[2].y()), CGAL::to_double(points[2].z())};
            std::array<double, 3> p3_ = {CGAL::to_double(points[3].x()),
                CGAL::to_double(points[3].y()), CGAL::to_double(points[3].z())};
            for (int i = 0; i < sphere_count; i++) {
                for (int j = 0; j <= sphere_count - i; j++) {
                    for (int k = 0; k <= sphere_count - i - j; k++) {
                        std::array<double, 3> point;
                        int l = n - i - j - k;

                        point[0] = (i * p0_[0] + j * p1_[0] + k * p2_[0] +
                                       l * p3_[0]) /
                            n;
                        point[1] = (i * p0_[1] + j * p1_[1] + k * p2_[1] +
                                       l * p3_[1]) /
                            n;
                        point[2] = (i * p0_[2] + j * p1_[2] + k * p2_[2] +
                                       l * p3_[2]) /
                            n;

                        if (i <= corner || j <= corner || k <= corner) {
                            // Check if any sphere is in conflict with a vtx
                            double const dist0 =
                                std::sqrt(pow((point[0] - p0_[0]), 2.0) +
                                    pow((point[1] - p0_[1]), 2.0) +
                                    pow((point[2] - p0_[2]), 2.0));
                            if (dist0 < in_vtces_radii[ii][0]) {
                                continue;
                            }
                        }
                        ana_void_top.add_atom(chemfiles::Atom("H", "HB"));
                        ana_void_frame.add_atom(chemfiles::Atom("H"),
                            chemfiles::Vector3D(point[0], point[1], point[2]));
                        ++atom_cnt;
                    }
                }
            }
            break;
        }
        case 2: {
            // This polyhedron comprises 3 tetrahedrons from the intersection of
            // a
            // cell with 2 vtces inside the included area
            for (std::size_t jj = 0; jj < 3; ++jj) {
                std::array<double, 3> p0_ = {
                    CGAL::to_double(points[(jj * 4) + 0].x()),
                    CGAL::to_double(points[(jj * 4) + 0].y()),
                    CGAL::to_double(points[(jj * 4) + 0].z())};
                std::array<double, 3> p1_ = {
                    CGAL::to_double(points[(jj * 4) + 1].x()),
                    CGAL::to_double(points[(jj * 4) + 1].y()),
                    CGAL::to_double(points[(jj * 4) + 1].z())};
                std::array<double, 3> p2_ = {
                    CGAL::to_double(points[(jj * 4) + 2].x()),
                    CGAL::to_double(points[(jj * 4) + 2].y()),
                    CGAL::to_double(points[(jj * 4) + 2].z())};
                std::array<double, 3> p3_ = {
                    CGAL::to_double(points[(jj * 4) + 3].x()),
                    CGAL::to_double(points[(jj * 4) + 3].y()),
                    CGAL::to_double(points[(jj * 4) + 3].z())};
                for (int i = 0; i < sphere_count; i++) {
                    for (int j = 0; j <= sphere_count - i; j++) {
                        for (int k = 0; k <= sphere_count - i - j; k++) {
                            std::array<double, 3> point;
                            int l = n - i - j - k;

                            point[0] = (i * p0_[0] + j * p1_[0] + k * p2_[0] +
                                           l * p3_[0]) /
                                n;
                            point[1] = (i * p0_[1] + j * p1_[1] + k * p2_[1] +
                                           l * p3_[1]) /
                                n;
                            point[2] = (i * p0_[2] + j * p1_[2] + k * p2_[2] +
                                           l * p3_[2]) /
                                n;

                            if (i <= corner || j <= corner || k <= corner) {
                                // Check if any sphere is in conflict with a vtx
                                switch (jj) {
                                case 0: {
                                    // 1st and 4th vtces correspond to atoms
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][0]) {
                                        continue;
                                    }

                                    double const dist1 = std::sqrt(
                                        pow((point[0] - p3_[0]), 2.0) +
                                        pow((point[1] - p3_[1]), 2.0) +
                                        pow((point[2] - p3_[2]), 2.0));
                                    if (dist1 < in_vtces_radii[ii][1]) {
                                        continue;
                                    }
                                    break;
                                }
                                case 1: {
                                    // 1st vtx corresponds to an atom
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][1]) {
                                        continue;
                                    }
                                    [[fallthrough]];
                                }
                                case 2: {
                                    // 1st vtx corresponds to an atom
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][1]) {
                                        continue;
                                    }
                                    break;
                                }
                                }
                            }

                            ana_void_top.add_atom(chemfiles::Atom("H", "HB"));
                            ana_void_frame.add_atom(chemfiles::Atom("H"),
                                chemfiles::Vector3D(
                                    point[0], point[1], point[2]));
                            ++atom_cnt;
                        }
                    }
                }
            }
            break;
        }
        case 1: {
            // This polyhedron comprises 3 tetrahedrons from the intersection of
            // a
            // cell with 3 vtces inside the included area
            for (std::size_t jj = 0; jj < 3; ++jj) {
                std::array<double, 3> p0_ = {
                    CGAL::to_double(points[(jj * 4) + 0].x()),
                    CGAL::to_double(points[(jj * 4) + 0].y()),
                    CGAL::to_double(points[(jj * 4) + 0].z())};
                std::array<double, 3> p1_ = {
                    CGAL::to_double(points[(jj * 4) + 1].x()),
                    CGAL::to_double(points[(jj * 4) + 1].y()),
                    CGAL::to_double(points[(jj * 4) + 1].z())};
                std::array<double, 3> p2_ = {
                    CGAL::to_double(points[(jj * 4) + 2].x()),
                    CGAL::to_double(points[(jj * 4) + 2].y()),
                    CGAL::to_double(points[(jj * 4) + 2].z())};
                std::array<double, 3> p3_ = {
                    CGAL::to_double(points[(jj * 4) + 3].x()),
                    CGAL::to_double(points[(jj * 4) + 3].y()),
                    CGAL::to_double(points[(jj * 4) + 3].z())};
                for (int i = 0; i < sphere_count; i++) {
                    for (int j = 0; j <= sphere_count - i; j++) {
                        for (int k = 0; k <= sphere_count - i - j; k++) {
                            std::array<double, 3> point;
                            int l = n - i - j - k;

                            point[0] = (i * p0_[0] + j * p1_[0] + k * p2_[0] +
                                           l * p3_[0]) /
                                n;
                            point[1] = (i * p0_[1] + j * p1_[1] + k * p2_[1] +
                                           l * p3_[1]) /
                                n;
                            point[2] = (i * p0_[2] + j * p1_[2] + k * p2_[2] +
                                           l * p3_[2]) /
                                n;

                            if (i <= corner || j <= corner || k <= corner) {
                                // Check if any sphere is in conflict with a vtx
                                switch (jj) {
                                case 0: {
                                    // 1st, 2nd and 3rd vtces correspond to
                                    // atoms
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][0]) {
                                        continue;
                                    }

                                    double const dist1 = std::sqrt(
                                        pow((point[0] - p1_[0]), 2.0) +
                                        pow((point[1] - p1_[1]), 2.0) +
                                        pow((point[2] - p1_[2]), 2.0));
                                    if (dist1 < in_vtces_radii[ii][1]) {
                                        continue;
                                    }

                                    double const dist2 = std::sqrt(
                                        pow((point[0] - p2_[0]), 2.0) +
                                        pow((point[1] - p2_[1]), 2.0) +
                                        pow((point[2] - p2_[2]), 2.0));
                                    if (dist2 < in_vtces_radii[ii][2]) {
                                        continue;
                                    }
                                    break;
                                }
                                case 1: {
                                    // 1st, and 2nd vtces correspond to atoms
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][0]) {
                                        continue;
                                    }
                                    double const dist1 = std::sqrt(
                                        pow((point[0] - p1_[0]), 2.0) +
                                        pow((point[1] - p1_[1]), 2.0) +
                                        pow((point[2] - p1_[2]), 2.0));
                                    if (dist1 < in_vtces_radii[ii][2]) {
                                        continue;
                                    }
                                    break;
                                }
                                case 2: {
                                    // 1st vtx corresponds to an atom
                                    double const dist0 = std::sqrt(
                                        pow((point[0] - p0_[0]), 2.0) +
                                        pow((point[1] - p0_[1]), 2.0) +
                                        pow((point[2] - p0_[2]), 2.0));
                                    if (dist0 < in_vtces_radii[ii][2]) {
                                        continue;
                                    }
                                    break;
                                }
                                }
                            }

                            ana_void_top.add_atom(chemfiles::Atom("H", "HB"));
                            ana_void_frame.add_atom(chemfiles::Atom("H"),
                                chemfiles::Vector3D(
                                    point[0], point[1], point[2]));
                            ++atom_cnt;
                        }
                    }
                }
            }
            break;
        }
        }
        // Collect all this points as a residue
        ana_void_top.add_residue(
            make_grid_residue(atom_cnt_old, atom_cnt, res_cnt));
    }

    return atom_cnt;
}

// Draw pockets as points in a PDB. Using CGAL Cell data structure. Static
// version.
void draw_grid_pdb(NA_Vector const &pocket,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<int> &intersecting_total, const Poly_Vector &polys,
    std::string &out_filename, int const sphere_count, int const precision) {
    // Write PDB header.
    int res_cnt = 0;
    chemfiles::Topology ana_void_top;
    chemfiles::Frame ana_void_frame;
    auto out_traj = chemfiles::Trajectory(out_filename, 'w');

    auto atom_cnt = make_grid_pdb(
        pocket, ana_void_top, ana_void_frame, sphere_count, res_cnt);

    if (precision == 1) {
        make_grid_pdb_polyhedrons(ana_void_top, ana_void_frame, in_vtces_radii,
            intersecting_total, polys, sphere_count, atom_cnt, res_cnt);
    }

    // Applying the topology to the frame will set the atom names and bonds
    ana_void_frame.set_topology(ana_void_top);
    // Write
    out_traj.write(ana_void_frame);
    return;
}

// Draw cells as dots in a PDB. Using CGAL Cell data structure.
int make_grid_pdb(NA_Vector const &cells_to_draw,
    chemfiles::Topology &ana_void_top, chemfiles::Frame &ana_void_frame,
    int const sphere_count, int &res_cnt) {

    int n = sphere_count, atom_cnt_old, atom_cnt = 0;
    int const corner = (int)(sphere_count / 3);

    for (auto const &ac_ite : cells_to_draw) {
        atom_cnt_old = atom_cnt;
        ++res_cnt;
        CPoint p0 = ac_ite->vertex(0)->point();
        CPoint p1 = ac_ite->vertex(1)->point();
        CPoint p2 = ac_ite->vertex(2)->point();
        CPoint p3 = ac_ite->vertex(3)->point();

        double const VdW_0 = ac_ite->vertex(0)->info()._radius;
        double const VdW_1 = ac_ite->vertex(1)->info()._radius;
        double const VdW_2 = ac_ite->vertex(2)->info()._radius;
        double const VdW_3 = ac_ite->vertex(3)->info()._radius;

        std::array<double, 3> p0_ = {CGAL::to_double(p0.x()),
            CGAL::to_double(p0.y()), CGAL::to_double(p0.z())};
        std::array<double, 3> p1_ = {CGAL::to_double(p1.x()),
            CGAL::to_double(p1.y()), CGAL::to_double(p1.z())};
        std::array<double, 3> p2_ = {CGAL::to_double(p2.x()),
            CGAL::to_double(p2.y()), CGAL::to_double(p2.z())};
        std::array<double, 3> p3_ = {CGAL::to_double(p3.x()),
            CGAL::to_double(p3.y()), CGAL::to_double(p3.z())};

        for (int i = 0; i < sphere_count; i++) {
            for (int j = 0; j <= sphere_count - i; j++) {
                for (int k = 0; k <= sphere_count - i - j; k++) {
                    std::array<double, 3> point;
                    int l = n - i - j - k;

                    point[0] =
                        (i * p0_[0] + j * p1_[0] + k * p2_[0] + l * p3_[0]) / n;
                    point[1] =
                        (i * p0_[1] + j * p1_[1] + k * p2_[1] + l * p3_[1]) / n;
                    point[2] =
                        (i * p0_[2] + j * p1_[2] + k * p2_[2] + l * p3_[2]) / n;

                    if (i <= corner || j <= corner || k <= corner) {
                        // Check if any sphere is in conflict with a vtx
                        double const dist0 =
                            std::sqrt(pow((point[0] - p0_[0]), 2.0) +
                                pow((point[1] - p0_[1]), 2.0) +
                                pow((point[2] - p0_[2]), 2.0));
                        if (dist0 < VdW_0) {
                            continue;
                        }
                        double const dist1 =
                            std::sqrt(pow((point[0] - p1_[0]), 2.0) +
                                pow((point[1] - p1_[1]), 2.0) +
                                pow((point[2] - p1_[2]), 2.0));
                        if (dist1 < VdW_1) {
                            continue;
                        }
                        double const dist2 =
                            std::sqrt(pow((point[0] - p2_[0]), 2.0) +
                                pow((point[1] - p2_[1]), 2.0) +
                                pow((point[2] - p2_[2]), 2.0));
                        if (dist2 < VdW_2) {
                            continue;
                        }
                        double const dist3 =
                            std::sqrt(pow((point[0] - p3_[0]), 2.0) +
                                pow((point[1] - p3_[1]), 2.0) +
                                pow((point[2] - p3_[2]), 2.0));
                        if (dist3 < VdW_3) {
                            continue;
                        }
                    }

                    ana_void_top.add_atom(chemfiles::Atom("H", "H"));
                    ana_void_frame.add_atom(chemfiles::Atom("H"),
                        chemfiles::Vector3D(point[0], point[1], point[2]));
                    ++atom_cnt;
                }
            }
        }
        // Collect all this points as a residue
        ana_void_top.add_residue(
            make_grid_residue(atom_cnt_old, atom_cnt, res_cnt));
    }

    return atom_cnt;
}

/////////
//
////////

// Write .PDB header
inline void header_PDB(
    std::string const &out_pdb_filename, std::string const &in_pdb_filename) {
    [[maybe_unused]] auto borrame = out_pdb_filename.end();
    [[maybe_unused]] auto borrame1 = in_pdb_filename.end();
    return;
}

// Construct a residue object with 4 atoms starting at index cell_cnt
inline chemfiles::Residue make_cell_residue(
    int const cell_cnt, int const atom_cnt) {
    // Each cell will be considered as a separate residue
    chemfiles::Residue pocket_res("ANA", cell_cnt);
    // Create all cell "atoms"
    pocket_res.add_atom(atom_cnt - 4);
    pocket_res.add_atom(atom_cnt - 3);
    pocket_res.add_atom(atom_cnt - 2);
    pocket_res.add_atom(atom_cnt - 1);
    return pocket_res;
}

// Construct a residue object with atoms from 'first_atom' to 'atom_cnt'
inline chemfiles::Residue make_polyhedron_residue(
    int const res_cnt, int const first_atom, int const atom_cnt) {
    // Each polyhedron will be considered as a separate residue
    chemfiles::Residue pocket_res("ANA", res_cnt);

    for (int i = first_atom; i < atom_cnt; i++) {
        pocket_res.add_atom(i);
    }

    return pocket_res;
}

// Construct a residue object for grid output
inline chemfiles::Residue make_grid_residue(
    int const atom_cnt_old, int const atom_cnt, int const res_cnt) {
    // Each cell will be considered as a separate residue
    chemfiles::Residue pocket_res("ANA", res_cnt);
    // Create all cell "atoms"
    for (int i = atom_cnt_old; i < atom_cnt; ++i) {
        pocket_res.add_atom(i);
    }
    return pocket_res;
}

// Connect a residue object with the atoms starting at index cell_cnt
inline void connect_cell_residue(
    chemfiles::Topology &topology, int const atom_index) {
    topology.add_bond(atom_index - 4, atom_index - 3);
    topology.add_bond(atom_index - 4, atom_index - 2);
    topology.add_bond(atom_index - 4, atom_index - 1);
    topology.add_bond(atom_index - 3, atom_index - 2);
    topology.add_bond(atom_index - 3, atom_index - 1);
    topology.add_bond(atom_index - 2, atom_index - 1);
    return;
}

// Connect a facet of a polyhedron with the atoms starting at index atom_index
inline void connect_facet(chemfiles::Topology &topology, int const atom_index) {

    topology.add_bond(atom_index, atom_index + 1);
    topology.add_bond(atom_index, atom_index + 2);
    topology.add_bond(atom_index + 1, atom_index + 2);
    return;
}

//////////
//
/////////

// Draw a whole convex hull contained in a polyhedron. Static version.
void draw_CH(const Polyhedron &CH, std::string &out_filename) {

    P_Facet_const_iterator fc_ite_end = CH.facets_end();
    out_filename.append(".pdb");
    chemfiles::Trajectory out_trj(out_filename, 'w');
    chemfiles::Frame out_frm;
    chemfiles::Topology out_top;
    int atom_cnt = 0;

    for (P_Facet_const_iterator fc_ite = CH.facets_begin();
         fc_ite != fc_ite_end; ++fc_ite) {
        P_Halfedge_around_facet_const_circulator he_circ =
            fc_ite->facet_begin();
        do {
            out_top.add_atom(chemfiles::Atom("CH", "C"));
            out_frm.add_atom(chemfiles::Atom("C"),
                chemfiles::Vector3D(
                    CGAL::to_double(he_circ->vertex()->point().x()),
                    CGAL::to_double(he_circ->vertex()->point().y()),
                    CGAL::to_double(he_circ->vertex()->point().z())));
            ++atom_cnt;
        } while (++he_circ != fc_ite->facet_begin());
        out_top.add_bond(atom_cnt - 3, atom_cnt - 2);
        out_top.add_bond(atom_cnt - 3, atom_cnt - 1);
        out_top.add_bond(atom_cnt - 2, atom_cnt - 1);
    }
    out_frm.set_topology(out_top);
    out_trj.write(out_frm);

    return;
}

// Draw a whole convex hull contained in vector of triangles. Static version.
void draw_CH(Triang_Vector const &CH_triang, std::string &out_filename) {

    out_filename.append(".pdb");
    chemfiles::Trajectory out_trj(out_filename, 'w');
    chemfiles::Frame out_frm;
    chemfiles::Topology out_top;
    int atom_cnt = 0;

    for (auto const &triangle : CH_triang) {
        for (size_t i = 0; i < 3; i++) {
            out_top.add_atom(chemfiles::Atom("CH", "C"));
            out_frm.add_atom(chemfiles::Atom("C"),
                chemfiles::Vector3D(CGAL::to_double(triangle.vertex(i).x()),
                    CGAL::to_double(triangle.vertex(i).y()),
                    CGAL::to_double(triangle.vertex(i).z())));
            ++atom_cnt;
        }
        out_top.add_bond(atom_cnt - 3, atom_cnt - 2);
        out_top.add_bond(atom_cnt - 3, atom_cnt - 1);
        out_top.add_bond(atom_cnt - 2, atom_cnt - 1);
    }
    out_frm.set_topology(out_top);
    out_trj.write(out_frm);

    return;
}

// Draw a whole convex hull contained in vector of triangles. MD version.
void draw_CH(Triang_Vector const &CH_triang, chemfiles::Trajectory &out_traj) {

    chemfiles::Frame out_frm;
    chemfiles::Topology out_top;
    int atom_cnt = 0;

    for (auto const &triangle : CH_triang) {
        for (size_t i = 0; i < 3; i++) {
            out_top.add_atom(chemfiles::Atom("CH", "C"));
            out_frm.add_atom(chemfiles::Atom("C"),
                chemfiles::Vector3D(CGAL::to_double(triangle.vertex(i).x()),
                    CGAL::to_double(triangle.vertex(i).y()),
                    CGAL::to_double(triangle.vertex(i).z())));
            ++atom_cnt;
        }
        out_top.add_bond(atom_cnt - 3, atom_cnt - 2);
        out_top.add_bond(atom_cnt - 3, atom_cnt - 1);
        out_top.add_bond(atom_cnt - 2, atom_cnt - 1);
    }
    out_frm.set_topology(out_top);
    out_traj.write(out_frm);

    return;
}

// Write wall amino acids and atoms
void wall_atom_output(std::ofstream &wall_out, NA_Vector const &in_cells,
    NA_Vector const &in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    bool const requested_CH, int const precision, int const pock_cnt,
    int const frame_cnt, std::string const &list_wall_separator) {

    std::vector<int> wall_aa_idx, wall_atom_idx;
    std::vector<std::string> wall_aa_id;
    // Pocket ID
    std::string pock_out_filename = "pocket_";
    pock_out_filename.append(std::to_string(pock_cnt));

    // Get wall amino acids and atoms
    get_info_cell(in_cells, wall_aa_idx, wall_aa_id, wall_atom_idx);
    if (requested_CH && precision == 1) {
        get_info_cell(in_intersecting_cells, intersecting_bool, wall_aa_idx,
            wall_aa_id, wall_atom_idx);
    }

    write_wall_file(wall_out, pock_out_filename, wall_aa_idx, wall_aa_id,
        wall_atom_idx, frame_cnt, list_wall_separator);

    return;
}

// Write wall amino acids.
void wall_aa_output(std::ofstream &wall_out, NA_Vector const &in_cells,
    NA_Vector const &in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    bool const requested_CH, int const precision, int const pock_cnt,
    int const frame_cnt, std::string const &list_wall_separator) {

    std::vector<int> wall_aa_idx;
    std::vector<std::string> wall_aa_id;
    std::string pock_out_filename = "pocket_";
    pock_out_filename.append(std::to_string(pock_cnt));

    // Get wall amino acids.
    get_info_cell(in_cells, wall_aa_idx, wall_aa_id);
    if (requested_CH && (precision == 1)) {
        get_info_cell(
            in_intersecting_cells, intersecting_bool, wall_aa_idx, wall_aa_id);
    }

    // Write output.
    write_wall_file(wall_out, pock_out_filename, wall_aa_idx, wall_aa_id,
        frame_cnt, list_wall_separator);

    return;
}

// Get names and indices of participating atoms and amino acids
void get_info_cell(NA_Vector const &null_areas_vtor,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id,
    std::vector<int> &wall_atom_idx) {

    std::vector<int>::iterator a_idx;
    int atom_idx, a;
    VertexInfo v_info;

    for (auto const &cell : null_areas_vtor) {
        for (std::size_t i = 0; i < 3; i++) {
            // Get Info data structure
            v_info = cell->vertex(i)->info();
            // Get Atom Index
            atom_idx = v_info._index;
            // Get position of the new index
            a_idx = std::lower_bound(
                wall_atom_idx.begin(), wall_atom_idx.end(), atom_idx);

            if (a_idx == wall_atom_idx.end()) {
                // Store new atom (with largest index) info
                wall_atom_idx.push_back(v_info._index);
                wall_aa_idx.push_back(v_info._resn);
                wall_aa_id.push_back(v_info._resi);
                continue;
            } else if (atom_idx == *a_idx || atom_idx == *(++a_idx)) {
                // This atom is already present
                continue;
            }

            a = (--a_idx) - wall_atom_idx.begin();
            // Store new atom info
            wall_atom_idx.insert(a_idx, v_info._index);
            wall_aa_idx.insert(wall_aa_idx.begin() + a, v_info._resn);
            wall_aa_id.insert(wall_aa_id.begin() + a, v_info._resi);
        }
    }

    return;
}

// Get names and indices of participating atoms and amino acids from
// intersecting cells
void get_info_cell(NA_Vector const &cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id,
    std::vector<int> &wall_atom_idx) {

    std::vector<int>::iterator a_idx;
    int atom_idx, a, ii = 0;
    VertexInfo v_info;

    for (auto const &cell : cavity_intersecting_cells) {
        for (std::size_t i = 0; i < 3; i++) {
            if (intersecting_bool[ii][i]) {
                // This vtx is outside the included area
                continue;
            }
            // Get Info data structure
            v_info = cell->vertex(i)->info();
            // Get Atom Index
            atom_idx = v_info._index;
            // Get position of the new index
            a_idx = std::lower_bound(
                wall_atom_idx.begin(), wall_atom_idx.end(), atom_idx);

            if (a_idx == wall_atom_idx.end()) {
                // Store new atom (with largest index) info
                wall_atom_idx.push_back(v_info._index);
                wall_aa_idx.push_back(v_info._resn);
                wall_aa_id.push_back(v_info._resi);
                continue;
            } else if (atom_idx == *a_idx || atom_idx == *(++a_idx)) {
                // This atom is already present
                continue;
            }

            a = (--a_idx) - wall_atom_idx.begin();
            // Store new atom info
            wall_atom_idx.insert(a_idx, v_info._index);
            wall_aa_idx.insert(wall_aa_idx.begin() + a, v_info._resn);
            wall_aa_id.insert(wall_aa_id.begin() + a, v_info._resi);
        }
        ++ii;
    }

    return;
}

// Get names and indices of participating amino acids.
void get_info_cell(NA_Vector const &null_areas_vtor,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id) {

    std::vector<int>::iterator r_idx;
    int res_idx, r;
    VertexInfo v_info;

    for (auto const &cell : null_areas_vtor) {
        for (std::size_t i = 0; i < 3; i++) {
            // Get Info data structure
            v_info = cell->vertex(i)->info();
            // Get residue Index
            res_idx = v_info._resn;
            // Get position of the new index
            r_idx = std::lower_bound(
                wall_aa_idx.begin(), wall_aa_idx.end(), res_idx);

            if (r_idx == wall_aa_idx.end()) {
                // Store new residue (with largest index) info
                wall_aa_idx.push_back(res_idx);
                wall_aa_id.push_back(v_info._resi);
            } else if (res_idx == *r_idx || res_idx == *(++r_idx)) {
                // This residue is already present
                continue;
            } else {
                // Store new residue info
                r = (--r_idx) - wall_aa_idx.begin();
                wall_aa_idx.insert(r_idx, res_idx);
                wall_aa_id.insert(wall_aa_id.begin() + r, v_info._resi);
            }
        }
    }

    return;
}

// Get names and indices of participating amino acids from intersecting cells
void get_info_cell(NA_Vector const &cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id) {

    std::vector<int>::iterator r_idx;
    int res_idx, r, ii = 0;
    VertexInfo v_info;

    for (auto const &cell : cavity_intersecting_cells) {
        for (std::size_t i = 0; i < 3; i++) {
            if (intersecting_bool[ii][i]) {
                // This vtx is outside the included area
                continue;
            }
            // Get Info data structure
            v_info = cell->vertex(i)->info();
            // Get residue Index
            res_idx = v_info._resn;
            // Get position of the new index
            r_idx = std::lower_bound(
                wall_aa_idx.begin(), wall_aa_idx.end(), res_idx);

            if (r_idx == wall_aa_idx.end()) {
                // Store new residue (with largest index) info
                wall_aa_idx.push_back(res_idx);
                wall_aa_id.push_back(v_info._resi);
                continue;
            } else if (res_idx == *r_idx || res_idx == *(++r_idx)) {
                // This residue is already present
                continue;
            }

            r = (--r_idx) - wall_aa_idx.begin();
            // Store new residue info
            wall_aa_idx.insert(r_idx, res_idx);
            wall_aa_id.insert(wall_aa_id.begin() + r, v_info._resi);
        }
        ++ii;
    }

    return;
}

// Write wall amino acids and atoms
void write_wall_file(std::ofstream &pock_out_file,
    std::string const &pock_out_filename, const std::vector<int> &wall_aa_idx,
    const std::vector<std::string> &wall_aa_id,
    const std::vector<int> &wall_atom_idx, int const frame_cnt,
    std::string const &list_wall_separator) {

    std::size_t atom_length = wall_atom_idx.size();

    if (pock_out_file.is_open()) {
        // Header
        pock_out_file << pock_out_filename << "\tFrame: " << frame_cnt << '\n';

        pock_out_file << '\n' << "ATOM |\t";
        pock_out_file << wall_atom_idx[0];
        for (std::size_t i = 1; i < atom_length; ++i) {
            pock_out_file << list_wall_separator << wall_atom_idx[i];
        }
        pock_out_file << '\n' << "RESN |\t";
        pock_out_file << wall_aa_id[0];
        for (std::size_t i = 1; i < atom_length; ++i) {
            pock_out_file << list_wall_separator << wall_aa_id[i];
        }
        pock_out_file << '\n' << "RESI |\t";
        pock_out_file << wall_aa_idx[0];
        for (std::size_t i = 1; i < atom_length; ++i) {
            pock_out_file << list_wall_separator << wall_aa_idx[i];
        }

        pock_out_file << '\n';
        pock_out_file << "-------------" << '\n' << '\n';
    } else
        throw std::runtime_error(
            "Unable to open output file for Wall amino acids & atoms.");

    return;
}

// Write wall amino acids
void write_wall_file(std::ofstream &pock_out_file,
    std::string const &pock_out_filename, const std::vector<int> &wall_aa_idx,
    const std::vector<std::string> &wall_aa_id, int const frame_cnt,
    std::string const &list_wall_separator) {

    std::size_t resi_length = wall_aa_idx.size();

    if (pock_out_file.is_open()) {
        // Header
        pock_out_file << pock_out_filename << "\tFrame: " << frame_cnt << '\n';

        pock_out_file << '\n' << "RESN |\t";
        pock_out_file << wall_aa_id[0];
        for (std::size_t i = 1; i < resi_length; ++i) {
            pock_out_file << list_wall_separator << wall_aa_id[i];
        }

        pock_out_file << '\n';
        pock_out_file << "RESI |\t";

        pock_out_file << wall_aa_idx[0];
        for (std::size_t i = 1; i < resi_length; ++i) {
            pock_out_file << list_wall_separator << wall_aa_idx[i];
        }

        pock_out_file << '\n';
        pock_out_file << "-------------" << '\n' << '\n';
    } else
        throw std::runtime_error(
            "Unable to open output file for Wall amino acids & atoms.");

    return;
}
// Open output volume file, if requested.
void open_vol_file(std::string const &out_vol) {
    if (out_vol != "none") {
        out_vol_stream.open(out_vol);
        if (!(out_vol_stream)) {
            std::cerr << "Coulnd't open volume output file. Redirecting to "
                         "stdout"
                      << '\n';
        }
    }
    return;
}
// Final function to output volume. NA_Matrix (static) version.
void write_output_volume(NA_Matrix const &null_areas_vt_mt, double poly_vol) {

    int pock_cnt = 1;
    if (out_vol_stream.is_open()) {
        for (NA_Vector const &null_areas_vtor : null_areas_vt_mt) {
            double const volume = ANA::get_void_volume(null_areas_vtor);
            out_vol_stream << "Pocket " << pock_cnt << '\t' << volume + poly_vol
                           << '\n';
            ++pock_cnt;
        }
    } else {
        for (NA_Vector const &null_areas_vtor : null_areas_vt_mt) {
            double const volume = ANA::get_void_volume(null_areas_vtor);
            std::cout << "Pocket " << '\t' << volume + poly_vol << '\n';
            ++pock_cnt;
        }
    }
    return;
}

// Final function to output volume. NA_Vector (NDD) version.
void write_output_volume(
    NA_Vector const &null_areas_vtor, double const poly_vol) {

    double const volume = ANA::get_void_volume(null_areas_vtor);
    if (out_vol_stream.is_open()) {
        out_vol_stream << "Pocket " << '\t' << volume + poly_vol << '\n';
    } else {
        std::cout << "Pocket " << '\t' << volume + poly_vol << '\n';
    }
    return;
}

// Final function to output volume. NA_Vector (MD) version.
void write_output_volume(NA_Vector const &null_areas_vtor,
    double const poly_vol, int const frame_cnt) {

    double const volume = ANA::get_void_volume(null_areas_vtor);
    if (out_vol_stream.is_open()) {
        out_vol_stream << "Frame " << frame_cnt << '\t' << volume + poly_vol
                       << '\n';
    } else {
        std::cout << "Frame " << frame_cnt << '\t' << volume + poly_vol << '\n';
    }
    return;
}

} // namespace ANA
