#include <ANA/MD.hpp>

namespace ANA {

int MD_ANA(std::string const &in_filename, std::string const &in_md_filename,
    std::string &AA_indices_proto, std::string const &ASA_method,
    std::string const &only_side_ASA,
    std::string &exclude_ca_for_ASA_indices_proto, std::string const &list_wall,
    std::string const &list_wall_separator, std::string &include_CH_aa_proto,
    std::string &include_CH_atom_proto, std::string &sphere_proto,
    std::string &cylinder_proto, std::string &prism_proto,
    std::string const &include_CH_filename, std::string &out_filename,
    std::string const &out_type, bool const triangulate_only_included_aas,
    bool const atom_only, CellFilteringOptions const cell_opts,
    double const max_probe, double const max_probe_length,
    int const &sphere_count, int const nbr_of_vertices_to_include,
    int const precision, int const &md_start, int const &md_step, int &md_end) {

    // Read trajectory.
    chemfiles::Trajectory in_traj(in_md_filename);
    // Set end.
    if (md_end == 0) {
        md_end = in_traj.nsteps();
    }

    // Get ready.
    double poly_vol = 0;
    int frame_cnt = 1;
    int max_atom_cnt = 0;
    auto const nbr_of_frames =
        ceil((md_end - md_start + 1) / static_cast<double>(md_step));

    std::vector<int> atom_cnt_md(nbr_of_frames);
    std::vector<int> CA_indices, AA_indices, include_CH_atoms;
    std::vector<int> hetatm_atoms;
    CPoint cm;
    Polyhedron CH;
    Triang_Vector CH_triangs;
    std::vector<CPoint> CAs_Points;
    ANA_molecule molecule_points;
    std::vector<std::array<double, 3>> in_vtces_radii;
    NDD_Matrix voids_md(nbr_of_frames);
    Poly_Matrix polys_md(nbr_of_frames);
    std::vector<std::vector<int>> list_intersecting_total(nbr_of_frames);

    // Handle input trajectory file format.
    std::string in_md_format =
        in_md_filename.substr(in_md_filename.length() - 2, 2);
    bool const in_md_nc = (in_md_format == "nc");
    if (md_step != 1 && !in_md_nc) {
        std::cerr << "Warning: netcdf is the only format that supports an "
                     "\"md_step\" value other than 1. Reading every step."
                  << '\n';
    }

    // Remove ".pdb" and create file.
    std::string filename = in_filename.substr(0, in_filename.size() - 4);
    filename.insert(0, "wall_");
    std::ofstream wall_out;
    if (list_wall == "atom" || list_wall == "residue") {
        wall_out.open(filename);
    }
    // Get topology
    bool const requested_CH = ANA::read_static(in_filename,
        triangulate_only_included_aas, atom_only, AA_indices_proto,
        exclude_ca_for_ASA_indices_proto, include_CH_aa_proto,
        include_CH_atom_proto, sphere_proto, cylinder_proto, prism_proto,
        include_CH_filename, molecule_points, cm, AA_indices, CA_indices,
        CAs_Points, include_CH_atoms, CH_triangs, hetatm_atoms);

    while (!in_traj.done()) {
        int atom_cnt_poly = 0;
        NA_Vector cavity_cells;
        NA_Vector cavity_included_cells;
        NA_Vector cavity_void_cells;
        NA_Vector cavity_intersecting_cells;
        Poly_Vector border_poly;
        std::vector<std::array<bool, 4>> intersecting_bool;
        std::vector<int> intersecting_total;

        // Set next step.
        int const current_step = (frame_cnt - 1) * md_step + (md_start - 1);
        if (current_step >= md_end) {
            // Done.
            break;
        }

        // Read next frame.
        chemfiles::Frame in_frame;
        if (in_md_nc) {
            in_frame = in_traj.read_step(current_step);
        } else {
            in_frame = in_traj.read();
        }

        // Get frame coordinates
        ANA::read_MD(in_frame, requested_CH, sphere_proto, cylinder_proto,
            prism_proto, hetatm_atoms, include_CH_atoms, include_CH_filename,
            CH_triangs, ASA_method, CA_indices, CAs_Points, molecule_points);

        // Triangulate frame
        Delaunay T = ANA::triangulate(molecule_points);

        ANA::get_all_voids(T, cavity_cells, cell_opts);

        if (!requested_CH) {

            if (triangulate_only_included_aas == false && AA_indices[0] != 0) {
                ANA::keep_included_aa_cells(cavity_cells, AA_indices,
                    nbr_of_vertices_to_include, cavity_included_cells);
            } else {
                cavity_included_cells = cavity_cells;
            }

            if (ASA_method == "none") {
                cavity_void_cells = cavity_included_cells;
            } else if (ASA_method == "cm") {
                ANA::discard_ASA_dot_pdt_cm(cm, CAs_Points, max_probe,
                    max_probe_length, only_side_ASA, cavity_included_cells,
                    cavity_void_cells);
            } else if (ASA_method == "backbone") {
                ANA::discard_ASA_CACH(CAs_Points, only_side_ASA,
                    cavity_included_cells, cavity_void_cells);
            } else if (ASA_method == "axes") {
                ANA::discard_ASA_dot_pdt_axes(CAs_Points, max_probe,
                    max_probe_length, only_side_ASA, cavity_included_cells,
                    cavity_void_cells);
            }

        } else {
            if (precision == 1) {
                discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells,
                    cavity_intersecting_cells, intersecting_bool,
                    intersecting_total);
                poly_vol = discard_CH_1(cavity_intersecting_cells, CH_triangs,
                    intersecting_bool, intersecting_total, border_poly,
                    in_vtces_radii, atom_cnt_poly);
                // Store this frame's polyhedrons and the number of
                // intersections.
                polys_md[frame_cnt - 1] = std::move(border_poly);
                list_intersecting_total[frame_cnt - 1] =
                    std::move(intersecting_total);
            } else { // assume precision = 0
                discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells);
            }
        }

        if (list_wall == "atom") {
            ANA::wall_atom_output(wall_out, cavity_void_cells,
                cavity_intersecting_cells, intersecting_bool, requested_CH,
                precision, 1, frame_cnt, list_wall_separator);

        } else if (list_wall == "residue") {
            // Write output.
            ANA::wall_aa_output(wall_out, cavity_void_cells,
                cavity_intersecting_cells, intersecting_bool, requested_CH,
                precision, 1, frame_cnt, list_wall_separator);
        }

        ANA::write_output_volume(cavity_void_cells, poly_vol, frame_cnt);

        // Get the number of vertices (atoms) to write, look for the highest nbr
        // of "voids_md" them, and store it. Then, store voids for later output.
        int atom_cnt_step = cavity_void_cells.size() * 4 + atom_cnt_poly;
        atom_cnt_md[frame_cnt - 1] = atom_cnt_step;
        if (atom_cnt_step > max_atom_cnt) {
            max_atom_cnt = atom_cnt_step;
        }

        // Store this frame's cells.
        NDD_Vector temp_md_vec;
        na_vector_into_ndd_vector(cavity_void_cells, temp_md_vec);
        voids_md[frame_cnt - 1] = std::move(temp_md_vec);
        ++frame_cnt;
    }

    // Now, if requested, write the output.
    std::cerr << "Writing output..." << '\n';
    if (out_filename != "none") {
        if (out_type == "raw_pdb") {
            ANA::draw_raw_PDB(voids_md, polys_md, out_filename, max_atom_cnt);
        } else if (out_type == "grid_pdb") {
            ANA::draw_grid_pdb(voids_md, polys_md, in_vtces_radii,
                list_intersecting_total, sphere_count, precision, out_filename);
        }
    }

    std::cerr << "Done" << '\n';

    return 0;
}

} // namespace ANA
