#include <ANA/Static.hpp>

namespace ANA {

int static_ANA(InOutOptions const &io_opts, std::string &AA_indices_proto,
    std::string const &ASA_method, std::string const &only_side_ASA,
    std::string &exclude_ca_for_ASA_indices_proto, std::string const &list_wall,
    std::string const &list_wall_separator, std::string const &clusters_method,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    std::string &out_filename, std::string const &out_type,
    bool const triangulate_only_included_aas, bool const atom_only,
    CellFilteringOptions const cell_opts, double const max_probe,
    double const max_probe_length, int const sphere_count,
    int const nbr_of_vertices_to_include, int const clusters_min_size,
    int const precision) {

    // atom_cnt_poly is for MD only.
    int atom_cnt_poly = 0;
    double poly_vol = 0;
    std::vector<CPoint> CAs_Points;
    ANA_molecule molecule_points;
    std::vector<int> AA_indices, CA_indices, include_CH_atoms;
    CPoint cm;
    Triang_Vector CH_triangs;
    NA_Vector cavity_cells, cavity_included_cells, cavity_void_cells,
        cavity_intersecting_cells;
    NA_Matrix null_areas_mtx;
    Poly_Vector border_poly;
    std::vector<std::array<bool, 4>> intersecting_bool;
    std::vector<int> intersecting_total;
    std::vector<int> wall_aa_idx, wall_atom_idx;
    std::vector<int> hetatm_atoms;
    std::vector<std::string> wall_aa_id;
    std::vector<std::array<double, 3>> in_vtces_radii;

    // Read input file
    bool const requested_CH = ANA::read_static(io_opts._in_filename,
        triangulate_only_included_aas, atom_only, AA_indices_proto,
        exclude_ca_for_ASA_indices_proto, include_CH_aa_proto,
        include_CH_atom_proto, sphere_proto, cylinder_proto, prism_proto,
        include_CH_filename, molecule_points, cm, AA_indices, CA_indices,
        CAs_Points, include_CH_atoms, CH_triangs, hetatm_atoms);

    Delaunay T = ANA::triangulate(molecule_points);

    ANA::get_all_voids(T, cavity_cells, cell_opts);

    if (!requested_CH) {
        // Discover pockets
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

        if (clusters_method == "facets") {
            ANA::cluster_cells_cgal(
                cavity_void_cells, null_areas_mtx, clusters_min_size);
        } else if (clusters_method == "boxes") {
            ANA::cluster_cells_boxes(cavity_void_cells, null_areas_mtx);
        } else {
            null_areas_mtx.push_back(cavity_void_cells);
        }

        if (out_filename != "none") {
            if (out_type == "raw_pdb") {
                // Get ready to write the pockets in different .pdbs
                int pock_cnt = 1;
                for (NA_Vector const &null_areas_vtor : null_areas_mtx) {

                    std::string const pock_out_filename =
                        get_output_pocket_filename(out_filename, pock_cnt);

                    ANA::draw_raw_PDB(
                        null_areas_vtor, border_poly, pock_out_filename);
                    ++pock_cnt;
                }
            } else if (out_type == "grid_pdb") {
                // Get ready to write the pockets in different .pdbs
                int pock_cnt = 1;
                for (NA_Vector const &null_areas_vtor : null_areas_mtx) {

                    std::string const pock_out_filename =
                        get_output_pocket_filename(out_filename, pock_cnt);

                    ANA::draw_grid_pdb(null_areas_vtor, in_vtces_radii,
                        intersecting_total, border_poly, pock_out_filename,
                        sphere_count, precision);
                    ++pock_cnt;
                }
            }
        }

    } else {
        // Included area was specified.
        if (precision == 1) {
            discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells,
                cavity_intersecting_cells, intersecting_bool,
                intersecting_total);
            poly_vol = discard_CH_1(cavity_intersecting_cells, CH_triangs,
                intersecting_bool, intersecting_total, border_poly,
                in_vtces_radii, atom_cnt_poly);
        } else { // assume precision = 0
            discard_CH_0(cavity_cells, CH_triangs, cavity_void_cells);
        }
        // null_areas_mtx.push_back(cavity_void_cells);

        if (out_filename != "none") {
            if (out_type == "raw_pdb") {
                std::string const pock_out_filename =
                    get_output_pocket_filename(out_filename);

                ANA::draw_raw_PDB(
                    cavity_void_cells, border_poly, pock_out_filename);

            } else if (out_type == "grid_pdb") {

                std::string const pock_out_filename =
                    get_output_pocket_filename(out_filename);

                ANA::draw_grid_pdb(cavity_void_cells, in_vtces_radii,
                    intersecting_total, border_poly, pock_out_filename,
                    sphere_count, precision);
            }
        }
    }

    if (list_wall == "atom") {
        int pock_cnt = 1;
        std::ofstream wall_out(io_opts._out_wall_filename);
        for (NA_Vector const &null_areas_vtor : null_areas_mtx) {
            ANA::wall_atom_output(wall_out, null_areas_vtor,
                cavity_intersecting_cells, intersecting_bool, requested_CH,
                precision, pock_cnt, 1, list_wall_separator);
            ++pock_cnt;
        }
    } else if (list_wall == "residue") {
        int pock_cnt = 1;
        std::ofstream wall_out(io_opts._out_wall_filename);
        for (NA_Vector const &null_areas_vtor : null_areas_mtx) {
            ANA::wall_aa_output(wall_out, null_areas_vtor,
                cavity_intersecting_cells, intersecting_bool, requested_CH,
                precision, pock_cnt, 1, list_wall_separator);
            ++pock_cnt;
        }
    }

    ANA::write_output_volume(null_areas_mtx, poly_vol);

    return 0;
}

} // namespace ANA
