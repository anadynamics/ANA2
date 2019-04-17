#include <ANA/MD.hpp>
#include <ANA/NDD.hpp>
#include <ANA/Options.hpp>
#include <ANA/ProgramOptions.hpp>
#include <ANA/Static.hpp>

int main(int argc, char *argv[]) {
    ANA::IncludedAreaOptions IA_opts;
    ANA::NDDOptions NDD_opts;
    ANA::CellFilteringOptions cell_opts;
    ANA::InOutOptions io_opts;

    std::string AA_indices_proto, exclude_ca_for_ASA_indices_proto,
        clusters_method = "facets", ASA_method = "dot_pdt", only_side_ASA,
        list_wall, list_wall_separator, tool_check_CH, tool_pdb_to_ch,
        tool_pdb_norm, tool_aa_to_ca;

    bool triangulate_only_included_aas = false, atom_only = true;

    int nbr_of_vertices_to_include = 2, min_cells_cluster = 2, estado = 1,
        md_start, md_step, md_end, precision, sphere_count;
    std::vector<int> include_CH_aa, hetatm_atoms;
    double max_probe, max_probe_length;

    estado = ANA::get_parameters(argc, argv, io_opts, IA_opts, AA_indices_proto,
        triangulate_only_included_aas, atom_only, precision, min_cells_cluster,
        nbr_of_vertices_to_include, md_start, md_step, md_end, cell_opts,
        max_probe, max_probe_length, sphere_count, list_wall,
        list_wall_separator, clusters_method, only_side_ASA, ASA_method,
        exclude_ca_for_ASA_indices_proto, NDD_opts, tool_check_CH,
        tool_pdb_to_ch, tool_pdb_norm, tool_aa_to_ca);
    if (estado == 1) {
        // Some error was found. Terminating execution.
        return 0;
    }

    //////////////////////////////////////////////////////////////////////////////
    // tools
    //////////////////////////////////////////////////////////////////////////////

    // Tool for CH included area visual check.
    if (tool_check_CH != "none") {
        Triang_Vector CH_triangs;
        std::vector<CPoint> CAs_Points;
        ANA_molecule molecule_points;
        std::vector<int> AA_indices, CA_indices, include_CH_atoms;
        CPoint cm;

        // Read input file
        bool const requested_CH = ANA::read_static(io_opts._in_filename,
            triangulate_only_included_aas, atom_only, AA_indices_proto,
            exclude_ca_for_ASA_indices_proto, IA_opts._resn_proto,
            IA_opts._atom_proto, IA_opts._sphere_proto, IA_opts._cylinder_proto,
            IA_opts._prism_proto, IA_opts._filename, molecule_points, cm,
            AA_indices, CA_indices, CAs_Points, include_CH_atoms, CH_triangs,
            hetatm_atoms);

        if (!requested_CH) {
            std::cerr << "No valid input for triangulation of included area. "
                         "Exiting now."
                      << '\n';
            return 0;
        }

        if (io_opts._in_md_filename == "none") {

            ANA::draw_CH(CH_triangs, tool_check_CH);

        } else {

            chemfiles::Trajectory in_traj(io_opts._in_md_filename);
            int frame_cnt = 1;
            // Set end.
            if (md_end == 0) {
                md_end = in_traj.nsteps();
            }
            // Set output file.
            tool_check_CH.append(".pdb");
            auto out_traj = chemfiles::Trajectory(tool_check_CH, 'w');

            while (!in_traj.done()) {
                // Set next step.
                int const current_step =
                    (frame_cnt - 1) * md_step + (md_start - 1);
                if (current_step >= md_end) {
                    // Done.
                    break;
                }
                // Read next frame.
                chemfiles::Frame in_frame = in_traj.read_step(current_step);
                // Update CH
                ANA::read_MD(in_frame, requested_CH, IA_opts._sphere_proto,
                    IA_opts._cylinder_proto, IA_opts._prism_proto, hetatm_atoms,
                    include_CH_atoms, IA_opts._filename, CH_triangs, ASA_method,
                    CA_indices, CAs_Points, molecule_points);
                // New CH model.
                ANA::draw_CH(CH_triangs, out_traj);
                ++frame_cnt;
            }
        }
        return 0;
    }

    // Tool to get convex hull of input PDB and writing its vertices to a file.
    if (tool_pdb_to_ch != "none") {
        ANA::tool_PDB_to_CH(io_opts._in_filename, tool_pdb_to_ch);
        return 0;
    }

    // Tool for normalizing PDB, by renumbering its atoms and residues.
    if (tool_pdb_norm != "none") {
        tool_pdb_norm.append(".pdb");
        ANA::tool_PDB_norm(io_opts._in_filename, tool_pdb_norm);
        return 0;
    }

    // Tool for getting the indices of the Calpha atoms from the
    // included_area_residues config.
    if (tool_aa_to_ca != "none") {

        Triang_Vector CH_triangs;
        std::vector<CPoint> CAs_Points;
        ANA_molecule molecule_points;
        std::vector<int> AA_indices, CA_indices, include_CH_atoms;
        CPoint cm;

        // Read input file
        ANA::read_static(io_opts._in_filename, triangulate_only_included_aas,
            atom_only, AA_indices_proto, exclude_ca_for_ASA_indices_proto,
            IA_opts._resn_proto, IA_opts._atom_proto, IA_opts._sphere_proto,
            IA_opts._cylinder_proto, IA_opts._prism_proto, IA_opts._filename,
            molecule_points, cm, AA_indices, CA_indices, CAs_Points,
            include_CH_atoms, CH_triangs, hetatm_atoms);
        std::sort(include_CH_atoms.begin(), include_CH_atoms.end());

        std::cout << "\t\t/// Calpha indices ///" << '\n';
        std::cout << "include_area_residues = ";
        for (auto const &each : include_CH_atoms) {
            std::cout << each + 1 << " ";
        }
        std::cout << '\n';

        return 0;
    }

    //////////////////////////////////////////////////////////////////////////////
    // end tools
    //////////////////////////////////////////////////////////////////////////////
    // Get output volume file ready, if requested.
    ANA::open_vol_file(io_opts._out_vol_filename);

    if (io_opts._in_md_filename != "none") {

        ANA::MD_ANA(io_opts._in_filename, io_opts._in_md_filename,
            AA_indices_proto, ASA_method, only_side_ASA,
            exclude_ca_for_ASA_indices_proto, list_wall, list_wall_separator,
            IA_opts._resn_proto, IA_opts._atom_proto, IA_opts._sphere_proto,
            IA_opts._cylinder_proto, IA_opts._prism_proto, IA_opts._filename,
            io_opts._out_pdb_filename, io_opts._out_type,
            triangulate_only_included_aas, atom_only, cell_opts, max_probe,
            max_probe_length, sphere_count, nbr_of_vertices_to_include,
            precision, md_start, md_step, md_end);
    } else if (NDD_opts._modes_ndd_filename != "none") {
        IA_opts.has_info = true;
        ANA::NDD_ANA(io_opts, IA_opts, NDD_opts, cell_opts, atom_only);
    } else {

        ANA::static_ANA(io_opts._in_filename, AA_indices_proto, ASA_method,
            only_side_ASA, exclude_ca_for_ASA_indices_proto, list_wall,
            list_wall_separator, clusters_method, IA_opts._resn_proto,
            IA_opts._atom_proto, IA_opts._sphere_proto, IA_opts._cylinder_proto,
            IA_opts._prism_proto, IA_opts._filename, io_opts._out_pdb_filename,
            io_opts._out_type, triangulate_only_included_aas, atom_only,
            cell_opts, max_probe, max_probe_length, sphere_count,
            nbr_of_vertices_to_include, min_cells_cluster, precision);
    }

    return 0;
}
