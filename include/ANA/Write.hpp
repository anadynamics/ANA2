#ifndef ANA_WRITE_H
#define ANA_WRITE_H
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/PrimitivesUtils.hpp>
#include <ANA/Utils.hpp>

namespace ANA {
extern std::ofstream out_vol_stream;

// Draw a vector of polyhedrons
void draw_raw_polyhedrons(std::ofstream &pymol_script,
    const Poly_Vector &CH_vec, std::string const &model_template);

// Draw pockets in .PDB format. Static version
void draw_raw_PDB(NA_Vector const &list_of_pockets, const Poly_Vector &polys,
    std::string const &out_filename);

// Draw pockets in .PDB format. MD version, whole trajectory.
void draw_raw_PDB(const NDD_Matrix &list_of_pockets,
    const Poly_Matrix &list_of_polys, std::string &out_filename,
    int const max_atom_cnt);

// Draw pockets as dots in a PDB. MD version, whole trajectory.
void draw_grid_pdb(const NDD_Matrix &list_of_pockets,
    const Poly_Matrix &list_of_polys,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<std::vector<int>> &intersecting_total,
    int const sphere_count, int const precision, std::string &out_filename);

// Draw cells as dots in a PDB. Using "NDD_Vector" data structure.
int make_grid_pdb(NDD_Vector const &cells_to_draw,
    chemfiles::Topology &ana_void_top, chemfiles::Frame &ana_void_frame,
    int const sphere_count, int &res_cnt);

// Draw a vector of polyhedrons in a PDB.
int make_grid_pdb_polyhedrons(chemfiles::Topology &ana_void_top,
    chemfiles::Frame &ana_void_frame,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<int> &intersecting_total, const Poly_Vector &CH_vec,
    double const sphere_count, int atom_cnt, int &res_cnt);

// Draw pockets as dots in a PDB. Using CGAL Cell data
// structure. Static version.
void draw_grid_pdb(NA_Vector const &pocket,
    const std::vector<std::array<double, 3>> &in_vtces_radii,
    const std::vector<int> &intersecting_total, const Poly_Vector &polys,
    std::string &out_filename, int const sphere_count, int const precision);

// Draw cells as dots in a PDB. Using CGAL Cell data structure. Static version.
int make_grid_pdb(NA_Vector const &cells_to_draw,
    chemfiles::Topology &ana_void_top, chemfiles::Frame &ana_void_frame,
    int const sphere_count, int &res_cnt);

//////////
//
/////////

// Write .PDB header
void header_PDB(
    std::string const &out_pdb_filename, std::string const &in_pdb_filename);

// Construct a residue object with 4 atoms starting at index cell_cnt
chemfiles::Residue make_cell_residue(int const cell_cnt, int const atom_cnt);

// Construct a residue object with atoms from 'first_atom' to 'atom_cnt'
inline chemfiles::Residue make_polyhedron_residue(
    int const res_cnt, int const first_atom, int const atom_cnt);

// Construct a residue object for grid output
chemfiles::Residue make_grid_residue(
    int const atom_cnt_old, int const atom_cnt, int const resi);

// Connect a residue object with the atoms starting at index cell_cnt
void connect_cell_residue(chemfiles::Topology &topology, int const atom_index);

// Connect a residue object with the atoms starting at index atom_index
inline void connect_facet(chemfiles::Topology &topology, int const atom_index);

//////////
//
/////////

// Draw a whole convex hull contained in a polyhedron. Static version.
void draw_CH(const Polyhedron &CH, std::string &out_filename);

// Draw a whole convex hull contained in vector of triangles. Static version.
void draw_CH(Triang_Vector const &CH_triang, std::string &out_filename);

// Draw a whole convex hull contained in vector of triangles. MD version.
void draw_CH(Triang_Vector const &CH_triang, chemfiles::Trajectory &out_traj);

// Draw a whole triangulation
void draw_triangulation(Delaunay const &T, std::string &script_filename);

// Write wall amino acids and atoms.
void wall_atom_output(std::ofstream &wall_out, NA_Vector const &in_cells,
    NA_Vector const &in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    bool const requested_CH, int const precision, int const pock_cnt,
    int const frame_cnt, std::string const &list_wall_separator);

// Write wall amino acids.
void wall_aa_output(std::ofstream &wall_out, NA_Vector const &in_cells,
    NA_Vector const &in_intersecting_cells,
    const std::vector<std::array<bool, 4>> intersecting_bool,
    bool const requested_CH, int const precision, int const pock_cnt,
    int const frame_cnt, std::string const &list_wall_separator);

// Get names and indices of participating atoms and amino acids.
void get_info_cell(NA_Vector const &null_areas_vtor,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id,
    std::vector<int> &wall_atom_idx);

// Get names and indices of participating atoms and amino acids for intersecting
// cells.
void get_info_cell(NA_Vector const &cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id,
    std::vector<int> &wall_atom_idx);

// Get names and indices of participating amino acids.
void get_info_cell(NA_Vector const &null_areas_vtor,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id);

// Get names and indices of participating amino acids from intersecting cells.
void get_info_cell(NA_Vector const &cavity_intersecting_cells,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &wall_aa_idx, std::vector<std::string> &wall_aa_id);

// Write wall amino acids and atoms.
void write_wall_file(std::ofstream &pock_out_file,
    std::string const &pock_out_filename, const std::vector<int> &wall_aa_idx,
    const std::vector<std::string> &wall_aa_id,
    const std::vector<int> &wall_atom_idx, int const frame_cnt,
    std::string const &list_wall_separator);

// Write wall amino acids.
void write_wall_file(std::ofstream &pock_out_file,
    std::string const &pock_out_filename, const std::vector<int> &wall_aa_idx,
    const std::vector<std::string> &wall_aa_id, int const frame_cnt,
    std::string const &list_wall_separator);

// Final function to output volume. NA_Matrix (static) version.
void write_output_volume(
    NA_Matrix const &null_areas_vt_mt, double const poly_vol);

// Final function to output volume. NA_Vector (NDD) version.
void write_output_volume(
    NA_Vector const &null_areas_vt_mt, double const poly_vol);

// Final function to output volume. NA_Vector (MD) version.
void write_output_volume(NA_Vector const &null_areas_vtor,
    double const poly_vol, int const frame_cnt);

// Open output volume file, if requested.
void open_vol_file(std::string const &out_vol);

// Final function to output volume. NA_Vector (MD) version.
void write_output_volume(NA_Vector const &null_areas_vtor,
    double const poly_vol, int const frame_cnt);
}
namespace ANA {
namespace NDD {
    // Write file with volumes of each input PDB
    inline void ndd_write_out_file(const std::vector<double> &output_volumes,
        std::string const &out_file) {
        std::ofstream output(out_file);

        int i = 1;

        if (output.is_open()) {
            output << "Frame\tVolume" << '\n';

            for (double volume : output_volumes) {
                output << i << "\t" << volume << '\n';

                ++i;
            }
        } else
            throw std::runtime_error("Unable to open output file for NDD");

        return;
    }

} // namespace NDD
} // namespace ANA
#endif // _H
