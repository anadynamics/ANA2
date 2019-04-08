#ifndef ANA_NDD_UTILS_H
#define ANA_NDD_UTILS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/ConvexHullFunctions.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Modes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/PrimitivesUtils.hpp>
#include <ANA/Read.hpp>
#include <ANA/Utils.hpp>
// DELETE after testing.
#include <ANA/PDB.hpp>

using std::size_t;

namespace ANA::NDD {

// NDD Specific function for PDB input
void ndd_read_PDB_get_cells(std::string const &filename,
    NDD_IVector const &in_cells_indices, NDD_Vector &output_cells);

// NDD Specific function for PDB input. Hi precision method
void ndd_read_PDB_get_cells(std::string const &filename,
    NDD_IVector const &in_void_cells_indices,
    const std::vector<int> &include_CH_atoms, NDD_Vector &output_cells,
    Triang_Vector &CH_triangs);

// On-site NDD.
void ndd(Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts);

void write_vector(std::vector<double> vec, std::string const &filename);

// Perform Non Delaunay Dynamics.
void ndd_nondelaunay_dynamics_old(NA_Vector const &cavity_void_cells,
    std::string const &pdb_list, bool const precision,
    const std::vector<int> include_CH_atoms, std::string const &out_file);

// Get the indices of the atoms involved in the given cells
NDD_IVector get_vertices(NA_Vector const &cavity_void_cells);

// Get the volume occupied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(CPoint const &p_0, CPoint const &p_1,
    CPoint const &p_2, CPoint const &p_3, double const radius);

// Calc volume of the input cells. Reedited for array container.
double ndd_get_void_volume(NDD_Vector const &cavity_void_cells);

// Substract the volume filled with the 4 atoms from the total volume of
// the corresponding cell. Reedited for array container.
inline double ndd_refine_cell_volume(
    double const entire_cell_vol, const NDD_Element &cell) {
    CPoint const p_0 = cell[0].first;

    CPoint const p_1 = cell[1].first;

    CPoint const p_2 = cell[2].first;

    CPoint const p_3 = cell[3].first;

    double const vtx_0_sphere_sector =
        sphere_sector_vol(p_0, p_1, p_2, p_3, cell[0].second);

    double const vtx_1_sphere_sector =
        sphere_sector_vol(p_3, p_0, p_1, p_2, cell[3].second);

    double const vtx_2_sphere_sector =
        sphere_sector_vol(p_2, p_3, p_0, p_1, cell[2].second);

    double const vtx_3_sphere_sector =
        sphere_sector_vol(p_1, p_2, p_3, p_0, cell[1].second);

    return (entire_cell_vol - vtx_0_sphere_sector - vtx_1_sphere_sector -
        vtx_2_sphere_sector - vtx_3_sphere_sector);
}

// Discard cells without a vertex inside the specified convex hull. Hi
// precision, NDD version
void ndd_discard_CH_0(NDD_Vector const &in_cells,
    Triang_Vector const &CH_triangs, NDD_Vector &out_cells,
    NDD_Vector &out_intersecting_cells,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &intersecting_total);

// Discard parts of cells outside the specified triangulation using
// intersecitons
double ndd_discard_CH_1(NDD_Vector const &in_intersecting_coords,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<int> &intersecting_total, Poly_Vector &border_poly);

} // namespace NDD:: ANA

#endif // _H
