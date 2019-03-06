#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(in_filename, atom_only);

    // IA_opts._resn_proto = "1 12 21 51 68 97 111 119 71 130 133 128 16 24 84";
    // IA_opts._opt = IncludedAreaOptions::residue;

    ANA::CConvexHull CH = create_cconvex_hull(protein, IA_opts);

    ANA::CCavity hueco(protein, cell_opts);

    ANA::carve_CCH_into_cavity(hueco, CH);

    // ANA::NDD::ndd(hueco, CH, NDD_opts);

    // ANA::write_output_volume(cavity_void_cells, poly_vol);

    ANA::write_PDB(hueco, "sal.pdb");
    ANA::write_PDB(CH, "hull.pdb");
    printf("Volumen:  %f\n", hueco._volume + hueco._outer_volume);

    ANA::NDD::Modes modo(NDD_opts._modes_ndd_filename);
    modo.calpha_to_full_atom(in_filename);

    return 0;
}

} // namespace ANA
