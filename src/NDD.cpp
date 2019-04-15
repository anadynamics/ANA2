#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(InOutOptions const &io_opts, IncludedAreaOptions const &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(io_opts._in_filename, atom_only);

    ANA::ConvexHull const CH = create_convex_hull(protein, IA_opts);

    ANA::Cavity hueco = Cavity(protein, cell_opts);

    ANA::carve_CH_into_cavity(hueco, CH);

    // ANA::write_PDB(hueco, "cav.pdb");
    // ANA::write_PDB(CH, "ch.pdb");

    ANA::NDD::ndd(hueco, CH, NDD_opts);

    if (io_opts._out_vol_filename == "none") {
        printf("Volumen:  %f\n", hueco._volume + hueco._outer_volume);
    } else {
        write_volume(hueco, io_opts._out_vol_filename);
    }

    return 0;
}

} // namespace ANA
