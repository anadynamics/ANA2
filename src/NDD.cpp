#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(in_filename, atom_only);

    ANA::ConvexHull const CH = create_convex_hull(protein, IA_opts);

    ANA::Cavity hueco = Cavity(protein, cell_opts);

    ANA::carve_CH_into_cavity(hueco, CH);

    ANA::NDD::ndd(hueco, CH, NDD_opts);

    ANA::write_PDB(CH, "hull.pdb");
    ANA::write_PDB(hueco, "sal.pdb");

    printf("Volumen:  %f\n", hueco._volume + hueco._outer_volume);

    return 0;
}

} // namespace ANA
