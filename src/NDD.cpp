#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(InOutOptions const &io_opts, IncludedAreaOptions const &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(io_opts._in_filename, atom_only);

    ANA::ConvexHull const CH = create_convex_hull(protein, IA_opts);

    ANA::Cavity hueco = Cavity(protein, cell_opts);

    ANA::carve_CH_into_cavity(hueco, CH);

    ANA::NDD::ndd(hueco, CH, NDD_opts, io_opts._in_filename);

    return 0;
}

} // namespace ANA
