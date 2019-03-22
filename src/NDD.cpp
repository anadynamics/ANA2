#include <ANA/NDD.hpp>

namespace ANA {

int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only) {

    Molecule const protein = ANA::Molecule(in_filename, atom_only);

    // IA_opts._resn_proto = "1 12 21 51 68 97 111 119 71 130 133 128 16 24 84";
    // IA_opts._opt = IncludedAreaOptions::residue;

    ANA::ConvexHull const CH = create_convex_hull(protein, IA_opts);

    ANA::Cavity hueco = Cavity(protein, cell_opts);

    ANA::carve_CCH_into_cavity(hueco, CH);

    ANA::NDD::ndd(protein, hueco, CH, IA_opts, NDD_opts, in_filename);

    // ANA::write_output_volume(cavity_void_cells, poly_vol);

    Point p0{1., 0., 0.}, p1{0., 1., 0.}, p2{0., 0., 1.}, q{1.1, 1.3, 3.5},
        r{0.2, .3, -.5};
    TTriangle t(p0, p1, p2);

    Vector u = p1 - p0, v = p2 - p0;
    Vector d = q - p0, s = r - q;
    Vector N = normalize(cross_product(u, v));

    double const det = determinant(s, u, v);
    double const det1 = determinant(d, u, v);
    double const det2 = determinant(s, d, v);
    double const det3 = determinant(s, u, d);

    double const l2 = det2 / det;
    double const l3 = det3 / det;
    double const l1 = 1 - l2 - l3;

    Point i(l1, l2, l3);
    std::cout << "l1: " << l1 << " l2: " << l2 << " l3:" << l3 << '\n';

    FILE *out_file = std::fopen("sal.pdb", "w");
    if (out_file) {

        std::pair<int, int> idx_resid{1, 1};
        idx_resid = draw_lines(t, out_file, idx_resid);
        idx_resid.first = draw(q, out_file, idx_resid, "AAA");
        idx_resid.first = draw(r, out_file, idx_resid, "AAA");
        idx_resid.first = draw(i, out_file, idx_resid, "AAA");
        connect_triangle(out_file, 1, idx_resid.second);

    } else {
        printf("Could not open %s.\n", "sal.pdb");
    }
    std::fclose(out_file);

    // ANA::write_PDB(hueco, "sal.pdb");
    // ANA::write_PDB(CH, "hull.pdb");
    printf("Volumen:  %f\n", hueco._volume + hueco._outer_volume);

    return 0;
}

} // namespace ANA
