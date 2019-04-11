#include <ANA/Cavity.hpp>

namespace ANA {

Cavity::Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts) {

    Delaunay const triangulation(molecule._data.begin(), molecule._data.end());

    Finite_cells_iterator fc_ite_end = triangulation.finite_cells_end();
    for (auto fc_ite = triangulation.finite_cells_begin(); fc_ite != fc_ite_end;
         ++fc_ite) {
        if (volume(fc_ite) > cell_opts._min_CV) {
            auto const v0 = fc_ite->vertex(0);
            auto const v1 = fc_ite->vertex(1);
            auto const v2 = fc_ite->vertex(2);
            auto const v3 = fc_ite->vertex(3);

            _all_cells.emplace_back(
                v0->point(), v1->point(), v2->point(), v3->point());
            _all_info.emplace_back(
                v0->info(), v1->info(), v2->info(), v3->info());
        }
    }

    return;
}

void Cavity::add_inner_cell(Tetrahedron const &cell, TetraInfo const &info) {
    _volume += volume(cell) - occupied_cell_vol(cell, info);
    _inner_cells.push_back(cell);
    _inner_info.push_back(info);
    return;
}

void Cavity::add_outer_cell(Tetrahedron const &cell, TetraInfo const &info) {
    _outer_cells.push_back(cell);
    _outer_info.push_back(info);
    return;
}

void Cavity::add_border_tetra(Point const &p0, Point const &ip1,
    Point const &ip2, Point const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - sphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    _tetra_border.emplace_back(p0, ip1, ip2, ip3);
    return;
}

void Cavity::add_border_penta_2_2(Point const &p0, Point const &p1,
    Point const &ip2, Point const &ip3, Point const &ip4, Point const &ip5,
    double const vdw0, double const vdw1) {

    double const vol1 = volume(p0, p1, ip2, ip3) -
        sphere_sector_vol(p0, p1, ip2, ip3, vdw0) -
        sphere_sector_vol(p1, p0, ip2, ip3, vdw1);

    double const vol2 =
        volume(p1, ip2, ip3, ip4) - sphere_sector_vol(p1, ip2, ip3, ip4, vdw1);

    double const vol3 =
        volume(p1, ip3, ip4, ip5) - sphere_sector_vol(p1, ip3, ip4, ip5, vdw1);
    _outer_volume += vol1 + vol2 + vol3;

    _penta_border_2_2.emplace_back(p0, p1, ip2, ip3, ip4, ip5);

    return;
}

void Cavity::add_border_penta_3_1(Point const &p0, Point const &p1,
    Point const &p2, Point const &ip3, Point const &ip4, Point const &ip5,
    double const vdw0, double const vdw1, double const vdw2) {

    double const vol1 = volume(p0, p1, p2, ip3) -
        sphere_sector_vol(p0, p1, p2, ip3, vdw0) -
        sphere_sector_vol(p1, p0, p2, ip3, vdw1) -
        sphere_sector_vol(p2, p1, p0, ip3, vdw2);

    double const vol2 = volume(p1, p2, ip3, ip4) -
        sphere_sector_vol(p1, p2, ip3, ip4, vdw1) -
        sphere_sector_vol(p2, p1, ip3, ip4, vdw2);

    double const vol3 =
        volume(p2, ip3, ip4, ip5) - sphere_sector_vol(p2, ip3, ip4, ip5, vdw2);

    _outer_volume += vol1 + vol2 + vol3;

    _penta_border_3_1.emplace_back(p0, p1, p2, ip3, ip4, ip5);

    return;
}

// Returns an updated Cavity displacing the input Cavity along the input
// vector scaled by the step_size. The vector must be alfa carbon mode.
Cavity::Cavity(Cavity const &hueco, std::vector<double> const &evector,
    double const step_size) {

    _all_cells.reserve(hueco._inner_cells.size() + hueco._outer_cells.size());
    _all_info.reserve(hueco._inner_info.size() + hueco._outer_info.size());

    move_cells(hueco._inner_cells, hueco._inner_info, evector, step_size);
    move_cells(hueco._outer_cells, hueco._outer_info, evector, step_size);

    return;
}

void Cavity::move_cells(std::vector<Tetrahedron> const &cells,
    std::vector<TetraInfo> const &info, std::vector<double> const &evector,
    double const step_size) {

    for (size_t c = 0; c < cells.size(); ++c) {
        TetraInfo ndd_info = info[c];
        // _resn is 1-indexed.
        int const resi_0_x = (ndd_info._resn[0] - 1) * 3;
        int const resi_1_x = (ndd_info._resn[1] - 1) * 3;
        int const resi_2_x = (ndd_info._resn[2] - 1) * 3;
        int const resi_3_x = (ndd_info._resn[3] - 1) * 3;

        Point const p0{cells[c][0] +
            step_size *
                Vector(evector[resi_0_x], evector[resi_0_x + 1],
                    evector[resi_0_x + 2])};
        Point const p1{cells[c][1] +
            step_size *
                Vector(evector[resi_1_x], evector[resi_1_x + 1],
                    evector[resi_1_x + 2])};
        Point const p2{cells[c][2] +
            step_size *
                Vector(evector[resi_2_x], evector[resi_2_x + 1],
                    evector[resi_2_x + 2])};
        Point const p3{cells[c][3] +
            step_size *
                Vector(evector[resi_3_x], evector[resi_3_x + 1],
                    evector[resi_3_x + 2])};

        _all_cells.emplace_back(p0, p1, p2, p3);
        _all_info.push_back(ndd_info);
    }

    return;
}

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(Tetrahedron const &cell, TetraInfo const &info) {

    double const vtx_0_vol =
        sphere_sector_vol(cell[0], cell[1], cell[2], cell[3], info._radius[0]);
    double const vtx_1_vol =
        sphere_sector_vol(cell[3], cell[0], cell[1], cell[2], info._radius[3]);
    double const vtx_2_vol =
        sphere_sector_vol(cell[2], cell[3], cell[0], cell[1], info._radius[2]);
    double const vtx_3_vol =
        sphere_sector_vol(cell[1], cell[2], cell[3], cell[0], info._radius[1]);

    return (vtx_0_vol + vtx_1_vol + vtx_2_vol + vtx_3_vol);
}

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream) {
    std::string temp;
    double coord = 0;

    in_stream >> temp;
    try {
        coord = std::stof(temp);
    } catch (std::invalid_argument const &ia) {
        // some character present
        throw std::runtime_error("Can't parse input. There may be "
                                 "non-numerical characters. Aborting.");
    } catch (...) {
        // some other exception.
        throw std::runtime_error("Can't parse input. Aborting.");
    }

    return coord;
}

void write_volume(Cavity const &hueco, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        fmt::print(out_file, "{: <7}{: >5}{:8.3f}\n",
            "Volume:", hueco._volume + hueco._outer_volume);
    } else {
        printf("Could not write NDD output to: %s.\n", filename.c_str());
    }

    return;
}

} // namespace ANA