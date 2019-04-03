#include <ANA/Cavity.hpp>

namespace ANA {

Cavity::Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts) {

    _triangulation.insert(molecule._data.begin(), molecule._data.end());

    Finite_cells_iterator fc_ite_end = _triangulation.finite_cells_end();
    for (auto fc_ite = _triangulation.finite_cells_begin();
         fc_ite != fc_ite_end; ++fc_ite) {
        if (volume(fc_ite) > cell_opts._min_CV) {
            _all_cells.push_back(fc_ite);
        }
    }

    return;
}

void Cavity::add_inner_cell(Finite_cells_iterator const &cell) {
    _volume += volume(cell) - occupied_cell_vol(cell);
    _inner_cells.push_back(cell);
    return;
}

void Cavity::add_border_tetra(CPoint const &p0, CPoint const &ip1,
    CPoint const &ip2, CPoint const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - sphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    _tetra_border.emplace_back(p0, ip1, ip2, ip3);
    return;
}

void Cavity::add_border_penta(CPoint const &p0, CPoint const &p1,
    CPoint const &ip2, CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
    double const vdw0, double const vdw1) {

    double const vol1 = volume(p0, p1, ip2, ip3) -
        sphere_sector_vol(p0, p1, ip2, ip3, vdw0) -
        sphere_sector_vol(p1, p0, ip2, ip3, vdw1);

    double const vol2 =
        volume(p1, ip2, ip3, ip4) - sphere_sector_vol(p1, ip2, ip3, ip4, vdw1);

    double const vol3 =
        volume(p1, ip3, ip4, ip5) - sphere_sector_vol(p1, ip3, ip4, ip5, vdw1);
    _outer_volume += vol1 + vol2 + vol3;

    _penta_border.emplace_back(p0, p1, ip2, ip3, ip4, ip5);

    return;
}

void Cavity::add_border_penta(CPoint const &p0, CPoint const &p1,
    CPoint const &p2, CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
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

    _penta_border.emplace_back(p0, p1, p2, ip3, ip4, ip5);

    return;
}

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(Finite_cells_iterator const cell_iterator) {

    double const rdW_0 = double(cell_iterator->vertex(0)->info()._radius);
    double const rdW_1 = double(cell_iterator->vertex(1)->info()._radius);
    double const rdW_2 = double(cell_iterator->vertex(2)->info()._radius);
    double const rdW_3 = double(cell_iterator->vertex(3)->info()._radius);
    CPoint const p_0 = cell_iterator->vertex(0)->point();
    CPoint const p_1 = cell_iterator->vertex(1)->point();
    CPoint const p_2 = cell_iterator->vertex(2)->point();
    CPoint const p_3 = cell_iterator->vertex(3)->point();

    double const vtx_0_vol = sphere_sector_vol(p_0, p_1, p_2, p_3, rdW_0);
    double const vtx_1_vol = sphere_sector_vol(p_3, p_0, p_1, p_2, rdW_3);
    double const vtx_2_vol = sphere_sector_vol(p_2, p_3, p_0, p_1, rdW_2);
    double const vtx_3_vol = sphere_sector_vol(p_1, p_2, p_3, p_0, rdW_1);

    return (vtx_0_vol + vtx_1_vol + vtx_2_vol + vtx_3_vol);
}
/////////////
// Own data structures
/////////////

TCavity::TCavity(
    Molecule const &molecule, CellFilteringOptions const cell_opts) {

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

void TCavity::add_inner_cell(TTetrahedron const &cell, TetraInfo const &info) {
    _volume += volume(cell) - occupied_cell_vol(cell, info);
    _inner_cells.push_back(cell);
    _inner_info.push_back(info);
    return;
}

void TCavity::add_outer_cell(TTetrahedron const &cell, TetraInfo const &info) {
    _outer_cells.push_back(cell);
    _inner_info.push_back(info);
    return;
}

void TCavity::add_border_tetra(Point const &p0, Point const &ip1,
    Point const &ip2, Point const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - sphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    _tetra_border.emplace_back(p0, ip1, ip2, ip3);
    return;
}

void TCavity::add_border_penta(Point const &p0, Point const &p1,
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

    _penta_border.emplace_back(p0, p1, ip2, ip3, ip4, ip5);

    return;
}

void TCavity::add_border_penta(Point const &p0, Point const &p1,
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

    _penta_border.emplace_back(p0, p1, p2, ip3, ip4, ip5);

    return;
}

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(TTetrahedron const &cell, TetraInfo const &info) {

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

} // namespace ANA