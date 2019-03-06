#include <ANA/Cavity.hpp>

namespace ANA {

CCavity::CCavity(
    Molecule const &molecule, CellFilteringOptions const cell_opts) {

    _triangulation.insert(molecule._data.begin(), molecule._data.end());

    Finite_cells_iterator fc_ite_end = _triangulation.finite_cells_end();
    for (auto fc_ite = _triangulation.finite_cells_begin();
         fc_ite != fc_ite_end; ++fc_ite) {
        if (volume(fc_ite) > cell_opts._min_CV) {
            _all_cells.push_back(fc_ite);
        }
    }
}

void CCavity::cadd_border_tetra(CPoint const &p0, CPoint const &ip1,
    CPoint const &ip2, CPoint const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - csphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    _tetra_border.emplace_back(p0, ip1, ip2, ip3);
    return;
}

void CCavity::cadd_border_penta(CPoint const &p0, CPoint const &p1,
    CPoint const &ip2, CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
    double const vdw0, double const vdw1) {

    double const vol1 = volume(p0, p1, ip2, ip3) -
        csphere_sector_vol(p0, p1, ip2, ip3, vdw0) -
        csphere_sector_vol(p1, p0, ip2, ip3, vdw1);

    double const vol2 =
        volume(p1, ip2, ip3, ip4) - csphere_sector_vol(p1, ip2, ip3, ip4, vdw1);

    double const vol3 =
        volume(p1, ip3, ip4, ip5) - csphere_sector_vol(p1, ip3, ip4, ip5, vdw1);
    _outer_volume += vol1 + vol2 + vol3;

    _penta_border.emplace_back(p0, p1, ip2, ip3, ip4, ip5);

    return;
}

void CCavity::cadd_border_penta(CPoint const &p0, CPoint const &p1,
    CPoint const &p2, CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
    double const vdw0, double const vdw1, double const vdw2) {

    double const vol1 = volume(p0, p1, p2, ip3) -
        csphere_sector_vol(p0, p1, p2, ip3, vdw0) -
        csphere_sector_vol(p1, p0, p2, ip3, vdw1) -
        csphere_sector_vol(p2, p1, p0, ip3, vdw2);

    double const vol2 = volume(p1, p2, ip3, ip4) -
        csphere_sector_vol(p1, p2, ip3, ip4, vdw1) -
        csphere_sector_vol(p2, p1, ip3, ip4, vdw2);

    double const vol3 =
        volume(p2, ip3, ip4, ip5) - csphere_sector_vol(p2, ip3, ip4, ip5, vdw2);

    _outer_volume += vol1 + vol2 + vol3;

    _penta_border.emplace_back(p0, p1, p2, ip3, ip4, ip5);

    return;
}

Cavity::Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts) {

    _triangulation.insert(molecule._data.begin(), molecule._data.end());

    Finite_cells_iterator fc_ite_end = _triangulation.finite_cells_end();
    for (auto fc_ite = _triangulation.finite_cells_begin();
         fc_ite != fc_ite_end; ++fc_ite) {
        if (volume(fc_ite) > cell_opts._min_CV) {
            _all_cells.emplace_back(fc_ite);

        }
    }
}

void Cavity::add_border_tetra(Point const &p0, Point const &ip1,
    Point const &ip2, Point const &ip3, double const vdw0) {

    double const vol =
        volume(p0, ip1, ip2, ip3) - sphere_sector_vol(p0, ip1, ip2, ip3, vdw0);
    _outer_volume += vol;

    _tetra_border.emplace_back(p0, ip1, ip2, ip3);
    return;
}

void Cavity::add_border_penta(Point const &p0, Point const &p1,
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

void Cavity::add_border_penta(Point const &p0, Point const &p1, Point const &p2,
    Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
    double const vdw1, double const vdw2) {

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

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream) {
    std::string temp;
    double coord = 0;

    in_stream >> temp;
    try {
        coord = std::stof(temp);
    } catch (const std::invalid_argument &ia) {
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