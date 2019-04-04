#ifndef ANA_CAVITY_H
#define ANA_CAVITY_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/PrimitivesUtils.hpp>

namespace ANA {

class Cavity {
public:
    Cavity() = default;

    Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

    void add_inner_cell(Finite_cells_iterator const &cell);

    void add_border_tetra(CPoint const &p0, CPoint const &ip1,
        CPoint const &ip2, CPoint const &ip3, double const vdw0);

    void add_border_penta(CPoint const &p0, CPoint const &p1, CPoint const &ip2,
        CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
        double const vdw0, double const vdw1);

    void add_border_penta(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
        double const vdw0, double const vdw1, double const vdw2);

    friend void write_PDB(Cavity const &hueco, std::string const &filename);

    std::vector<Finite_cells_iterator> _all_cells, _inner_cells, _outer_cells;
    double _volume = 0, _outer_volume = 0;

private:
    Delaunay _triangulation;
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<CCTetrahedron> _tetra_border;
    std::vector<CTriangularPrism> _penta_border;
    // Number of polyhedron vertices. Needed for output.
    std::vector<int> _poly_vtx_cnt;
};

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(Finite_cells_iterator const cell_iterator);

class TCavity {
public:
    TCavity() = default;

    TCavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

    void add_inner_cell(TTetrahedron const &cell, TetraInfo const &info);

    void add_outer_cell(TTetrahedron const &cell, TetraInfo const &info);

    void add_border_tetra(Point const &p0, Point const &ip1, Point const &ip2,
        Point const &ip3, double const vdw0);

    void add_border_penta_2_2(Point const &p0, Point const &p1,
        Point const &ip2, Point const &ip3, Point const &ip4, Point const &ip5,
        double const vdw0, double const vdw1);

    void add_border_penta_3_1(Point const &p0, Point const &p1, Point const &p2,
        Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
        double const vdw1, double const vdw2);

    friend void write_PDB(TCavity const &hueco, std::string const &filename);

    std::vector<TTetrahedron> _all_cells, _inner_cells, _outer_cells;
    std::vector<TetraInfo> _all_info, _inner_info, _outer_info;
    double _volume = 0, _outer_volume = 0;

private:
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<TTetrahedron> _tetra_border;
    std::vector<TTriangularPrism> _penta_border_2_2, _penta_border_3_1;
    // Number of polyhedron vertices. Needed for output.
    std::vector<int> _poly_vtx_cnt;
};

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(TTetrahedron const &cell, TetraInfo const &info);

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

} // namespace ANA

#endif // _H