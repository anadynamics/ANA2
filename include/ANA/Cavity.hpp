#ifndef ANA_CAVITY_H
#define ANA_CAVITY_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>

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

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

} // namespace ANA

#endif // _H