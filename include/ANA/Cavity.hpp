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

    void add_inner_cell(Tetrahedron const &cell, TetraInfo const &info);

    void add_outer_cell(Tetrahedron const &cell, TetraInfo const &info);

    void add_border_tetra(Point const &p0, Point const &ip1, Point const &ip2,
        Point const &ip3, double const vdw0);

    void add_border_penta_2_2(Point const &p0, Point const &p1,
        Point const &ip2, Point const &ip3, Point const &ip4, Point const &ip5,
        double const vdw0, double const vdw1);

    void add_border_penta_3_1(Point const &p0, Point const &p1, Point const &p2,
        Point const &ip3, Point const &ip4, Point const &ip5, double const vdw0,
        double const vdw1, double const vdw2);

    // Constructor for NDD. Returns an updated Cavity displacing the input
    // Cavity along the input vector scaled by the step_size.
    Cavity(Cavity const &hueco, std::vector<double> const &evector,
        double const step_size);

    // Auxiliary function for the NDD constructor
    void move_cells(std::vector<Tetrahedron> const &cells,
        std::vector<TetraInfo> const &info, std::vector<double> const &evector,
        double const step_size);

    friend void write_PDB(Cavity const &hueco, std::string const &filename);
    friend void debug_write_PDB_1(
        Cavity const &hueco, std::string const &filename);
    friend void debug_write_PDB_2(
        Cavity const &hueco, std::string const &filename);

    std::vector<Tetrahedron> _all_cells, _inner_cells, _outer_cells;
    std::vector<TetraInfo> _all_info, _inner_info, _outer_info;
    double _volume = 0, _outer_volume = 0;

private:
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<Tetrahedron> _tetra_border;
    std::vector<TriangularPrism> _penta_border_2_2, _penta_border_3_1;
    // Number of polyhedron vertices. Needed for output.
    std::vector<int> _poly_vtx_cnt;
};

// Get the volume filled by the cells 4 atoms.
double occupied_cell_vol(Tetrahedron const &cell, TetraInfo const &info);

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

void write_volume(Cavity const &hueco, std::string const &filename);

} // namespace ANA

#endif // _H