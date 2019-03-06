#ifndef ANA_CAVITY_H
#define ANA_CAVITY_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

class CCavity {
public:
    CCavity() = default;

    CCavity(Delaunay T) noexcept : _triangulation(T) {}

    CCavity(Delaunay &&T) noexcept : _triangulation(T) {}

    CCavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

    void cadd_border_tetra(CPoint const &p0, CPoint const &ip1,
        CPoint const &ip2, CPoint const &ip3, double const vdw0);

    void cadd_border_penta(CPoint const &p0, CPoint const &p1,
        CPoint const &ip2, CPoint const &ip3, CPoint const &ip4,
        CPoint const &ip5, double const vdw0, double const vdw1);

    void cadd_border_penta(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &ip3, CPoint const &ip4, CPoint const &ip5,
        double const vdw0, double const vdw1, double const vdw2);

    friend void write_PDB(CCavity const &hueco, std::string const &filename);

    Delaunay _triangulation;
    std::vector<Finite_cells_iterator> _all_cells, _included_cells;
    double _volume = 0, _outer_volume = 0;

private:
    // Border polyhedrons from the cells that intersected the convex hull.
    std::vector<CCTetrahedron> _tetra_border;
    std::vector<CTriangularPrism> _penta_border;
    // Number of polyhedron vertices. Needed for output.
    std::vector<int> _poly_vtx_cnt;
};

// class Cavity {
// public:
//     Cavity() = default;

//     Cavity(Delaunay T) noexcept : _triangulation(T) {}

//     Cavity(Delaunay &&T) noexcept : _triangulation(T) {}

//     Cavity(Molecule const &molecule, CellFilteringOptions const cell_opts);

//     void add_border_tetra(Point const &p0, Point const &ip1, Point const
//     &ip2,
//         Point const &ip3, double const vdw0);

//     void add_border_penta(Point const &p0, Point const &p1, Point const &ip2,
//         Point const &ip3, Point const &ip4, Point const &ip5, double const
//         vdw0, double const vdw1);

//     void add_border_penta(Point const &p0, Point const &p1, Point const &p2,
//         Point const &ip3, Point const &ip4, Point const &ip5, double const
//         vdw0, double const vdw1, double const vdw2);

//     friend void draw_lines(Cavity const &hueco, std::string const &filename);

//     Delaunay _triangulation;
//     std::vector<Tetrahedron> _all_cells, _included_cells;
//     VtxInfo _all_info, _included_info;
//     double _volume = 0, _outer_volume = 0;

// private:
//     // Border polyhedrons from the cells that intersected the convex hull.
//     std::vector<Tetrahedron> _tetra_border;
//     std::vector<TriangularPrism> _penta_border;
//     // Number of polyhedron vertices. Needed for output.
//     std::vector<int> _poly_vtx_cnt;
// };

// Tool for parsing a double from input file stringstream
double parse_double(std::stringstream &in_stream);

} // namespace ANA

#endif // _H