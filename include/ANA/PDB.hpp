#ifndef ANA_PDB_H
#define ANA_PDB_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/PrimitivesUtils.hpp>

namespace ANA {

// CGAL types.

void write_PDB(Delaunay const &T, std::string const &filename);

auto draw_lines(Cell_iterator const cell, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int>;

auto draw(CPoint const &punto, FILE *out_file, std::pair<int, int> idx_resid,
    std::string const &name) -> int;

// Own types.

// TODO: replace `std::string const &` with `std::string_view` as soon as a new
// string_view supporting fmt version comes out.
void write_PDB(ConvexHull const &CH, std::string const &filename);

void write_PDB(Cavity const &hueco, std::string const &filename);

[[nodiscard]] auto draw_lines(Triangle const &t, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int>;

[[nodiscard]] auto draw_lines(Tetrahedron const cell, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int>;

template <class T>
[[nodiscard]] auto draw_polyhedron(T const &poly, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int> {

    for (auto const &punto : poly._data) {
        idx_resid.first = draw(punto, out_file, idx_resid, resname);
    }
    ++idx_resid.second;

    return idx_resid;
}

[[nodiscard]] auto draw(Point const &punto, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &name) -> int;

void connect_triangle(FILE *out_file, int const first_t, int const last_t);

void connect_tetrahedra(
    FILE *out_file, int const first_atm, int const last_atm);

void connect_pentahedra_2_2(
    FILE *out_file, int const first_atm, int const last_atm);

void connect_pentahedra_3_1(
    FILE *out_file, int const first_atm, int const last_atm);
} // namespace ANA
#endif