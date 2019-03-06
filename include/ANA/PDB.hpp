#ifndef ANA_PDB_H
#define ANA_PDB_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Primitives.hpp>

namespace ANA {

void draw(Point const &punto, FILE *out_file, int &idx, int const resid,
    std::string const &name);

void draw_lines(CPoint const &punto, FILE *out_file, int constidx,
    int const resid, std::string const &name);

void draw_lines(CTriangle const &t, FILE *out_file, int &idx, int &resid);

void draw_lines(Triangle const &t, FILE *out_file, int &idx, int &resid);

void draw_lines(
    Finite_cells_iterator const cell, FILE *out_file, int &idx, int &resid);

template <class T>
void draw_polyhedron_lines(
    T const &poly, FILE *out_file, int &idx, int &resid) {

    for (auto const &each : poly._data) {
        draw_lines(each, out_file, idx++, resid, "POL");
    }
    ++resid;
    return;
}

void draw_lines(CConvexHull const &CH, std::string const &filename);

void draw_lines(CCavity const &hueco, std::string const &filename);

void connect_triangle(FILE *out_file, int const first_t, int const last_t);

void connect_tetrahedra(
    FILE *out_file, int const first_tetra, int const last_tetra);

void connect_pentahedra(
    FILE *out_file, int const first_penta, int const last_penta);

} // namespace PDB
// namespace ANA
#endif