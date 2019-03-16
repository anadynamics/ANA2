#include <ANA/PDB.hpp>

namespace ANA {

void write_PDB(ConvexHull const &CH, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        std::pair<int, int> idx_resid {1, 1};
        for (auto const &triangle : CH._triangles) {
            idx_resid = draw_lines(triangle, out_file, idx_resid);
        }
        connect_triangle(out_file, 1, idx_resid.second);
    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);

    return;
}

void write_PDB(Cavity const &hueco, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        std::pair<int, int> idx_resid {1, 1};
        for (auto const &cell : hueco._inner_cells) {
            idx_resid = draw_lines(cell, out_file, idx_resid);
        }

        for (auto const &tetra : hueco._tetra_border) {
            idx_resid = draw_polyhedron(tetra, out_file, idx_resid);
        }

        int const pre_penta = idx_resid.second + 1;
        for (auto const &penta : hueco._penta_border) {
            idx_resid = draw_polyhedron(penta, out_file, idx_resid);
        }

        connect_tetrahedra(out_file, 1, pre_penta);
        connect_pentahedra(out_file, pre_penta, idx_resid.second);
    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);

    return;
}

auto draw(CPoint const &punto, FILE *out_file, std::pair<int, int> idx_resid,
    std::string const &name) -> int {
    fmt::print(out_file,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx_resid.first++, "H", name, "A", idx_resid.second,
        punto.x(), punto.y(), punto.z(), 1.0, 0.0, "H");
    return idx_resid.first;
}

auto draw_lines(CTriangle const &t, FILE *out_file,
    std::pair<int, int> idx_resid) -> std::pair<int, int> {

    idx_resid.first = draw(t.vertex(0), out_file, idx_resid, " IA");
    idx_resid.first = draw(t.vertex(1), out_file, idx_resid, " IA");
    idx_resid.first = draw(t.vertex(2), out_file, idx_resid, " IA");
    ++idx_resid.second;

    return idx_resid;
}

auto draw_lines(Finite_cells_iterator const cell, FILE *out_file,
    std::pair<int, int> idx_resid) -> std::pair<int, int> {

    idx_resid.first =
        draw(cell->vertex(0)->point(), out_file, idx_resid, "CEL");
    idx_resid.first =
        draw(cell->vertex(1)->point(), out_file, idx_resid, "CEL");
    idx_resid.first =
        draw(cell->vertex(2)->point(), out_file, idx_resid, "CEL");
    idx_resid.first =
        draw(cell->vertex(3)->point(), out_file, idx_resid, "CEL");
    ++idx_resid.second;

    return idx_resid;
}

void connect_triangle(FILE *out_file, int const first_t, int const last_t) {

    for (auto r = first_t; r < last_t; ++r) {
        auto const i = (r - 1) * 3 + 1;
        auto const j = i + 1;
        auto const k = i + 2;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, i, k);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, i, j);
    }
    return;
}

void connect_tetrahedra(
    FILE *out_file, int const first_tetra, int const last_tetra) {
    assert(first_tetra < last_tetra);

    for (auto r = first_tetra; r < last_tetra; ++r) {
        auto const i = (r - 1) * 4 + 1;
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", j, k, l, i);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, l, i, j);
        // fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", l, i, j, k);
    }
    return;
}

void connect_pentahedra(
    FILE *out_file, int const first_penta, int const last_penta) {
    assert(first_penta < last_penta);

    int const first = 4 * (first_penta - 1) + 1;
    for (auto p = first_penta; p < last_penta; ++p) {
        auto const i = 6 * (p - first_penta) + first;
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;
        auto const m = i + 4;
        auto const n = i + 5;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", j, k, l, i);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, l, i, j);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", l, k, m, n);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", m, n, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", n, k, l, m);
    }
    return;
}

} // namespace ANA