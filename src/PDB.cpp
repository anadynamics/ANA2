#include <ANA/PDB.hpp>

namespace ANA {

void write_PDB(ConvexHull const &CH, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        std::pair<int, int> idx_resid{1, 1};
        for (auto const &triangle : CH._triangles) {
            idx_resid = draw_lines(triangle, out_file, idx_resid, " IA");
        }
        connect_triangle(out_file, 1, idx_resid.first);
    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);

    return;
}

void write_PDB(Cavity const &hueco, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        std::pair<int, int> idx_resid{1, 1};

        for (auto const &cell : hueco._inner_cells) {
            idx_resid = draw_lines(cell, out_file, idx_resid, "CEL");
        }

        for (auto const &tetra : hueco._tetra_border) {
            idx_resid = draw_polyhedron(tetra, out_file, idx_resid, "TET");
        }
        int const pre_penta_2_2 = idx_resid.first - 1;

        for (auto const &penta : hueco._penta_border_2_2) {
            idx_resid = draw_polyhedron(penta, out_file, idx_resid, "P22");
        }
        int const pre_penta_3_1 = idx_resid.first - 1;

        for (auto const &penta : hueco._penta_border_3_1) {
            idx_resid = draw_polyhedron(penta, out_file, idx_resid, "P31");
        }

        connect_tetrahedra(out_file, 1, pre_penta_2_2);
        connect_pentahedra_2_2(out_file, pre_penta_2_2 + 1, pre_penta_3_1);
        connect_pentahedra_3_1(
            out_file, pre_penta_3_1 + 1, idx_resid.first - 1);

    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);

    return;
}

auto draw_lines(Triangle const &t, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int> {

    idx_resid.first = draw(t[0], out_file, idx_resid, resname);
    idx_resid.first = draw(t[1], out_file, idx_resid, resname);
    idx_resid.first = draw(t[2], out_file, idx_resid, resname);
    ++idx_resid.second;

    return idx_resid;
}

auto draw_lines(Tetrahedron const cell, FILE *out_file,
    std::pair<int, int> idx_resid, std::string const &resname)
    -> std::pair<int, int> {

    idx_resid.first = draw(cell[0], out_file, idx_resid, resname);
    idx_resid.first = draw(cell[1], out_file, idx_resid, resname);
    idx_resid.first = draw(cell[2], out_file, idx_resid, resname);
    idx_resid.first = draw(cell[3], out_file, idx_resid, resname);
    ++idx_resid.second;

    return idx_resid;
}

auto draw(Point const &punto, FILE *out_file, std::pair<int, int> idx_resid,
    std::string const &name) -> int {
    fmt::print(out_file,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx_resid.first++, "H", name, "A", idx_resid.second, punto[0],
        punto[1], punto[2], 1.0, 0.0, "H");
    return idx_resid.first;
}

void connect_triangle(FILE *out_file, int const first_t, int const last_t) {

    for (auto i = first_t; i < last_t; i += 3) {
        auto const j = i + 1;
        auto const k = i + 2;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
        fmt::print(out_file, "CONECT {:>4} {:>4}\n", j, k);
    }
    return;
}

void connect_tetrahedra(
    FILE *out_file, int const first_atm, int const last_atm) {
    assert(first_atm <= last_atm);

    for (auto i = first_atm; i < last_atm; i += 4) {
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", j, k, l, i);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, l, i, j);
    }
    return;
}

void connect_pentahedra_2_2(
    FILE *out_file, int const first_atm, int const last_atm) {
    assert(first_atm <= last_atm);

    for (auto i = first_atm; i < last_atm; i += 6) {
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;
        auto const m = i + 4;
        auto const n = i + 5;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", n, j, l, m);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", m, j, k);
        fmt::print(out_file, "CONECT {:>4} {:>4}\n", k, l);
    }
    return;
}

void connect_pentahedra_3_1(
    FILE *out_file, int const first_atm, int const last_atm) {
    assert(first_atm <= last_atm);

    for (auto i = first_atm; i < last_atm; i += 6) {
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;
        auto const m = i + 4;
        auto const n = i + 5;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", n, m, l, k);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", m, l, j);
        fmt::print(out_file, "CONECT {:>4} {:>4}\n", j, k);
    }
    return;
}

} // namespace ANA