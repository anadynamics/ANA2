#include <ANA/Molecule.hpp>

namespace ANA {

Molecule::Molecule(std::string const &pdb_filename, bool const atom_only) {

    // Read PDB
    chemfiles::Trajectory input_pdb_traj(pdb_filename);
    auto const input_pdb_frame = input_pdb_traj.read();
    auto const in_xyz = input_pdb_frame.positions();
    auto const input_pdb_top = input_pdb_frame.topology();

    int const natoms = input_pdb_top.natoms();
    int const nres = input_pdb_top.residues().size();

    _data.reserve(natoms);
    _alphaCarbons.reserve(nres);

    for (auto const &residuo : input_pdb_top.residues()) {
        auto const res_name = residuo.name();
        for (auto const &i : residuo) {
            if (atom_only &&
                input_pdb_top[i].get("is_hetatm").value_or(false).as_bool()) {
                // Save the HETATM indices to discard them during MD and
                // NDD runs.
                _hetatms.push_back(i);
                continue;
            }

            auto const vdw_opt = input_pdb_top[i].vdw_radius();
            double vdw = 1.5;
            if (vdw_opt) {
                vdw = vdw_opt.value();
            } else {
                printf("Element for atom %i not available. Using Van Der Walls "
                       "radius of 1.5.\n",
                    static_cast<int>(i) + 1);
            }

            VertexInfo const vi1(i, vdw, residuo.id().value(), res_name);
            CPoint const p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
            _data.emplace_back(p1, vi1);

            if (input_pdb_top[i].name() == "CA") {
                // Will use to draw Convex Hull or whatever. It's always
                // useful.
                _alphaCarbons.push_back(i);
            }
        }
    }

    return;
}


} // namespace ANA