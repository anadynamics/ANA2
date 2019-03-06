#ifndef ANA_MOLECULE_H
#define ANA_MOLECULE_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>
#include <ANA/SoAPrimitives.hpp>

namespace ANA {

class Molecule {
public:
    Molecule() = default;

    Molecule(std::string const &pdb_filename, bool const atom_only);

    int _natoms, _nres;
    std::vector<std::pair<CPoint, VertexInfo>> _data;
    // Alpha carbon indices
    std::vector<int> _alphaCarbons;
    // Hetero-atoms indices
    std::vector<int> _hetatms;
};

} // namespace ANA

#endif // _H