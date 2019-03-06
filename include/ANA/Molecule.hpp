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

    Molecule(std::string const &filename, bool const atom_only);

    int _natoms, _nres;
    std::vector<std::pair<CPoint, VertexInfo>> _data;
    // Alpha carbon indices
    std::vector<int> _alphaCarbons;
    // Hetero-atoms indices
    std::vector<int> _hetatms;
};

class SoAMolecule {
public:
    SoAMolecule() = default;

    SoAMolecule(std::string const &filename, bool const atom_only);

    // Residue count.
    int _nres;
    // xyz atoms coordinates.
    Points _data;
    // Atom's index, Van der Waals radius, and its residue name and index.
    VtxInfo _info;
    // Alpha carbon indices.
    std::vector<int> _alphaCarbons;
    // Hetero-atoms indices.
    std::vector<int> _hetatms;
};

} // namespace ANA

#endif // _H