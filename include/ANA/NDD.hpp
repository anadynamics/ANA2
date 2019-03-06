#ifndef ANA_NDD_H
#define ANA_NDD_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/ConvexHullFunctions.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Modes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/NDDUtils.hpp>
#include <ANA/Options.hpp>
#include <ANA/PDB.hpp>
#include <ANA/Primitives.hpp>
#include <ANA/Read.hpp>
#include <ANA/Utils.hpp>
#include <ANA/Write.hpp>

namespace ANA {

// Main function for NDD version of ANA.
int NDD_ANA(std::string const &in_filename, IncludedAreaOptions &IA_opts,
    NDDOptions const &NDD_opts, CellFilteringOptions const cell_opts,
    bool const atom_only);

} // namespace ANA

#endif // _H
