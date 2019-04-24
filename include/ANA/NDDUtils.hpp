#ifndef ANA_NDD_UTILS_H
#define ANA_NDD_UTILS_H
#include <ANA/Cavity.hpp>
#include <ANA/ConvexHull.hpp>
#include <ANA/ConvexHullFunctions.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Modes.hpp>
#include <ANA/Options.hpp>
#include <ANA/PrimitivesUtils.hpp>
#include <ANA/TextUtils.hpp>

// DELETE after debugging.
#include <ANA/PDB.hpp>

namespace ANA::NDD {

// On-site NDD.
void ndd(Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts);

auto initialize_scaling_factors(NDDOptions const &NDD_opts, std::size_t j)
    -> std::vector<double>;

} // namespace NDD:: ANA

#endif // _H
