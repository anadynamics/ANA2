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

namespace ANA::NDD {

double constexpr cte = 5.0219059911006245e-5;
double constexpr TO_CM1 = 3.548143227025099e-4;

// On-site NDD.
void ndd(Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts,
    std::string const &pdb_filename);

auto initialize_scaling_factors(Modes const &modos, NDDOptions const &NDD_opts)
    -> std::vector<double>;

void write_result(NDDOptions const &NDD_opts, Modes const &modos,
    std::vector<double> const &neg_vols_ndd,
    std::vector<double> const &pos_vols_ndd, std::vector<double> const &vgv);

void barletta_index(Modes const &modos, NDDOptions const &NDD_opts,
    std::vector<double> const &vgv);

} // namespace NDD:: ANA

#endif // _H
