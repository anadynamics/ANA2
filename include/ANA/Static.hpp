#ifndef ANA_STATIC_H
#define ANA_STATIC_H
#include <ANA/ConvexHullFunctions.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>
#include <ANA/Read.hpp>
#include <ANA/Utils.hpp>
#include <ANA/Write.hpp>

namespace ANA {

// Main function for static version of ANA.
int static_ANA(InOutOptions const &io_opts, std::string &AA_indices_proto,
    std::string const &ASA_method, std::string const &only_side_ASA,
    std::string &exclude_ca_for_ASA_indices_proto, std::string const &list_wall,
    std::string const &list_wall_separator, std::string const &clusters_method,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    std::string &out_filename, std::string const &out_type,
    bool const triangulate_only_included_aas, bool const atom_only,
    CellFilteringOptions const cell_opts, double const max_probe,
    double const max_probe_length, int const sphere_count,
    int const nbr_of_vertices_to_include, int const clusters_min_size,
    int const precision);

} // namespace ANA

#endif
