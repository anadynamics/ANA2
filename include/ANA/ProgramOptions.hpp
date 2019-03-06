#ifndef ANA_PROGRAM_OPTIONS_H
#define ANA_PROGRAM_OPTIONS_H
#include <ANA/Options.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace ANA {

namespace PO = boost::program_options;

using namespace std;

std::string const help_header = "\t\t\t----- ANA -----";

int get_parameters(int ac, char *av[], std::string &input_struct_filename,
    std::string &input_md_filename, ANA::IncludedAreaOptions &IA_opts,
    std::string &AA_indices_proto, bool &triangulate_only_included_aas,
    bool &atom_only, int &precision, int &clusters_min_size,
    int &nbr_of_vertices_to_include, int &md_start, int &md_step, int &md_end,
    CellFilteringOptions &cell_opts, double &max_probe,
    double &max_probe_length, int &sphere_count, std::string &list_wall,
    std::string &list_wall_separator, std::string &clusters_method,
    std::string &only_side_ASA, std::string &ASA_method,
    std::string &exclude_ca_for_ASA, NDDOptions &NDD_opts,
    std::string &out_filename, std::string &out_vol, std::string &output_type,
    std::string &tool_check_CH, std::string &tool_pdb_to_ch,
    std::string &tool_pdb_norm, std::string &tool_aa_to_ca);

} // namespace ANA

#endif // _H
