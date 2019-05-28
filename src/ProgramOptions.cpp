#include <ANA/ProgramOptions.hpp>

namespace ANA {

int get_parameters(int ac, char *av[], ANA::InOutOptions &io_opts,
    ANA::IncludedAreaOptions &IA_opts, std::string &AA_indices_proto,
    bool &triangulate_only_included_aas, bool &atom_only, int &precision,
    int &clusters_min_size, int &nbr_of_vertices_to_include, int &md_start,
    int &md_step, int &md_end, CellFilteringOptions &cell_opts,
    double &max_probe, double &max_probe_length, int &sphere_count,
    std::string &list_wall, std::string &list_wall_separator,
    std::string &clusters_method, std::string &only_side_ASA,
    std::string &ASA_method, std::string &exclude_ca_for_ASA,
    NDDOptions &NDD_opts, std::string &tool_check_CH,
    std::string &tool_pdb_to_ch, std::string &tool_pdb_norm,
    std::string &tool_aa_to_ca) {
    // clang-format off
    
    try { 
    std::string config_filename;
    // Declare options that are allowed only in the command line
    PO::options_description CLI_only_options(help_header);

    CLI_only_options.add_options()
///
// ## Terminal options
///
    ("input_struct,s",
    PO::value<std::string>(&io_opts._in_filename)->default_value("none"),
    "Input structure. This is the only positional argument so you can type it first and skip the flag.\n")

    ("input_md,d",
    PO::value<std::string>(&io_opts._in_md_filename)->default_value("none"),
    "Input file with MD simulation.\n")

    ("config_file,c",
    PO::value<std::string>(&config_filename)->default_value("none"),
    "Filename of the configuration file.\n")
    
    ("include,i",
    PO::value<std::string>(&IA_opts._filename)->default_value("none")->composing(),
    "Coordinates of the included area in PDB format.\n")

    ("output_draw,f",
    PO::value<std::string>(&io_opts._out_pdb_filename)->default_value("none")->composing(),
    "PDB output filename.\n")
    
    ("output_vol,o",
    PO::value<std::string>(&io_opts._out_vol_filename)->default_value("none")->composing(),
    "Volume output filename.\n")

    ("output_wall,w",
    PO::value<std::string>(&io_opts._out_wall_filename)->default_value("none")->composing(),
    "Cavity's wall atoms/residues output filename.\n")

    ("NDD_modes,M",
    PO::value<std::string>(&NDD_opts._modes_ndd_filename)->default_value("none")->composing(),
    "Input vectors for Non-Delaunay dynamics (NDD).\n")
    
    ("NDD_frequencies,F",
    PO::value<std::string>(&NDD_opts._freqs_ndd_filename)->default_value("none")->composing(),
    "Input frequencies (in cm^-1) to calculate the flexibility index by NDD.\n")

    ("NDD_scaling,S",
    PO::value<std::string>(&NDD_opts._scaling_ndd_filename)->default_value("none")->composing(),
    "Input scaling factors for NDD. See ANA's manual to understand why they are necessary.\n")

    ("NDD_size,Z",
    PO::value<int>(&NDD_opts._size)->default_value(1),
    "Scaling magnitude for the input vectors for NDD. Default: 1.\n")

    ("NDD_output,O",
    PO::value<std::string>(&NDD_opts._out_ndd_filename)->default_value("none")->composing(),
    "Suffix for the NDD output file. The scaling magnitude will be the prefix.\n")

    ("tool_check_CH,t",
    PO::value<std::string>(&tool_check_CH)->default_value("none")->composing(),
    "Filename of output PDB displaying the included area.\n")
    
    ("tool_pdb_to_ch,p",
    PO::value<std::string>(&tool_pdb_to_ch)->default_value("none")->composing(),
    "Read the input PDB and write 'include.ANA' file with the vertices of its convex hull.\n")
    
    ("tool_pdb_norm,n",
    PO::value<std::string>(&tool_pdb_norm)->default_value("none")->composing(),
    "Read the input PDB and renumber its atoms and residues. "
    "Be aware that his tool won't fix every error in your PDB.\n"
    "Writes the output PDB to \"tool_pdb_norm\".\n")
    
    ("tool_aa_to_ca,a",
    PO::value<std::string>(&tool_aa_to_ca)->default_value("none")->composing(),
    "Write Calpha atoms indices for the included residues to stdout.\n")
    
    ("ver,v", "Output version number.\n")
    
    ("help,h", "Output help message.\n");

    // Input structures is a positional argument
    PO::positional_options_description input_struct;
    input_struct.add("input_struct", 1);

    // Declare options that are only allowed in the config file
    PO::options_description cfg_only_options(help_header);
    cfg_only_options.add_options()

///
// ## Configuration file options
///
///
// ### Included area options
///
    ("included_area_residues",
    PO::value<std::string>(&IA_opts._resn_proto)->default_value("none"),
    "Amino acids that delimit the convex hull of the included area. "
    "Numbers from 1 to 9999, separated by commas, spaces or tabs.\n")

    ("included_area_atoms",
    PO::value<std::string>(&IA_opts._atom_proto)->default_value("none"),
    "Atoms that delimit the convex hull of the included area."
    "Numbers from 1 to 9999, separated by commas, spaces or tabs.\n")

    ("included_area_precision",
    PO::value<int>(&precision)->default_value(0),
    "0: keep all cells that intersect in the included area.\n"
    "1: prune the intersecting cells to keep only the voids that lie inside the included area.\n"
    "Default: 0.\n")

    ("sphere", 
    PO::value<std::string>(&IA_opts._sphere_proto)->default_value("none")->composing(),
    "Read the input coordinates and write 'include_sphere.ANA' file with the requested pseudo sphere.\n")
    ("cylinder",
    PO::value<std::string>(&IA_opts._cylinder_proto)->default_value("none")->composing(),
    "Read the input coordinates and write 'include_cylinder.ANA' file with the requested pseudo cylinder.\n")
    ("prism",
    PO::value<std::string>(&IA_opts._prism_proto)->default_value("none")->composing(),
    "Read the input coordinates and write 'include_prism.ANA' file with the requested prism.\n")

///
// ### Useful options for outlining the desired cavity with ANA Static
///

    ("included_residues",
    PO::value<std::string>(&AA_indices_proto)->default_value("none"),
    "Amino acids that line the desired cavity. Useful for discovering new cavities.\n")

    ("minimum_number_of_vertices_to_include",
    PO::value<int>(&nbr_of_vertices_to_include)->default_value(2),
    "Minimum number of wall atoms of the included amino acids.\n"
    "Default: 2.\n")

///
// ### Clusters options
///

    ("clusters_method",
    PO::value<std::string>(&clusters_method)->default_value("boxes"),
    "Define the method to clusterize cells.\n"
    "none: don't group null areas.\n"
    "facets: group null areas if the tetrahedrons share facets.\n"
    "boxes: group null areas if the tetrahedrons bounding boxes intersect.\n"
    "Default: boxes.\n")

    ("clusters_min_size",
    PO::value<int>(&clusters_min_size)->default_value(2),
    "Minimum number of cells a cluster has to have to be included.\n"
    "Default: 2.\n")

///
// #### Accessible Surface Area options
///

    ("ASA_discard_method",
    PO::value<std::string>(&ASA_method)->default_value("cm"),
    "Method for discarding cells according to solvent exposure.\n"
    " none: don't.\n"
    "cm: determine cells surrondings using the global CM as reference.\n"
    "backbone: draw a convex hull between the Calphas and discard every cell with its centroid outside(inside) the hull.\n"
    "axes: determine cells surrondings using its centroid and cartesian axes as references.\n"
    "Default: cm.\n")

    ("ASA_only_side",
    PO::value<std::string>(&only_side_ASA)->default_value("inside"),
    "Which side to keep.\n"
    "inside: keep inside nulls.\n"
    "outside: keep outside nulls.\n" 
    "Default: inside.\n")

    ("ASA_exclude_residues",
    PO::value<std::string>(&exclude_ca_for_ASA)->default_value("none"),
    "Residues to exclude from the ASA algorithms.\n")

    ("ASA_min_dot_pdt",
    PO::value<double>(&max_probe)->default_value(0.7)->composing(),
    "Minimum dot product to keep(discard) cell. As it increases more cells are classified as being \"outside\".\n"
    "Default: 0.7.\n")

    ("ASA_max_dist",
    PO::value<double>(&max_probe_length)->default_value(15)->composing(),
    "Maximum distance between cell and Calpha used to keep cell. The bigger this number is the better (but slower) the process gets.\n"
    "Default: 15.\n")
///
// ### MD options
///
    ("start",
    PO::value<int>(&md_start)->default_value(1),
    "Frame to begin reading at.\n"
    "Default: 1.\n")

    ("step",
    PO::value<int>(&md_step)->default_value(1),
    "Step count to read trajectory frames. NetCDF format only.\n"
    "Default: 1.\n")

    ("stop",
    PO::value<int>(&md_end)->default_value(0),
    "Frame to stop reading. Reads the whole trajectory when set to 0.\n"
    "Default: 0.\n")

///
// ### NDD options
///
    ("NDD_step",
    PO::value<int>(&NDD_opts._step)->default_value(3),
    "Number of steps in the NDD pipeline ANA will perform.\n"
    "1: ANA will not perform the derivative and instead output 2 files with the volumes of the displaced cavity in the positive and the negative direction.\n"
    "2: ANA will output the VGV.\n"
    "3: ANA will output the flexibility index.\n"
    "Default: 3.\n")

    ("NDD_modes_format",
    PO::value<std::string>(&NDD_opts._modes_format)->default_value("row"),
    "Format of the input vectors.\n"
    "amber: vectors will be read as Amber PCA modes.\n"
    "row: vectors will be read in row major order.\n"
    "column: vectors will be read in column major order.\n"
    "Default: row.\n")

    ("NDD_frequences_scaling",
    PO::value<bool>(&NDD_opts._scale_w_freqs)->default_value(false),
    "If true ANA will scale the input vectors by their frequencies instead of automatically generating its own scaling factors.\n"
    "Default: false.\n")

///
// ### Output options
///
    ("output_type",
    PO::value<std::string>(&io_opts._out_type)->default_value("grid_pdb")->composing(),
    "How to visualize output cavities."
    "raw_pdb: null areas as tetrahedrons.\n"
    "raw_cgo: null areas as tetrahedrons formed with CGO lines in a pymol script.\n"
    "grid_pdb: null areas filled with points.\n"
    "grid_cgo: null areas filled with CGO spheres in a pymol script.\n"
    "Default: grid_pdb.\n")

    ("sphere_count",
    PO::value<int>(&sphere_count)->default_value(5)->composing(),
    "Sphere count for grid output.\n"
    "Default: 5.\n")
///
// ### Misc options
///
    ("min_vol_radius",
    PO::value<double>(&cell_opts._minVR)->default_value(1.4)->composing(),
    "Radius of the sphere with the minimum volume to be taken into account.\n"
    "Default: 1.4.\n")

    ("max_area_radius",
    PO::value<double>(&cell_opts._maxSR)->default_value(99)->composing(),
    "Radius of the sphere with the maximum surface to be taken into account.\n"
    "Default: 99.\n")

    ("atom_only",
    PO::value<bool>(&atom_only)->default_value(true),
    "Triangulate only ATOM records.\n"
    "Default: true.\n")

    ("list_wall",
    PO::value<std::string>(&list_wall)->default_value("none"),
    "atom: output a file with the wall atoms.\n"
    "residue: output a file with the wall amino acids.\n" 
    "none: don't.\n"
    "Default: none.\n")

    ("separator", PO::value<std::string>(&list_wall_separator)->default_value("\t"),
    "Separator character for wall residues (atoms). \"+\": might useful for pymol.\n"
    "Default: \t.\n")
    
    ("triangulate_only_included_aas",
    PO::value<bool>(&triangulate_only_included_aas)->default_value(false),
    "Instead of triangulating the whole molecule triangulate only the included amino acids.\n"
    "Default: false.\n");

    // Now, map the variables
    PO::variables_map vm;
    PO::store(PO::command_line_parser(ac, av)
                  .options(CLI_only_options)
                  .positional(input_struct)
                  .run(),
              vm);
    if (vm.count("help")) {
      cout << CLI_only_options << "\n";
      return 1;
    }
    if (vm.count("ver")) {
      cout << "ANA v0.6" << "\n";
      return 1;
    }

    // Parse variables to see if a config file was specificated. If so, read it,
    // and parse again.
    PO::notify(vm);
    std::ifstream ifs(config_filename);
    store(PO::parse_config_file(ifs, cfg_only_options), vm);
    // Throws an error in case the first argument is not specified
    PO::notify(vm);

  } catch (exception &e) {
    cerr << "error: " << e.what() << "\n\n";
    return 1;
  } catch (...) {
    cerr << "Exception of unknown type when reading config file.\n\n";
  }

  if(io_opts._in_filename == "none" && tool_pdb_to_ch == "none" &&
     tool_pdb_norm == "none"         && tool_aa_to_ca == "none") {
    std::cerr << "Error: option '--input_struct' is required but missing."
    << "\n\n";
    return 1;
  }

  if (triangulate_only_included_aas == true && AA_indices_proto == "none") {
      std::cerr << "ERROR: 'triangulate_only_included_aas' is TRUE, but "
      "'included_amino_acids' was not set." << "\n\n";
      return 1;
  }

  if (IA_opts._resn_proto != "none" && IA_opts._atom_proto != "none") {
    IA_opts._atom_proto = "none";
    std::cerr << "Input warning: Both 'included_area_residues' and "
    "'included_area_atoms' were set. Using the former. " << "\n\n";
  }


  bool const defined_included_area = (IA_opts._resn_proto != "none") || 
      (IA_opts._atom_proto != "none") || (IA_opts._sphere_proto != "none") || 
      (IA_opts._cylinder_proto != "none") || (IA_opts._prism_proto != "none") || 
      (IA_opts._filename != "none");
  bool const ndd_requested = NDD_opts._modes_ndd_filename != "none";
  bool const modes_or_ndd = ndd_requested || (io_opts._in_md_filename != "none");

  if (!defined_included_area && modes_or_ndd) {
    std::cerr << "Input error: You are running ANA MD/NDD without an included "
    "area. Check ANA's manual." << "\n\n";
    return 1;
  }

  if( ((io_opts._out_wall_filename != "none") && (list_wall == "none")) ||
      ((io_opts._out_wall_filename == "none") && (list_wall != "none")) ) {
    std::cerr << "Input warning: both \"output_wall,w\"  and \"list_wall\" options " 
    " need to be defined to obtain the cavity lining residues/atoms. ANA will not " 
    " write an \"output_wall,w\" file." << "\n\n";
  }

  if(ndd_requested){
      bool const no_frequencies = (NDD_opts._freqs_ndd_filename == "none") && 
        (NDD_opts._modes_format != "amber");
      if ((NDD_opts._step == 3) && no_frequencies) {
        std::cerr << "Input error: You are trying to get a cavity's flexibility index "
        "without the eigenvectors frequencies. Check ANA's manual." << "\n\n";
        return 1;
      }

      if ((NDD_opts._step < 3) && NDD_opts._out_ndd_filename == "none") {
        std::cerr << "Input error: NDD requested but no NDD_output filename was not set." << "\n\n";
        return 1;
      }
  
      if (NDD_opts._scale_w_freqs && no_frequencies) {
        std::cerr << "Input error: NDD_frequences_scaling is true but no "
        "eigenvectors frequencies could be found." << "\n\n";
        return 1;
      }

      if (NDD_opts._scaling_ndd_filename != "none" && NDD_opts._scale_w_freqs) {
        std::cerr << "Input warning: NDD_frequences_scaling is true and input file "
        "NDD_scaling was also set. The later will override the former." << "\n\n";
      }

      if ((NDD_opts._freqs_ndd_filename != "none") && 
      ((NDD_opts._modes_format != "row") && (NDD_opts._modes_format != "column"))) {
        std::cerr << "Input error: NDD_modes_format should be set to row or column." << "\n\n";
        return 1;
      }
  }
  cell_opts.update();
  
  return 0;
}

} // namespace ANAs
