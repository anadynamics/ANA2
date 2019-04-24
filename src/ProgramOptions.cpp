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
    ("input_struct,i", PO::value<std::string>(&io_opts._in_filename)
	  ->default_value("none"), "Input structure (pdb). Positional argument.\n")
    ("input_md,d", PO::value<std::string>(&io_opts._in_md_filename)
    ->default_value("none"), "Input file with MD simulation.\n")
    ("NDD_modes,M", PO::value<std::string>(&NDD_opts._modes_ndd_filename)
    ->default_value("none")->composing(), "Input vectors for non Delaunay dynamics.\n")
    ("NDD_evals,E", PO::value<std::string>(&NDD_opts._evalues_ndd_filename)
    ->default_value("none")->composing(), "Input eigenvalues for non Delaunay dynamics.\n")
    ("NDD_scaling,S", PO::value<std::string>(&NDD_opts._scaling_ndd_filename)
    ->default_value("none")->composing(), "Input scaling factors for non Delaunay dynamics.\n")
    ("config_file,c", PO::value<std::string>(&config_filename)
    ->default_value("ANA.cfg"), "Filename of the configuration file. Default: \"ANA.cfg\".\n")
    ("output_draw,o", PO::value<std::string>(&io_opts._out_pdb_filename)
    ->default_value("none")->composing(),"Output filename.\n")
    ("out_vol,v", PO::value<std::string>(&io_opts._out_vol_filename)
    ->default_value("none")->composing(),"Volume output filename.\n")
    ("NDD_output,O", PO::value<std::string>(&NDD_opts._out_ndd_filename)
    ->default_value("ANA_NDD.out")->composing(),
    "Name of the non Delaunay dynamics output file. Default: \"ANA_NDD.out\".\n")
    ("include,f", PO::value<std::string>(&IA_opts._filename)
    ->default_value("none")->composing(),
    "Coordinates of the included area in PDB format.\n")

    ("tool_check_CH,C", PO::value<std::string>(&tool_check_CH)->default_value("none")
    ->composing(), "Filename of the CGO pymol script displaying the included "
    "area. If set, must be at least 3 characters long.\n")
    ("tool_pdb_to_ch,P", PO::value<std::string>(&tool_pdb_to_ch)->default_value("none")
    ->composing(), "Read the input PDB and write 'include.ANA' file with the "
    "vertices of its convex hull.\n")
    ("tool_pdb_norm,N", PO::value<std::string>(&tool_pdb_norm)->default_value("none")
    ->composing(), "Read the input PDB and renumber its atoms and residues."
    "Write the output PDB to \"tool_pdb_norm\".\n")
    ("tool_aa_to_ca,T", PO::value<std::string>(&tool_aa_to_ca)->default_value("none")
      ->composing(), "Write Calpha atoms indices for the included residues "
      "to stdout.\n")
    ("ver", "Output version number.\n")("help,h", "Output help message.\n");

    // Input structures is a positional argument
    PO::positional_options_description input_struct;
    input_struct.add("input_struct", 1);

    // Declare options that are only allowed in the config file
    PO::options_description cfg_only_options(help_header);
    cfg_only_options.add_options()

    ("included_amino_acids",
    PO::value<std::string>(&AA_indices_proto)->default_value("none"),
    "Amino acids that are part of a cell.\n")

    ("included_area_residues", PO::value<std::string>(&IA_opts._resn_proto)
    ->default_value("none"), "Amino acids that delimit the convex hull of the"
    " included area.\n")

    ("included_area_atoms", PO::value<std::string>(&IA_opts._atom_proto)
    ->default_value("none"), "Atoms that delimit the convex hull of the"
    " included area.\n")

    ("included_area_precision", PO::value<int>(&precision)
    ->default_value(0), "0: keep all cells that intercede in the included"
    " area. 1: only keep null areas that are within the included area. "
    "Default: 0.\n")

    ("sphere", PO::value<std::string>(&IA_opts._sphere_proto)->default_value("none")
    ->composing(), "Read the input coordinates and write 'include_sphere.ANA' "
    "file with the requested pseudo sphere.\n")
    ("cylinder", PO::value<std::string>(&IA_opts._cylinder_proto)->default_value("none")
    ->composing(), "Read the input coordinates and write 'include_cylinder.ANA' "
    "file with the requested pseudo cylinder.\n")
    ("prism", PO::value<std::string>(&IA_opts._prism_proto)->default_value("none")
    ->composing(), "Read the input coordinates and write 'include_prism.ANA' "
    "file with the requested prism.\n")

    ("triangulate_only_included_aas",
    PO::value<bool>(&triangulate_only_included_aas)
    ->default_value(false), "Instead of triangulating the whole molecule. "
    "Triangulate only the included amino acids. Default: false.\n")

    ("atom_only", PO::value<bool>(&atom_only)
    ->default_value(true), "Triangulate only ATOM records. Default: true.\n")

    ("clusters_method", PO::value<std::string>(&clusters_method)->default_value("boxes"),
    "\t none: don't group null areas.\n\t"
    "facets: group null areas if the tetrahedrons share facets.\n\t"
    "boxes: group null areas if the tetrahedrons bounding boxes intersect."
    "Default: \"boxes\"")

    ("clusters_min_size", PO::value<int>(&clusters_min_size)
    ->default_value(2), "Minimum number of cells a cluster has to have to be "
    "included. Default: 2. \n")

    ("minimum_number_of_vertices_to_include", PO::value<int>(
    &nbr_of_vertices_to_include)->default_value(2), "Minimum number of wall "
    "atoms of the included amino acids. Default: 2.\n")

    ("list_wall", PO::value<std::string>(&list_wall)->default_value("none"),
    "atom: output a file with the wall atoms; residue: output a file with the wall "
    "amino acids; none: don't.\n")

    ("separator", PO::value<std::string>(&list_wall_separator)->default_value("\t"),
    "Separator character for wall amino acids / atoms. \"+\": useful for pymol.\n"
    "Default: \"\t\".\n")


    ("ASA_discard_method", PO::value<std::string>(&ASA_method)
    ->default_value("cm"), "\t none: don't discard \n\t"
    "cm: determine cells surrondings using the global CM as reference.\n\t"
    "backbone: draw a convex hull between the Calphas and discard every cell "
    "with its centroid outside(inside) the hull.\n\t"
    "axes: determine cells surrondings using its centroid and cartesian "
    " axes as references. Default: \"cm\".\n")

    ("ASA_only_side", PO::value<std::string>(&only_side_ASA)
    ->default_value("inside"), "\t inside: keep inside nulls \n\t"
    " outside: keep outside nulls. Default: \"inside\".\n")

    ("ASA_exclude_amino_acids", PO::value<std::string>(&exclude_ca_for_ASA)
    ->default_value("none"), "Residues to exclude from the ASA algorithms.\n")

    ("ASA_min_dot_pdt", PO::value<double>(&max_probe)->default_value(0.7)
    ->composing(), "Minimum dot product to keep(discard) cell. As it increases "
    "more cells are classified as being \"outside\". Default: 0.7.\n")

    ("ASA_max_dist", PO::value<double>(&max_probe_length)->default_value(15)
    ->composing(), "Maximum distance between cell and Calpha used to keep cell."
    "The bigger this number is the better (but slower) the process gets."
    "Default: 15.\n")

    ("start", PO::value<int>(&md_start)->default_value(1),
    "Frame to begin reading at. Default: 1.\n")

    ("step", PO::value<int>(&md_step)->default_value(1),
    "Step count to read trajectory frames. NetCDF format only."
    "Default: 1.\n")

    ("stop", PO::value<int>(&md_end)->default_value(0),
    "Frame to stop reading. If set to 0, read all. Default: 0.\n")

    ("NDD_derivative", PO::value<bool>(&NDD_opts._derivative)->default_value(true),
    "If set to false, ANA will not perform the derivative and instead output"
    "2 files with the volumes of the displaced cavity in the positive and"
    "the negative direction. Default: true\n")

    ("NDD_step", PO::value<int>(&NDD_opts._step)->default_value(5),
    "Scaling number for input vectors in non Delaunay dynamics. Default: 5\n")

    ("NDD_modes_format", PO::value<std::string>(&NDD_opts._modes_format)->default_value("row"),
    "amber: vectors will be read as Amber PCA modes. "
    "row(column): vectors will be read in row(column) major order. Default: row\n")

    ("min_vol_radius", PO::value<double>(&cell_opts._minVR)->default_value(1.4)
    ->composing(), "Radius of the sphere with the minimum volume to be taken "
    "into account. Default: 1.4.\n")

    ("max_area_radius", PO::value<double>(&cell_opts._maxSR)->default_value(99)
    ->composing(), "Radius of the sphere with the maximum surface to be taken "
    "into account. Default: 99.\n")

    ("sphere_count", PO::value<int>(&sphere_count)->default_value(5)
    ->composing(), "Sphere count for grid output. Default: 5.\n")

    ("output_type", PO::value<std::string>(&io_opts._out_type)
    ->default_value("grid_pdb")->composing(), "raw_pdb: null areas as "
    "tetrahedrons (residues). raw_cgo: null areas as tetrahedrons formed with "
    "CGO lines in a pymol script. grid_pdb: null areas filled with points. "
    "grid_cgo: null areas filled with CGO spheres in a pymol script. "
    "Default: grid_pdb\n");

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
    cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    cerr << "Exception of unknown type when reading config file.\n";
  }

  if(io_opts._in_filename == "none" && tool_pdb_to_ch == "none" &&
     tool_pdb_norm == "none"         && tool_aa_to_ca == "none") {
    std::cerr << "Error: option '--input_struct' is required but missing."
    << '\n';
    return 1;
  }

  if (triangulate_only_included_aas == true && AA_indices_proto == "none") {
      std::cerr << "ERROR: 'triangulate_only_included_aas' is TRUE, but "
      "'included_amino_acids' was not set." << "\n";
      return 1;
  }

  if (IA_opts._resn_proto != "none" && IA_opts._atom_proto != "none") {
    IA_opts._atom_proto = "none";
    std::cerr << "Input warning: Both 'included_area_residues' and "
    "'included_area_atoms' were set. Using the former. " << "\n";
  }


  bool const defined_included_area = (IA_opts._resn_proto != "none") || 
      (IA_opts._atom_proto != "none") || (IA_opts._sphere_proto != "none") || 
      (IA_opts._cylinder_proto != "none") || (IA_opts._prism_proto != "none") || 
      (IA_opts._filename != "none");
  bool const modes_or_ndd = (NDD_opts._modes_ndd_filename != "none") || 
      (io_opts._in_md_filename != "none");

  if (!defined_included_area && modes_or_ndd) {
    std::cerr << "Input error: You are running ANA MD/NDD without an inclusion "
    "area. Check ANA's manual." << "\n";
    return 1;
  }

  // if (NDD_opts._modes_format != "amber") {
  //   IA_opts._opt = NDDOptions::amber;
  // } else if (NDD_opts._modes_format != "column") {
  //   IA_opts._opt = NDDOptions::column;
  // } else if (NDD_opts._modes_format != "row") {
  //   IA_opts._opt = NDDOptions::row ;
  // }

  if (NDD_opts._modes_ndd_filename == "none" && NDD_opts._out_ndd_filename == "none") {
    std::cerr << "Input error: NDD_input/NDD_output filename was not set." << "\n";
    return 1;
  }

  if ((NDD_opts._evalues_ndd_filename != "none") && 
  ((NDD_opts._modes_format != "row") && (NDD_opts._modes_format != "column"))) {
    std::cerr << "Input error: NDD_modes_format should be set to row or column." << "\n";
    return 1;
  }

  cell_opts.update();
  
  return 0;
}

} // namespace ANAs
