#include <ANA/Read.hpp>

namespace ANA {

// Read coordinates in pdb format using chemfiles.
bool read_static(std::string const &filename,
    bool const triangulate_only_included_aas, bool const atom_only,
    std::string &aa_list_proto, std::string &exclude_ca_for_ASA_proto,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    ANA_molecule &molecule_points, CPoint &cm, std::vector<int> &aa_list,
    std::vector<int> &CA_indices, std::vector<CPoint> &CAs_Points,
    std::vector<int> &include_CH_atoms, Triang_Vector &CH_triangs,
    std::vector<int> &hetatm_atoms) {

    bool requested_CH = false;
    std::vector<int> exclude_ca_for_ASA, include_CH_aa;
    std::string elmnt;
    std::vector<CPoint> incl_area_points;
    Finite_vertices_iterator fv_ite;
    // For: // Sort "molecule_points" according to atom indices. Needed for
    // wall_output in MD
    std::vector<int> atom_indices_to_sort;

    // Read PDB
    chemfiles::Trajectory input_pdb_traj(filename);
    chemfiles::Frame input_pdb_frame = input_pdb_traj.read();
    auto in_xyz = input_pdb_frame.positions();
    chemfiles::Topology input_pdb_top = input_pdb_frame.topology();
    // get center of mass
    int const natoms = input_pdb_top.natoms();
    cm = getCM(in_xyz, natoms);
    // Turn list of aminoacids and atoms from strings to vectors of ints.
    bool const listed_included_aa = adapt_AA_list(aa_list_proto, aa_list);
    auto const AA_end = aa_list.end();
    bool const listed_excluded_ASA =
        adapt_AA_list(exclude_ca_for_ASA_proto, exclude_ca_for_ASA);
    auto const CA_ASA_end = exclude_ca_for_ASA.end();
    // Residue selection takes precedence over atom selection, for included
    // area.
    bool const listed_incl_res_CH =
        ANA::adapt_AA_list(include_CH_aa_proto, include_CH_aa);
    auto const CH_AA_end = include_CH_aa.end();
    bool listed_incl_atom_CH = false;
    if (!listed_incl_res_CH) {
        listed_incl_atom_CH =
            ANA::adapt_AA_list(include_CH_atom_proto, include_CH_atoms);
    }
    bool listed_incl_CH = false;

    // If requested, use atoms indices to get included area.
    if (listed_incl_atom_CH) {
        // Check the input atoms don't go out of bounds. This array is sorted.
        if (include_CH_atoms[include_CH_atoms.size() - 1] > natoms) {
            std::cerr
                << "included_area_atoms goes out of bounds. Check this input "
                   "list and your input PDB atom count. Quiting now."
                << '\n';
            exit(0);
        }
        for (auto &each : include_CH_atoms) {
            // 0-index normalization.
            each = each - 1;
            incl_area_points.push_back(
                CPoint(in_xyz[each][0], in_xyz[each][1], in_xyz[each][2]));
            listed_incl_CH = true;
        }
    }
    // Store the atoms points along with information in the variable
    // molecule_points
    if (triangulate_only_included_aas && listed_included_aa) {
        // Only triangulate the included amino acids
        if (aa_list.size() < 4) {
            throw std::runtime_error(
                "Not possible to triangulate less than 4 points.");
        }

        for (auto const &residuo : input_pdb_top.residues()) {
            auto resid = residuo.id().value();
            if (std::lower_bound(aa_list.begin(), AA_end, resid) != AA_end) {
                // The current residue was specified
                auto res_name = residuo.name();
                for (auto const &i : residuo) {
                    if (atom_only &&
                        input_pdb_top[i]
                            .get("is_hetatm")
                            .value_or(false)
                            .as_bool()) {
                        hetatm_atoms.push_back(i);
                        continue;
                    }
                    // Iterate over each atom, construct the vertex info and get
                    // coords
                    // GetVdwRad or GetCovalentRad?
                    VertexInfo vi1;
                    vi1._resi = res_name;
                    vi1._index = i;
                    vi1._resn = resid;
                    // For later sorting
                    atom_indices_to_sort.push_back(i);
                    auto vdw = input_pdb_top[i].vdw_radius();
                    if (vdw) {
                        vi1._radius = vdw.value();
                    } else {
                        vi1._radius = 1.5;
                        std::cerr << "Element from atom " << i + 1
                                  << " not available. "
                                  << " Using covalent radius of 1.5." << '\n';
                    }
                    CPoint p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
                    molecule_points.push_back(std::make_pair(p1, vi1));

                    if (input_pdb_top[i].name() == "CA") {
                        // If its a Calpha atom store its coordinates.
                        if (listed_excluded_ASA == false ||
                            std::lower_bound(exclude_ca_for_ASA.begin(),
                                CA_ASA_end, resid) == CA_ASA_end) {
                            CAs_Points.push_back(p1);
                            CA_indices.push_back(i);
                        }
                        // Use it to draw the included convex hull, if requested
                        if (listed_incl_res_CH) {
                            if (std::binary_search(
                                    include_CH_aa.begin(), CH_AA_end, resid)) {
                                // The current residue was specified
                                incl_area_points.push_back(p1);
                                include_CH_atoms.push_back(i);
                            }
                        }
                    }
                }
            } else {
                // The current residue was not specified
                for (auto const &i : residuo) {
                    if (atom_only &&
                        input_pdb_top[i]
                            .get("is_hetatm")
                            .value_or(false)
                            .as_bool()) {
                        hetatm_atoms.push_back(i);
                        continue;
                    }
                    if (input_pdb_top[i].name() == "CA") {
                        auto resid = residuo.id().value();
                        CPoint p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
                        if (listed_excluded_ASA == false ||
                            std::lower_bound(exclude_ca_for_ASA.begin(),
                                CA_ASA_end, resid) == CA_ASA_end) {
                            // If its a Calpha atom and wasn't excluded store
                            // its
                            // coordinates.
                            CAs_Points.push_back(p1);
                            CA_indices.push_back(i);
                        }
                        // Use it to draw the included convex hull, if requested
                        if (listed_incl_res_CH) {

                            if (std::binary_search(
                                    include_CH_aa.begin(), CH_AA_end, resid)) {
                                // The current residue was specified
                                incl_area_points.push_back(p1);
                                include_CH_atoms.push_back(i);
                            }
                        }
                    }
                }
            }
        }
    } else { // triangulate the whole molecule
        for (auto const &residuo : input_pdb_top.residues()) {
            auto res_name = residuo.name();
            for (auto const &i : residuo) {
                if (atom_only &&
                    input_pdb_top[i]
                        .get("is_hetatm")
                        .value_or(false)
                        .as_bool()) {
                    // Save the HEATM indices to discard them during MD and NDD
                    // runs.
                    hetatm_atoms.push_back(i);
                    continue;
                }
                // Iterate over each atom, construct the vertex info and get
                // coords
                // GetVdwRad or GetCovalentRad?
                VertexInfo vi1;
                vi1._resi = res_name;
                vi1._index = i;
                // For later sorting
                atom_indices_to_sort.push_back(i);
                auto vdw = input_pdb_top[i].vdw_radius();
                if (vdw) {
                    vi1._radius = vdw.value();
                } else {
                    vi1._radius = .5;
                    std::cerr << "Element from atom " << i + 1
                              << " not available. "
                              << "Using Van Der Walls radii of 1.5." << '\n';
                }
                auto resid = residuo.id().value();
                vi1._resi = resid;
                CPoint p1(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
                molecule_points.push_back(std::make_pair(p1, vi1));

                if (input_pdb_top[i].name() == "CA") {
                    if (listed_excluded_ASA == false ||
                        std::lower_bound(exclude_ca_for_ASA.begin(), CA_ASA_end,
                            resid) == CA_ASA_end) {
                        // If it wasn't excluded, store its coordinates.
                        CAs_Points.push_back(p1);
                        CA_indices.push_back(i);
                    }
                    // Use it to draw the included convex hull, if requested
                    if (listed_incl_res_CH) {
                        if (std::binary_search(
                                include_CH_aa.begin(), CH_AA_end, resid)) {
                            // The current residue was specified
                            incl_area_points.push_back(p1);
                            include_CH_atoms.push_back(i);
                        }
                    }
                }
            }
        }
    }

    // Get convex hull if requested.
    if (!(listed_incl_CH || listed_incl_res_CH) &&
        include_CH_filename != "none") {
        ANA::read_included_area(include_CH_filename, incl_area_points);
        if (incl_area_points.size() < 4) {
            throw std::runtime_error(
                "Not possible to triangulate less than 4 points.");
        }

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }
        requested_CH = true;

    } else if (!(listed_incl_CH || listed_incl_res_CH) &&
        sphere_proto != "none") {
        std::stringstream stream_sphere(sphere_proto);
        double x = parse_double(stream_sphere);
        double y = parse_double(stream_sphere);
        double z = parse_double(stream_sphere);
        double r = parse_double(stream_sphere);
        double const cos_30 = sqrt(3) / 2;
        double const sin_30 = 0.5;

        CPoint center(x, y, z);
        incl_area_points.push_back(center + CVector(r, 0, 0));
        incl_area_points.push_back(center + CVector(0, r, 0));
        incl_area_points.push_back(center + CVector(0, 0, r));
        incl_area_points.push_back(center + CVector(-r, 0, 0));
        incl_area_points.push_back(center + CVector(0, -r, 0));
        incl_area_points.push_back(center + CVector(0, 0, -r));
        // X-Y plane
        incl_area_points.push_back(center + CVector(r * cos_30, r * sin_30, 0));
        incl_area_points.push_back(center + CVector(r * sin_30, r * cos_30, 0));
        incl_area_points.push_back(
            center + CVector(r * cos_30, -r * sin_30, 0));
        incl_area_points.push_back(
            center + CVector(r * sin_30, -r * cos_30, 0));
        incl_area_points.push_back(
            center + CVector(-r * cos_30, r * sin_30, 0));
        incl_area_points.push_back(
            center + CVector(-r * sin_30, r * cos_30, 0));
        incl_area_points.push_back(
            center + CVector(-r * cos_30, -r * sin_30, 0));
        incl_area_points.push_back(
            center + CVector(-r * sin_30, -r * cos_30, 0));
        // X-Z plane
        incl_area_points.push_back(center + CVector(r * cos_30, 0, r * sin_30));
        incl_area_points.push_back(center + CVector(r * sin_30, 0, r * cos_30));
        incl_area_points.push_back(
            center + CVector(r * cos_30, 0, -r * sin_30));
        incl_area_points.push_back(
            center + CVector(r * sin_30, 0, -r * cos_30));
        incl_area_points.push_back(
            center + CVector(-r * cos_30, 0, r * sin_30));
        incl_area_points.push_back(
            center + CVector(-r * sin_30, 0, r * cos_30));
        incl_area_points.push_back(
            center + CVector(-r * cos_30, 0, -r * sin_30));
        incl_area_points.push_back(
            center + CVector(-r * sin_30, 0, -r * cos_30));
        // Y-Z plane
        incl_area_points.push_back(center + CVector(0, r * cos_30, r * sin_30));
        incl_area_points.push_back(center + CVector(0, r * sin_30, r * cos_30));
        incl_area_points.push_back(
            center + CVector(0, r * cos_30, -r * sin_30));
        incl_area_points.push_back(
            center + CVector(0, r * sin_30, -r * cos_30));
        incl_area_points.push_back(
            center + CVector(0, -r * cos_30, r * sin_30));
        incl_area_points.push_back(
            center + CVector(0, -r * sin_30, r * cos_30));
        incl_area_points.push_back(
            center + CVector(0, -r * cos_30, -r * sin_30));
        incl_area_points.push_back(
            center + CVector(0, -r * sin_30, -r * cos_30));

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }
        requested_CH = true;

    } else if (!(listed_incl_CH || listed_incl_res_CH) &&
        cylinder_proto != "none") {
        std::stringstream stream_cylinder(cylinder_proto);

        double x1 = parse_double(stream_cylinder);
        double y1 = parse_double(stream_cylinder);
        double z1 = parse_double(stream_cylinder);
        double x2 = parse_double(stream_cylinder);
        double y2 = parse_double(stream_cylinder);
        double z2 = parse_double(stream_cylinder);
        double r = parse_double(stream_cylinder);

        double const cos_30 = sqrt(3) / 2;
        double const sin_30 = 0.5;

        CPoint center_1(x1, y1, z1);
        CPoint center_2(x2, y2, z2);
        CVector vdiff(center_2 - center_1);
        CVector n1(-vdiff.y(), vdiff.x(), 0);
        CVector n2 = CGAL::cross_product(vdiff, n1);
        n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
        n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

        // tap 1
        incl_area_points.push_back(center_1 + r * n1);
        incl_area_points.push_back(center_1 + r * n2);
        incl_area_points.push_back(center_1 - r * n1);
        incl_area_points.push_back(center_1 - r * n2);
        incl_area_points.push_back(
            center_1 + r * cos_30 * n1 + r * sin_30 * n2);
        incl_area_points.push_back(
            center_1 + r * sin_30 * n1 + r * cos_30 * n2);
        incl_area_points.push_back(
            center_1 + r * cos_30 * n1 - r * sin_30 * n2);
        incl_area_points.push_back(
            center_1 + r * sin_30 * n1 - r * cos_30 * n2);
        incl_area_points.push_back(
            center_1 - r * cos_30 * n1 + r * sin_30 * n2);
        incl_area_points.push_back(
            center_1 - r * sin_30 * n1 + r * cos_30 * n2);
        incl_area_points.push_back(
            center_1 - r * cos_30 * n1 - r * sin_30 * n2);
        incl_area_points.push_back(
            center_1 - r * sin_30 * n1 - r * cos_30 * n2);

        // tap 2
        incl_area_points.push_back(center_2 + r * n1);
        incl_area_points.push_back(center_2 + r * n2);
        incl_area_points.push_back(center_2 - r * n1);
        incl_area_points.push_back(center_2 - r * n2);
        incl_area_points.push_back(
            center_2 + r * cos_30 * n1 + r * sin_30 * n2);
        incl_area_points.push_back(
            center_2 + r * sin_30 * n1 + r * cos_30 * n2);
        incl_area_points.push_back(
            center_2 + r * cos_30 * n1 - r * sin_30 * n2);
        incl_area_points.push_back(
            center_2 + r * sin_30 * n1 - r * cos_30 * n2);
        incl_area_points.push_back(
            center_2 - r * cos_30 * n1 + r * sin_30 * n2);
        incl_area_points.push_back(
            center_2 - r * sin_30 * n1 + r * cos_30 * n2);
        incl_area_points.push_back(
            center_2 - r * cos_30 * n1 - r * sin_30 * n2);
        incl_area_points.push_back(
            center_2 - r * sin_30 * n1 - r * cos_30 * n2);

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }
        requested_CH = true;

    } else if (!(listed_incl_CH || listed_incl_res_CH) &&
        prism_proto != "none") {
        std::stringstream stream_prism(prism_proto);
        double x1 = parse_double(stream_prism);
        double y1 = parse_double(stream_prism);
        double z1 = parse_double(stream_prism);
        double x2 = parse_double(stream_prism);
        double y2 = parse_double(stream_prism);
        double z2 = parse_double(stream_prism);
        double width = parse_double(stream_prism) / 2;
        double height = parse_double(stream_prism) / 2;

        CPoint center_1(x1, y1, z1);
        CPoint center_2(x2, y2, z2);
        CVector vdiff(center_2 - center_1);
        CVector n1(-vdiff.y(), vdiff.x(), 0);
        CVector n2 = CGAL::cross_product(vdiff, n1);
        n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
        n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

        // tap 1
        incl_area_points.push_back(center_1 + width * n1 + height * n2);
        incl_area_points.push_back(center_1 + width * n1 - height * n2);
        incl_area_points.push_back(center_1 - width * n1 + height * n2);
        incl_area_points.push_back(center_1 - width * n1 - height * n2);

        // tap 2
        incl_area_points.push_back(center_2 + width * n1 + height * n2);
        incl_area_points.push_back(center_2 + width * n1 - height * n2);
        incl_area_points.push_back(center_2 - width * n1 + height * n2);
        incl_area_points.push_back(center_2 - width * n1 - height * n2);

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }
        requested_CH = true;

    } else if (listed_incl_CH || listed_incl_res_CH) {
        if (incl_area_points.size() < 4) {
            throw std::runtime_error(
                "Not possible to triangulate less than 4 points.");
        }

        Polyhedron CH;
        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {

            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }
        requested_CH = true;
    }
    return requested_CH;
}

// Read coordinates in netcdf format.
void read_MD(const chemfiles::Frame &in_frame, bool const requested_CH,
    std::string const &sphere_proto, std::string const &cylinder_proto,
    std::string const &prism_proto, const std::vector<int> &hetatm_atoms,
    std::vector<int> &include_CH_atoms, std::string const &include_CH_filename,
    Triang_Vector &CH_triangs, std::string const &ASA_method,
    const std::vector<int> &CA_indices, std::vector<CPoint> &CAs_points,
    ANA_molecule &molecule_points) {
    int i;
    CPoint p1;

    // Get positions. Discard hetatms if so requested.
    int const natoms = molecule_points.size();
    auto in_xyz = in_frame.positions();
    if (hetatm_atoms.size() != 0) {
        auto ite_beg = hetatm_atoms.begin();
        auto ite_end = hetatm_atoms.end();
        for (i = 0; i < natoms; ++i) {
            if (std::binary_search(ite_beg, ite_end, i)) {
                // This is an HETATM
                continue;
            }
            molecule_points[i].first =
                CPoint(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
        }
    } else {
        for (i = 0; i < natoms; ++i) {
            molecule_points[i].first =
                CPoint(in_xyz[i][0], in_xyz[i][1], in_xyz[i][2]);
        }
    }

    if (requested_CH && include_CH_filename == "none" &&
        sphere_proto == "none" && cylinder_proto == "none" &&
        prism_proto == "none") {
        // Update the convex hull of the included area
        CH_triangs.clear();
        Polyhedron CH;
        std::vector<CPoint> incl_area_points;

        for (auto const &idx : include_CH_atoms) {
            incl_area_points.push_back(
                CPoint(in_xyz[idx][0], in_xyz[idx][1], in_xyz[idx][2]));
        }

        CGAL::convex_hull_3(
            incl_area_points.begin(), incl_area_points.end(), CH);
        P_Facet_const_iterator f_end = CH.facets_end();
        for (P_Facet_const_iterator f_ite = CH.facets_begin(); f_ite != f_end;
             ++f_ite) {
            // Fix around the weirdest CGAL bug.
            P_Halfedge_around_facet_const_circulator he_ite =
                f_ite->facet_begin();
            auto const he_ite_0 = he_ite++;
            auto const he_ite_1 = he_ite++;
            auto const he_ite_2 = he_ite;

            CH_triangs.push_back(CTriangle(he_ite_0->vertex()->point(),
                he_ite_1->vertex()->point(), he_ite_2->vertex()->point()));
        }

    } else {
        // Update Calphas positions if ASA method was specified instead of
        // included area
        if (ASA_method == "dot_pdt") {
            // Get CAs positions
            i = 0;
            for (auto const &CA_index : CA_indices) {
                CAs_points[i] = CPoint(in_xyz[CA_index][0], in_xyz[CA_index][1],
                    in_xyz[CA_index][2]);
                ++i;
            }
        }
    }

    return;
}

// Get the center of mass from a chemfiles molecule.
inline CPoint getCM(
    const chemfiles::span<chemfiles::Vector3D> &in_xyz, int const natoms) {
    double x_cm = 0;
    double y_cm = 0;
    double z_cm = 0;

    for (int i = 0; i < natoms; ++i) {
        x_cm += in_xyz[i][0];
        y_cm += in_xyz[i][1];
        z_cm += in_xyz[i][2];
    }
    return CPoint(x_cm / natoms, y_cm = y_cm / natoms, z_cm = z_cm / natoms);
}
// Read file with included area coordinates
inline void read_included_area(
    std::string const &filename, std::vector<CPoint> &area_points) {

    std::string linea, xs, ys, zs;
    std::ifstream infile(filename);

    if (infile.eof()) {
        throw std::runtime_error("Not enough parameters for "
                                 "pseudosphere/pseudocylinder/ prism input. "
                                 "Aborting.");
    }

    if (infile.is_open()) {
        while (std::getline(infile, linea)) {
            if (linea.substr(0, 1) == "#") {
                // a comment
                continue;
            }
            std::stringstream streamm(linea);
            double x, y, z;

            try {
                streamm >> x >> y >> z;
                area_points.push_back(CPoint(x, y, z));

            } catch (const std::invalid_argument &ia) {
                // some character present
                throw std::invalid_argument(
                    "Invalid input: " + xs + "  " + ys + "  " + zs);
                continue;

            } catch (...) {
                // some other exception.
                throw std::runtime_error(
                    "Invalid input: " + xs + "  " + ys + "  " + zs);
                continue;
            }
        }
    } else
        throw std::runtime_error("Failed to read included area");

    return;
}

} // namespace ANA
