#include <ANA/ConvexHull.hpp>

namespace ANA {

ConvexHull create_convex_hull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts) {

    if (IA_opts._resn_proto != "none") {
        return ConvexHull(protein, IA_opts, ResidueTag());
    } else if (IA_opts._atom_proto != "none") {
        return ConvexHull(protein, IA_opts, AtomTag());
    } else if (IA_opts._sphere_proto != "none") {
        return ConvexHull(IA_opts, SphereTag());
    } else if (IA_opts._cylinder_proto != "none") {
        return ConvexHull(IA_opts, CylinderTag());
    } else if (IA_opts._prism_proto != "none") {
        return ConvexHull(IA_opts, PrismTag());
    } else if (IA_opts._filename != "none") {
        return ConvexHull(IA_opts, FileTag());
    } else {
        throw std::invalid_argument(
            "No Convex Hull input could be parsed. This shouldn't happen.");
    }

    return {};
}

ConvexHull::ConvexHull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts, ResidueTag) {

    _included_resis = string_to_list(IA_opts._resn_proto, protein._nres);

    std::vector<CPoint> incl_area_points;
    incl_area_points.reserve(_included_resis.size());

    for (auto const i : _included_resis) {
        incl_area_points.push_back(
            protein._data[protein._alphaCarbons[i - 1]].first);
    }

    try {
        run_convex_hull(incl_area_points);
        if (IA_opts.has_info) {
            add_res_info(protein);
        }
        _dynamic = true;
    } catch (std::runtime_error const &e) {
        throw std::runtime_error(
            "runtime_error error when triangulating convex hull. Aborting.");
    } catch (...) {
        throw std::runtime_error(
            "Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts, AtomTag) {

    _included_atoms = string_to_list(IA_opts._atom_proto, protein._natoms);

    std::vector<CPoint> incl_area_points;
    incl_area_points.reserve(_included_atoms.size());

    for (auto const i : _included_atoms) {
        // 0-index normalization.
        incl_area_points.push_back(protein._data[i - 1].first);
    }

    try {
        run_convex_hull(incl_area_points);
        if (IA_opts.has_info) {
            add_atm_info(protein);
        }
        _dynamic = true;
    } catch (std::runtime_error const &e) {
        throw std::runtime_error(
            "runtime_error error when triangulating convex hull. Aborting.");
    } catch (...) {
        throw std::runtime_error(
            "Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(IncludedAreaOptions const &IA_opts, SphereTag) {

    std::stringstream stream_sphere(IA_opts._sphere_proto);
    double const x = parse_double(stream_sphere);
    double const y = parse_double(stream_sphere);
    double const z = parse_double(stream_sphere);
    double const r = parse_double(stream_sphere);
    double constexpr cos_30 = 0.86602540378;
    double constexpr sin_30 = 0.5;

    CPoint const center(x, y, z);
    // Drawin a pseudo-sphere. The first 6 correspondon the XYZ axes, the
    // next 8 to the X-Y plane, then 8 more for the X-Z plane and 8 for the
    // Y-Z plane.
    std::array<CPoint, 30> const incl_area_points{center + CVector(r, 0, 0),
        center + CVector(0, r, 0), center + CVector(0, 0, r),
        center + CVector(-r, 0, 0), center + CVector(0, -r, 0),
        center + CVector(0, 0, -r), center + CVector(r * cos_30, r * sin_30, 0),
        center + CVector(r * sin_30, r * cos_30, 0),
        center + CVector(r * cos_30, -r * sin_30, 0),
        center + CVector(r * sin_30, -r * cos_30, 0),
        center + CVector(-r * cos_30, r * sin_30, 0),
        center + CVector(-r * sin_30, r * cos_30, 0),
        center + CVector(-r * cos_30, -r * sin_30, 0),
        center + CVector(-r * sin_30, -r * cos_30, 0),
        center + CVector(r * cos_30, 0, r * sin_30),
        center + CVector(r * sin_30, 0, r * cos_30),
        center + CVector(r * cos_30, 0, -r * sin_30),
        center + CVector(r * sin_30, 0, -r * cos_30),
        center + CVector(-r * cos_30, 0, r * sin_30),
        center + CVector(-r * sin_30, 0, r * cos_30),
        center + CVector(-r * cos_30, 0, -r * sin_30),
        center + CVector(-r * sin_30, 0, -r * cos_30),
        center + CVector(0, r * cos_30, r * sin_30),
        center + CVector(0, r * sin_30, r * cos_30),
        center + CVector(0, r * cos_30, -r * sin_30),
        center + CVector(0, r * sin_30, -r * cos_30),
        center + CVector(0, -r * cos_30, r * sin_30),
        center + CVector(0, -r * sin_30, r * cos_30),
        center + CVector(0, -r * cos_30, -r * sin_30),
        center + CVector(0, -r * sin_30, -r * cos_30)};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(IncludedAreaOptions const &IA_opts, CylinderTag) {

    std::stringstream stream_cylinder(IA_opts._cylinder_proto);

    double const x1 = parse_double(stream_cylinder);
    double const y1 = parse_double(stream_cylinder);
    double const z1 = parse_double(stream_cylinder);
    double const x2 = parse_double(stream_cylinder);
    double const y2 = parse_double(stream_cylinder);
    double const z2 = parse_double(stream_cylinder);
    double const r = parse_double(stream_cylinder);
    double constexpr cos_30 = 0.86602540378;
    double const sin_30 = 0.5;

    CPoint const center_1(x1, y1, z1), center_2(x2, y2, z2);

    CVector const vdiff(center_2 - center_1);
    CVector n1(-vdiff.y(), vdiff.x(), 0);
    CVector n2 = CGAL::cross_product(vdiff, n1);
    n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
    n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

    // The first 12 correspond to the first tap, the other half correspond
    // to the 2nd tap.
    std::array<CPoint, 24> incl_area_points{center_1 + r * n1,
        center_1 + r * n2, center_1 - r * n1, center_1 - r * n2,
        center_1 + r * cos_30 * n1 + r * sin_30 * n2,
        center_1 + r * sin_30 * n1 + r * cos_30 * n2,
        center_1 + r * cos_30 * n1 - r * sin_30 * n2,
        center_1 + r * sin_30 * n1 - r * cos_30 * n2,
        center_1 - r * cos_30 * n1 + r * sin_30 * n2,
        center_1 - r * sin_30 * n1 + r * cos_30 * n2,
        center_1 - r * cos_30 * n1 - r * sin_30 * n2,
        center_1 - r * sin_30 * n1 - r * cos_30 * n2, center_2 + r * n1,
        center_2 + r * n2, center_2 - r * n1, center_2 - r * n2,
        center_2 + r * cos_30 * n1 + r * sin_30 * n2,
        center_2 + r * sin_30 * n1 + r * cos_30 * n2,
        center_2 + r * cos_30 * n1 - r * sin_30 * n2,
        center_2 + r * sin_30 * n1 - r * cos_30 * n2,
        center_2 - r * cos_30 * n1 + r * sin_30 * n2,
        center_2 - r * sin_30 * n1 + r * cos_30 * n2,
        center_2 - r * cos_30 * n1 - r * sin_30 * n2,
        center_2 - r * sin_30 * n1 - r * cos_30 * n2};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(IncludedAreaOptions const &IA_opts, PrismTag) {

    std::stringstream stream_prism(IA_opts._prism_proto);
    double const x1 = parse_double(stream_prism);
    double const y1 = parse_double(stream_prism);
    double const z1 = parse_double(stream_prism);
    double const x2 = parse_double(stream_prism);
    double const y2 = parse_double(stream_prism);
    double const z2 = parse_double(stream_prism);
    double const width = parse_double(stream_prism) / 2;
    double const height = parse_double(stream_prism) / 2;

    CPoint const center_1(x1, y1, z1), center_2(x2, y2, z2);
    CVector const vdiff(center_2 - center_1);
    CVector n1(-vdiff.y(), vdiff.x(), 0);
    CVector n2 = CGAL::cross_product(vdiff, n1);
    n1 = n1 / std::sqrt(CGAL::to_double(n1.squared_length()));
    n2 = n2 / std::sqrt(CGAL::to_double(n2.squared_length()));

    // 8 vertices of a prism.
    std::array<CPoint, 8> incl_area_points{center_1 + width * n1 + height * n2,
        center_1 + width * n1 - height * n2,
        center_1 - width * n1 + height * n2,
        center_1 - width * n1 - height * n2,
        center_2 + width * n1 + height * n2,
        center_2 + width * n1 - height * n2,
        center_2 - width * n1 + height * n2,
        center_2 - width * n1 - height * n2};

    try {
        run_convex_hull(incl_area_points);
    } catch (std::runtime_error const &e) {
        throw;
    } catch (...) {
        throw("Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(IncludedAreaOptions const &IA_opts, FileTag) {
    IA_opts._filename + "ss"; // just to avoid unused value warning.
    throw std::runtime_error(
        "Convex hull construction from input file not supported yet. "
        "Aborting.");
}

// Unfortunately I cant add Info on CGAL's Convex Hull Points, so I have to do
// this.
void ConvexHull::add_res_info(Molecule const &protein) {

    for (auto const &t : _triangles) {
        std::array<int, 3> indices{-666, -666, -666};
        for (auto const i : _included_resis) {
            for (int j = 0; j < 3; ++j) {
                if (equal(protein._data[protein._alphaCarbons[i - 1]].first,
                        t[j])) {
                    if (indices[j] != -666) {
                        throw std::runtime_error(
                            "add_res_info() error. A convex hull's atom "
                            "matched twice.");
                    } else {
                        indices[j] = i - 1;
                    }
                }
            }
            if (indices[0] != -666 && indices[1] != -666 &&
                indices[2] != -666) {
                break;
            }
        }
        _info.emplace_back(
            protein._data[protein._alphaCarbons[indices[0]]].second,
            protein._data[protein._alphaCarbons[indices[1]]].second,
            protein._data[protein._alphaCarbons[indices[2]]].second);
    }

    return;
}

// Unfortunately I cant add Info on CGAL's Convex Hull Points, so I have to do
// this.
void ConvexHull::add_atm_info(Molecule const &protein) {

    for (auto const &t : _triangles) {
        std::array<int, 3> indices{-666, -666, -666};
        for (auto const i : _included_resis) {
            for (int j = 0; j < 3; ++j) {
                if (equal(protein._data[i - 1].first, t[j])) {
                    if (indices[j] != -666) {
                        throw std::runtime_error(
                            "add_res_info() error. A convex hull's atom "
                            "matched twice.");
                    } else {
                        indices[j] = i - 1;
                    }
                }
            }
            if (indices[0] != -666 && indices[1] != -666 &&
                indices[2] != -666) {
                break;
            }
        }
        _info.emplace_back(
            protein._data[protein._alphaCarbons[indices[0]]].second,
            protein._data[protein._alphaCarbons[indices[1]]].second,
            protein._data[protein._alphaCarbons[indices[2]]].second);
    }

    return;
}

// Returns an updated Convex Hull displacing the input convex hull along the
// input vector scaled by the step_size. The vector must be alfa carbon mode.
ConvexHull::ConvexHull(ConvexHull const &CH, std::vector<double> const &evector,
    double const step_size) {

    _normals.reserve(CH._normals.size());
    _triangles.reserve(CH._triangles.size());

    for (size_t t = 0; t < CH._info.size(); ++t) {
        // _resn is 1-indexed.
        int const resi_0_x = (CH._info[t]._resn[0] - 1) * 3;
        int const resi_1_x = (CH._info[t]._resn[1] - 1) * 3;
        int const resi_2_x = (CH._info[t]._resn[2] - 1) * 3;

        Point const p0{CH._triangles[t][0] +
            step_size *
                Vector(evector[resi_0_x], evector[resi_0_x + 1],
                    evector[resi_0_x + 2])};
        Point const p1{CH._triangles[t][1] +
            step_size *
                Vector(evector[resi_1_x], evector[resi_1_x + 1],
                    evector[resi_1_x + 2])};
        Point const p2{CH._triangles[t][2] +
            step_size *
                Vector(evector[resi_2_x], evector[resi_2_x + 1],
                    evector[resi_2_x + 2])};

        Vector const v01 = p1 - p0, v02 = p2 - p0;

        _normals.emplace_back(normalize(cross_product(v01, v02)));
        _v01.push_back(v01);
        _v02.push_back(v02);
        _triangles.emplace_back(p0, p1, p2);
    }

    return;
}

// Refine the provided list of amino acids. Throws a lot.
std::vector<int> string_to_list(std::string const &list_proto, int const top) {

    std::vector<int> list;

    std::stringstream stream_aa(list_proto);
    std::string temp_aa;
    while (!stream_aa.eof()) {
        int aa;
        stream_aa >> temp_aa;
        try {
            aa = std::stoi(temp_aa);
        } catch (std::out_of_range const &oor) {
            // int is too large to be represented by int
            throw std::out_of_range(
                "Invalid atom / residue number in config file. Aborting.");
        } catch (...) {
            // some other exception.
            throw std::runtime_error(
                "Invalid atom / residue number. Aborting.");
        }
        list.push_back(aa);
    }
    // sort list of included amino acids
    std::sort(list.begin(), list.end());

    if (top != 0) {
        if (list[list.size() - 1] > top) {
            throw std::runtime_error(
                "Atom / residue list goes out of bounds. Check this input list "
                "and your input PDB atom / residue count. Aborting.");
        }
    }

    return list;
}

} // namespace ANA
