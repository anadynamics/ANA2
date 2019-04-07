#include <ANA/ConvexHull.hpp>

namespace ANA {

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

ConvexHull create_convex_hull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts) {

    switch (IA_opts._opt) {
    case IncludedAreaOptions::IAOption::residue:
        return ConvexHull(protein, IA_opts._resn_proto, ResidueTag());
    case IncludedAreaOptions::IAOption::atom:
        return ConvexHull(protein, IA_opts._atom_proto, AtomTag());
    case IncludedAreaOptions::IAOption::sphere:
        return ConvexHull(IA_opts._sphere_proto, SphereTag());
    case IncludedAreaOptions::IAOption::cylinder:
        return ConvexHull(IA_opts._cylinder_proto, CylinderTag());
    case IncludedAreaOptions::IAOption::prism:
        return ConvexHull(IA_opts._prism_proto, PrismTag());
    case IncludedAreaOptions::IAOption::file:
        return ConvexHull(IA_opts._filename, FileTag());
    case IncludedAreaOptions::IAOption::none:
        throw(std::invalid_argument(
            "No Convex Hull input could be parsed. This shouldn't happen."));
        break;
    }
    return {};
}

ConvexHull::ConvexHull(
    Molecule const &protein, std::string const &resn_proto, ResidueTag) {

    _included_resis = string_to_list(resn_proto, protein._nres);

    std::vector<CPoint> incl_area_points;
    incl_area_points.reserve(_included_resis.size());

    for (auto const i : _included_resis) {
        incl_area_points.push_back(
            protein._data[protein._alphaCarbons[i - 1]].first);
    }

    try {
        run_convex_hull(incl_area_points);
        add_res_info(protein);
    } catch (std::runtime_error const &e) {
        throw std::runtime_error(
            "runtime_error error when triangulating convex hull. Aborting.");
    } catch (...) {
        throw std::runtime_error(
            "Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(
    Molecule const &protein, std::string const &atom_proto, AtomTag) {

    _included_atoms = string_to_list(atom_proto, protein._natoms);

    std::vector<CPoint> incl_area_points;
    incl_area_points.reserve(_included_atoms.size());

    for (auto const i : _included_atoms) {
        // 0-index normalization.
        incl_area_points.push_back(protein._data[i - 1].first);
    }

    try {
        run_convex_hull(incl_area_points);
        add_atm_info(protein);
    } catch (std::runtime_error const &e) {
        throw std::runtime_error(
            "runtime_error error when triangulating convex hull. Aborting.");
    } catch (...) {
        throw std::runtime_error(
            "Uknown error when triangulating convex hull. Aborting.");
    }
}

ConvexHull::ConvexHull(std::string const &sphere_proto, SphereTag) {

    std::stringstream stream_sphere(sphere_proto);
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

ConvexHull::ConvexHull(std::string const &cylinder_proto, CylinderTag) {

    std::stringstream stream_cylinder(cylinder_proto);

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

ConvexHull::ConvexHull(std::string const &prism_proto, PrismTag) {

    std::stringstream stream_prism(prism_proto);
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

ConvexHull::ConvexHull(std::string const &filename, FileTag) {
    filename.end(); // just to avoid unused value warning.
    throw("Convex hull construction from input file not supported yet. "
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

} // namespace ANA
