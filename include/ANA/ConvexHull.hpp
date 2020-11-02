#ifndef ANA_CONVEX_HULL_H
#define ANA_CONVEX_HULL_H
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/PrimitivesUtils.hpp>
#include <ANA/SoAPrimitives.hpp>

namespace ANA {

struct ResidueTag { };
struct AtomTag { };
struct SphereTag { };
struct CylinderTag { };
struct PrismTag { };
struct FileTag { };

class ConvexHull {
public:
    ConvexHull() = default;

    ConvexHull(Molecule const &protein, IncludedAreaOptions const &resn_proto,
        ResidueTag);

    ConvexHull(Molecule const &protein, IncludedAreaOptions const &atom_proto,
        AtomTag);

    ConvexHull(IncludedAreaOptions const &sphere_proto, SphereTag);

    ConvexHull(IncludedAreaOptions const &cylinder_proto, CylinderTag);

    ConvexHull(IncludedAreaOptions const &prism_proto, PrismTag);

    ConvexHull(IncludedAreaOptions const &filename, FileTag);

    // Run the actual Convex Hull algorithm with CGAL.
    // Should only work with containers. Will fix w/ c++20. TODO
    template <class T>
    void run_convex_hull(T const &points) {

        if (points.size() < 4) {
            throw std::runtime_error(
                "Not possible to triangulate less than 4 points. Aborting.");
        }

        Polyhedron CH;
        CGAL::convex_hull_3(points.begin(), points.end(), CH);

        auto const triangles = CH.size_of_facets();
        _triangles.reserve(triangles);
        _normal_0.reserve(triangles);
        _normal_1.reserve(triangles);
        _normal_2.reserve(triangles);

        auto f_end = CH.facets_end();
        for (auto f_ite = CH.facets_begin(); f_ite != f_end; ++f_ite) {
            // Fix around the weirdest CGAL bug.
            auto he_ite = f_ite->facet_begin();
            Point const p0(he_ite->vertex()->point());
            he_ite++;
            Point const p1(he_ite->vertex()->point());
            he_ite++;
            Point const p2(he_ite->vertex()->point());
            Vector const v01 = p1 - p0, v02 = p2 - p0;
            Vector const v10 = p0 - p1, v12 = p2 - p1;
            Vector const v20 = p0 - p2, v21 = p1 - p2;

            _normal_0.emplace_back(normalize(cross_product(v01, v02)));
            _normal_1.emplace_back(normalize(cross_product(v10, v12)));
            _normal_2.emplace_back(normalize(cross_product(v20, v21)));

            _v01.push_back(v01);
            _v02.push_back(v02);
            _triangles.emplace_back(p0, p1, p2);
        }

        return;
    }

    // Unfortunately I cant add Info on CGAL's Convex Hull Points, so I have to
    // do this.
    void add_res_info(Molecule const &protein);
    // Unfortunately I cant add Info on CGAL's Convex Hull Points, so I have to
    // do this.
    void add_atm_info(Molecule const &protein);

    // Constructor for NDD. Returns an updated Convex Hull displacing the input
    // convex hull along the input vector scaled by the step_size.
    ConvexHull(ConvexHull const &CH, std::vector<double> const &vector,
        double const step_size);

    std::vector<Triangle> _triangles;
    std::vector<Vector> _normal_0, _normal_1, _normal_2, _v01, _v02;
    bool _dynamic = false;
    // Used in case the CH was constructed from residues or atoms.
    std::vector<int> _included_resis, _included_atoms;
    std::vector<TrianInfo> _info;
};

ConvexHull create_convex_hull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts);

// Refine the provided list of amino acids. Throws a lot.
std::vector<int> string_to_list(
    std::string const &list_proto, int const top = 0);

} // namespace ANA

#endif // _H