#ifndef ANA_CONVEX_HULL_H
#define ANA_CONVEX_HULL_H
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>
#include <ANA/SoAPrimitives.hpp>

namespace ANA {

// Refine the provided list of amino acids. Throws a lot.
std::vector<int> adapt_AA_list(std::string &aa_list_proto, int const top = 0);

struct ResidueTag {};
struct AtomTag {};
struct SphereTag {};
struct CylinderTag {};
struct PrismTag {};
struct FileTag {};

struct ConvexHull {
public:
    ConvexHull() = default;

    ConvexHull(
        Molecule const &protein, std::string const &resn_proto, ResidueTag);

    ConvexHull(Molecule const &protein, std::string const &atom_proto, AtomTag);

    ConvexHull(std::string const &sphere_proto, SphereTag);

    ConvexHull(std::string const &cylinder_proto, CylinderTag);

    ConvexHull(std::string const &prism_proto, PrismTag);

    ConvexHull(std::string const &filename, FileTag);

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
        _normals.reserve(triangles);

        auto f_end = CH.facets_end();
        for (auto f_ite = CH.facets_begin(); f_ite != f_end; ++f_ite) {
            // Fix around the weirdest CGAL bug.
            auto he_ite = f_ite->facet_begin();
            auto const p0 = he_ite->vertex()->point();
            he_ite++;
            auto const p1 = he_ite->vertex()->point();
            he_ite++;
            auto const p2 = he_ite->vertex()->point();
            CVector const v1 = p1 - p0;
            CVector const v2 = p2 - p1;
            CVector normal = CGAL::cross_product(v2, v1);
            normal =
                normal / std::sqrt(CGAL::to_double(normal.squared_length()));
            _normals.push_back(normal);

            _triangles.emplace_back(p0, p1, p2);
        }

        return;
    }

    std::vector<CTriangle> _triangles;
    std::vector<CVector> _normals;
    // Used in case the CH was constructed from residues or atoms.
    std::vector<int> _included_resis, _included_atoms;
};

ConvexHull create_convex_hull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts);

} // namespace ANA

#endif // _H