#ifndef ANA_CONVEX_HULL_H
#define ANA_CONVEX_HULL_H
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/Options.hpp>
#include <ANA/PrimitivesUtils.hpp>
#include <ANA/SoAPrimitives.hpp>

namespace ANA {

// Refine the provided list of amino acids. Throws a lot.
std::vector<int> adapt_AA_list(std::string &aa_list_proto, int const top = 0);

// Helper function. Returns all 'in_vec' elements not matching the query input.
inline std::vector<int> get_i_not_equal(
    std::vector<int> const &in_vec, const int query) {
    std::vector<int> out_vec;
    for (auto const &each : in_vec) {
        if (each != query) {
            out_vec.push_back(each);
        }
    }
    return out_vec;
}

struct ResidueTag {};
struct AtomTag {};
struct SphereTag {};
struct CylinderTag {};
struct PrismTag {};
struct FileTag {};

class ConvexHull {
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
            Point const p0(he_ite->vertex()->point());
            he_ite++;
            Point const p1(he_ite->vertex()->point());
            he_ite++;
            Point const p2(he_ite->vertex()->point());
            Vector const v01 = p1 - p0, v02 = p2 - p0;

            _normals.emplace_back(normalize(cross_product(v01, v02)));
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

    std::vector<Triangle> _triangles;
    std::vector<TrianInfo> _info;
    std::vector<Vector> _normals, _v01, _v02;
    // Used in case the CH was constructed from residues or atoms.
    std::vector<int> _included_resis, _included_atoms;
};

ConvexHull create_convex_hull(
    Molecule const &protein, IncludedAreaOptions const &IA_opts);

} // namespace ANA

#endif // _H