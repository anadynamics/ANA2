#include <ANA/PrimitivesUtils.hpp>

namespace ANA {
// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(Point const &p0, Point const &p1, Point const &p2,
    Point const &p3, double const radius) {

    // Get the mini tetrahedron's points.
    Vector v01 = normalize(p1 - p0) * radius;
    Point p01 = p0 + v01;

    Vector v02 = normalize(p2 - p0) * radius;
    Point p02 = p0 + v02;

    Vector v03 = normalize(p3 - p0) * radius;
    Point p03 = p0 + v03;

    double const mini_tetrahedron_vol = volume(p0, p01, p02, p03);

    // Get the distance between p0 and the plane formed by p1, p2 and p3.
    Vector const plane_normal = get_normal(p01, p02, p03);
    double const dist_to_plane = dot_product(plane_normal, v01);
    double const h = radius - dist_to_plane;
    // Now, get the volume of the sphere's slice.
    double const spherical_cap_volume =
        std::abs(M_PI3 * h * h * (3 * radius - h));

    // Add the 2 volumes that represent the space occupied by the atom with
    // coordinates p0.
    double const volume = mini_tetrahedron_vol + spherical_cap_volume;
    if (isnan(volume)) {
        return 0.;
    } else {
        return volume;
    }
}
} // namespace ANA
