#include <ANA/CGALUtils.hpp>

namespace ANA {

// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(CPoint const &p0, CPoint const &p1, CPoint const &p2,
    CPoint const &p3, double const radius) {

    // Get 1st point of the mini tetrahedron.
    CVector vec_1 = p1 - p0;
    vec_1 =
        (vec_1 / std::sqrt(CGAL::to_double(vec_1.squared_length()))) * radius;
    CPoint point_1 = p0 + vec_1;
    // Get 2nd point of the mini tetrahedron.
    CVector vec_2 = p2 - p0;
    vec_2 =
        (vec_2 / std::sqrt(CGAL::to_double(vec_2.squared_length()))) * radius;
    CPoint point_2 = p0 + vec_2;
    // Get 3rd point of the mini tetrahedron.
    CVector vec_3 = p3 - p0;
    vec_3 =
        (vec_3 / std::sqrt(CGAL::to_double(vec_3.squared_length()))) * radius;
    CPoint point_3 = p0 + vec_3;
    double const mini_tetrahedron_vol = volume(p0, point_1, point_2, point_3);

    // Get the distance between p0 and the plane formed by p1, p2 and p3.
    auto const plane_normal = normal(p1, p2, p3);
    double const dist_to_plane = CGAL::to_double(plane_normal * vec_1);
    double const h = radius - dist_to_plane;
    // Now, get the volume of the sphere's slice.
    double const spherical_cap_volume =
        std::abs(M_PI3 * h * h * (3 * radius - h));

    // Add the 2 volumes that represent the space occupied by the atom with
    // coordinates p0.
    double const volume = mini_tetrahedron_vol + spherical_cap_volume;
    if (isnan(volume)) {
        return 0;
    } else {
        return volume;
    }
}

} // namespace ANA