#ifndef ANA_PRIMITIVE_UTILS_H
#define ANA_PRIMITIVE_UTILS_H
#include <ANA/Primitives.hpp>

namespace ANA {

inline bool deq(double x, double y) {
    double const z = x - y;
    return (z > zero_bot and z < zero_top);
}

inline double norm(Vector const &v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline Vector normalize(Vector const &v) { return v / norm(v); }

inline double dot_product(Vector const &v, Vector const &w) {
    return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]);
}

inline Vector cross_product(Vector const &v, Vector const &w) {
    return {v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2],
        v[0] * w[1] - v[1] * w[0]};
}

inline double determinant(
    Vector const &v0, Vector const &v1, Vector const &v2) {
    // First, compute the det2x2.
    double const m01 = v0[0] * v1[1] - v0[1] * v1[0];
    double const m02 = v0[0] * v2[1] - v0[1] * v2[0];
    double const m12 = v1[0] * v2[1] - v1[1] * v2[0];
    // Now compute the minors of rank 3.
    return m01 * v2[2] - m02 * v1[2] + m12 * v0[2];
}

inline double distance(Point const &p0, Point const &p1) {
    double const dx = p0[0] - p1[0];
    double const dy = p0[1] - p1[1];
    double const dz = p0[2] - p1[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline bool equal(Point const &p, CPoint const &q) {
    return (deq(p[0], q.x()) && deq(p[1], q.y()) && deq(p[2], q.z()));
}

inline bool equal(CPoint const &q, Point const &p) {
    return (deq(p[0], q.x()) && deq(p[1], q.y()) && deq(p[2], q.z()));
}

inline double volume(
    const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
    return std::abs(determinant(p1 - p0, p2 - p0, p3 - p0)) / 6;
}

inline double volume(const Tetrahedron &t) {
    return volume(t[0], t[1], t[2], t[3]);
}

// Get the normal vector of the plane specified by p0, p1, p2.
inline Vector get_normal(Point const p0, Point const p1, Point const p2) {

    Vector const plane_vec_1 = p1 - p0;
    Vector const plane_vec_2 = p2 - p1;
    return normalize(cross_product(plane_vec_1, plane_vec_2));
}

// Get the volume ocuppied by the sector of the sphere inscribed in the
// incident cell.
double sphere_sector_vol(Point const &p0, Point const &p1, Point const &p2,
    Point const &p3, double const radius);

} // namespace ANA

#endif // _H