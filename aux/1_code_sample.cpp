#include <array>
#include <iostream>
#include <cmath>

class Vector {
public:
    Vector() = default;

    Vector(double const x, double const y, double const z) noexcept :
        _vxyz {x, y, z}, _origin {0., 0., 0.} { }

    Vector(double const x, double const y, double const z, double const ox,
        double const oy, double const oz) noexcept :
        _vxyz {x, y, z},
        _origin {ox, oy, oz} { }

    double operator[](int const idx) const { return _vxyz[idx]; }

    double &operator[](int const idx) { return _vxyz[idx]; }

    std::array<double, 3> _vxyz, _origin;
};

inline std::ostream &operator<<(std::ostream &stream, Vector const &v) {
    stream << v[0] << " " << v[1] << " " << v[2];
    return stream;
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator+(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
        lhs._origin[0], lhs._origin[1], lhs._origin[2]);
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator-(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
        lhs._origin[0], lhs._origin[1], lhs._origin[2]);
}

inline Vector operator/(Vector const &lhs, double const rhs) {
    if (rhs > 0.001 || rhs < -0.001) {
        double const v0 =
            (lhs[0] > 0.001 || lhs[0] < -0.001) ? lhs[0] / rhs : lhs[0];
        double const v1 =
            (lhs[1] > 0.001 || lhs[1] < -0.001) ? lhs[1] / rhs : lhs[1];
        double const v2 =
            (lhs[2] > 0.001 || lhs[2] < -0.001) ? lhs[2] / rhs : lhs[2];
        return {v0, v1, v2, lhs._origin[0], lhs._origin[1], lhs._origin[2]};
    } else {
        return lhs;
    }
}

inline Vector operator*(Vector const &lhs, double const rhs) {
    return {lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs, lhs._origin[0],
        lhs._origin[1], lhs._origin[2]};
}

inline Vector operator*(double const lhs, Vector const &rhs) {
    return {rhs[0] * lhs, rhs[1] * lhs, rhs[2] * lhs, rhs._origin[0],
        rhs._origin[1], rhs._origin[2]};
}

inline bool operator==(Vector const &lhs, Vector const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
        lhs._origin[0] == rhs._origin[0] && lhs._origin[1] == rhs._origin[1] &&
        lhs._origin[2] == rhs._origin[2]);
}

class Point {
public:
    Point() = default;

    Point(double const x, double const y, double const z) : _xyz {x, y, z} { }

    

    double operator[](int const idx) const { return _xyz[idx]; }

    std::array<double, 3> _xyz;
};

inline std::ostream &operator<<(std::ostream &stream, Point const &p) {
    stream << p[0] << " " << p[1] << " " << p[2];
    return stream;
}

// Returns a Vector starting on this Point coordinates.
inline Vector operator-(Point const &lhs, Point const &rhs) {
    return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], rhs[0], rhs[1],
        rhs[2]};
}

// Use point as vector.
inline Point operator+(Point const &p0, Point const &p1) {
    return Point(p0[0] + p1[0], p0[1] + p1[1], p0[2] + p1[2]);
}

// Displaces the Point along the Vector.
inline Point operator+(Point const &p, Vector const &v) {
    return Point(p[0] + v[0], p[1] + v[1], p[2] + v[2]);
}

// Displaces the Point along the Vector.
inline Point operator-(Point const &p, Vector const &v) {
    return Point(p[0] - v[0], p[1] - v[1], p[2] - v[2]);
}

// Displaces the Point along the Vector.
inline Point operator+(Vector const &v, Point const &p) { return p + v; }

// Displaces the Point along the Vector.
inline Point operator-(Vector const &v, Point const &p) { return p - v; }

inline Point operator*(Point const &p, double const &rhs) {
    return {p[0] * rhs, p[1] * rhs, p[2] * rhs};
}

inline Point operator*(double const &rhs, Point const &p) {
    return {p[0] * rhs, p[1] * rhs, p[2] * rhs};
}

inline bool operator==(Point const &lhs, Point const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
}

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

inline double volume(
    const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
    return std::abs(determinant(p1 - p0, p2 - p0, p3 - p0)) / 6;
}

int main() {
    Point p0(-2.0439999098885959, -3.062999963695026, 13.878000262654703);
    Point p1(-2.0550000667572021, -0.61100000143051147, 1.2979999780654907);
    Point ip2(-6.2881150941796635, -1.4487250377584024, 8.5971887436084913);
    Point ip3(-7.2606551334607587, -1.0200050533985316, 9.63917576630433);

    std::cout << volume(p0, p1, ip2, ip3) << '\n';



    return 0;
}