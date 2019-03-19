#ifndef ANA_PRIMITIVES_H
#define ANA_PRIMITIVES_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>

namespace ANA {

class CCTetrahedron {
public:
    CCTetrahedron() noexcept = default;

    CCTetrahedron(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &p3) :
        _data({p0, p1, p2, p3}) {}

    CCTetrahedron(CPoint &&p0, CPoint &&p1, CPoint &&p2, CPoint &&p3) :
        _data({p0, p1, p2, p3}) {}

    std::array<CPoint, 4> _data;
};

struct TetraInfo {
public:
    TetraInfo() = default;

    std::array<int, 4> _index;
    std::array<double, 4> _radius;
    // atom's residue number
    std::array<int, 4> _resn;
    // atom's residue name in 3 letter format
    std::array<std::string, 4> _resi;
};

class CTriangularPrism {
public:
    CTriangularPrism() noexcept = default;

    CTriangularPrism(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &p3, CPoint const &p4, CPoint const &p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    CTriangularPrism(CPoint &&p0, CPoint &&p1, CPoint &&p2, CPoint &&p3,
        CPoint &&p4, CPoint &&p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    CPoint const &operator[](int const idx) const { return _data[idx]; }

    CPoint &operator[](int const idx) { return _data[idx]; }

    std::array<CPoint, 6> _data;
};

class Vector {
public:
    Vector() = default;

    Vector(double const x, double const y, double const z) noexcept :
        _vxyz{x, y, z}, _origin{0., 0., 0.} {}

    Vector(double const x, double const y, double const z, double const ox,
        double const oy, double const oz) noexcept :
        _vxyz{x, y, z},
        _origin{ox, oy, oz} {}

    Vector(CVector const v) :
        _vxyz({CGAL::to_double(v.x()), CGAL::to_double(v.y()),
            CGAL::to_double(v.z())}) {}

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

    Point(double const x, double const y, double const z) : _xyz{x, y, z} {}

    Point(CPoint const p) :
        _xyz{CGAL::to_double(p.x()), CGAL::to_double(p.y()),
            CGAL::to_double(p.z())} {}

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
// Displaces the Point along the Vector.
inline Point operator+(Point const &p, Vector const &v) {
    return Point(p[0] + (v[0] - v._origin[0]), p[1] + (v[1] - v._origin[1]),
        p[2] + (v[2] - v._origin[2]));
}
// Displaces the Point along the Vector.
inline Point operator-(Point const &p, Vector const &v) {
    return Point(p[0] - (v[0] - v._origin[0]), p[1] - (v[1] - v._origin[1]),
        p[2] - (v[2] - v._origin[2]));
}

inline bool operator==(Point const &lhs, Point const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
}

inline Vector point_to_vector(Point const &in_point) {
    return Vector(in_point[0], in_point[1], in_point[2]);
}

class TTriangle {
public:
    TTriangle() noexcept = default;

    TTriangle(Point const &p0, Point const &p1, Point const &p2) :
        _data({p0, p1, p2}) {}

    TTriangle(CPoint const &p0, CPoint const &p1, CPoint const &p2) :
        _data({Point(p0), Point(p1), Point(p2)}) {}

    TTriangle(Point &&p0, Point &&p1, Point &&p2) : _data({p0, p1, p2}) {}

    Point const &operator[](int const idx) const { return _data[idx]; }

    Point &operator[](int const idx) { return _data[idx]; }

    std::array<Point, 3> _data;
};

class TTetrahedron {
public:
    TTetrahedron() noexcept = default;

    TTetrahedron(
        Point const &p0, Point const &p1, Point const &p2, Point const &p3) :
        _data({p0, p1, p2, p3}) {}

    TTetrahedron(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &p3) :
        _data({Point(p0), Point(p1), Point(p2), Point(p3)}) {}

    TTetrahedron(Point &&p0, Point &&p1, Point &&p2, Point &&p3) :
        _data({p0, p1, p2, p3}) {}

    TTetrahedron(CPoint const &&p0, CPoint const &&p1, CPoint const &&p2,
        CPoint const &&p3) :
        _data({Point(p0), Point(p1), Point(p2), Point(p3)}) {}

    TTetrahedron(Finite_cells_iterator const &cell) :
        _data({cell->vertex(0)->point(), cell->vertex(1)->point(),
            cell->vertex(2)->point(), cell->vertex(3)->point()}) {}

    Point const &operator[](int const idx) const { return _data[idx]; }

    Point &operator[](int const idx) { return _data[idx]; }

    std::array<Point, 4> _data;
};

class TTriangularPrism {
public:
    TTriangularPrism() noexcept = default;

    TTriangularPrism(Point const &p0, Point const &p1, Point const &p2,
        Point const &p3, Point const &p4, Point const &p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    TTriangularPrism(Point &&p0, Point &&p1, Point &&p2, Point &&p3, Point &&p4,
        Point &&p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    TTriangularPrism(CPoint const &p0, CPoint const &p1, CPoint const &p2,
        CPoint const &p3, CPoint const &p4, CPoint const &p5) :
        _data({Point(p0), Point(p1), Point(p2), Point(p3), Point(p4),
            Point(p5)}) {}

    TTriangularPrism(CPoint const &&p0, CPoint const &&p1, CPoint const &&p2,
        CPoint const &&p3, CPoint const &&p4, CPoint const &&p5) :
        _data({Point(p0), Point(p1), Point(p2), Point(p3), Point(p4),
            Point(p5)}) {}

    Point const &operator[](int const idx) const { return _data[idx]; }

    Point &operator[](int const idx) { return _data[idx]; }

    std::array<Point, 6> _data;
};

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

inline double volume(
    const Point &p0, const Point &p1, const Point &p2, const Point &p3) {
    return determinant(p1 - p0, p2 - p0, p3 - p0) / 6;
}

inline double volume(const TTetrahedron &t) {
    return volume(t[0], t[1], t[2], t[3]);
}
}

#endif // _H