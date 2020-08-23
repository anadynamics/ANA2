#ifndef ANA_SOA_PRIMITIVES_H
#define ANA_SOA_PRIMITIVES_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>

namespace ANA {

struct Points {
public:
    Points() = default;

    Points(int const i) : _cnt(i) {
        _x.reserve(_cnt);
        _y.reserve(_cnt);
        _z.reserve(_cnt);
    }

    void reserve(int const i) {
        _cnt = i;
        _x.reserve(_cnt);
        _y.reserve(_cnt);
        _z.reserve(_cnt);

        return;
    }

    void push(double const x, double const y, double const z) {
        _x.push_back(x);
        _y.push_back(y);
        _z.push_back(z);

        return;
    }

    int _cnt;
    std::vector<double> _x, _y, _z;
};

struct Triangles {
public:
    Triangles() = default;

    Triangles(int const i) : _cnt(i) {
        _x.reserve(_cnt);
        _y.reserve(_cnt);
        _z.reserve(_cnt);
    }

    void reserve(int const i) {
        _cnt = i;
        _x.reserve(_cnt);
        _y.reserve(_cnt);
        _z.reserve(_cnt);

        return;
    }

    void push(double const x, double const y, double const z) {
        _x.push_back(x);
        _y.push_back(y);
        _z.push_back(z);

        return;
    }

    int _cnt;
    std::vector<double> _x, _y, _z;
};
}

#endif // _H