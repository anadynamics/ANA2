#ifndef ANA_OPTIONS_H
#define ANA_OPTIONS_H
#include <ANA/Includes.hpp>

namespace ANA {

struct InOutOptions {
    std::string _in_filename;
    std::string _in_md_filename;
    std::string _out_vol_filename = "none";
    std::string _out_pdb_filename = "none";
    std::string _out_type = "grid_pdb";
};

struct IncludedAreaOptions {
public:
    // Would like to have only 1 string, but that won't play nicely with Boost
    // Program Options. Hopefully I'll fix it someday. TODO
    std::string _resn_proto = "none";
    std::string _atom_proto = "none";
    std::string _sphere_proto = "none";
    std::string _cylinder_proto = "none";
    std::string _prism_proto = "none";
    std::string _filename = "none";
    enum IAOption { none, residue, atom, sphere, cylinder, prism, file };
    IAOption _opt = IAOption::none;
};

struct NDDOptions {
public:
    std::string _modes_ndd_filename;
    std::string _evalues_ndd_filename;
    std::string _modes_format;
    std::string _out_ndd_filename;
    bool _derivative = true;
    int _step = 5;
};

class CellFilteringOptions {
public:
    CellFilteringOptions() = default;

    CellFilteringOptions(double const minVR, double const maxSR) :
        _minVR(minVR), _maxSR(maxSR) {
        _min_CV = (4 / 3) * M_PI * _minVR * _minVR * _minVR;
        _max_FA = M_PI * _maxSR * _maxSR;
    }

    void update() {
        _min_CV = (4 / 3) * M_PI * _minVR * _minVR * _minVR;
        _max_FA = M_PI * _maxSR * _maxSR;
    }

    // min cell volume
    double _minVR, _min_CV;
    // max facet area.
    double _maxSR, _max_FA;
};

} // namespace ANA

#endif // _H