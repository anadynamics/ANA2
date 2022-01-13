// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_UNIT_CELL_HPP
#define CHEMFILES_UNIT_CELL_HPP

#include "chemfiles/types.hpp"
#include "chemfiles/exports.hpp"
#include "chemfiles/config.hpp"

#ifdef CHEMFILES_WINDOWS
#undef INFINITE
#endif

namespace chemfiles {

/// An UnitCell represent the box containing the atoms, and its periodicity
///
/// A unit cell is fully represented by three lenghts (a, b, c); and three angles
/// (alpha, beta, gamma). The angles are stored in degrees, and the lenghts in
/// Angstroms.
///
/// A cell also has a matricial representation, by projecting the three base
/// vector into an orthonormal base. We choose to represent such matrix as an
/// upper triangular matrix:
///
/// ```
/// | a_x   b_x   c_x |
/// |  0    b_y   c_y |
/// |  0     0    c_z |
/// ```
class CHFL_EXPORT UnitCell {
public:
    /// Possible shapes for the unit cell
    enum CellShape {
        /// Orthorhombic cell, with the three angles equals to 90°
        ORTHORHOMBIC = 0,
        /// Triclinic cell, with any values for the angles.
        TRICLINIC = 1,
        /// Infinite cell, to use when there is no cell
        INFINITE = 2
    };

    /// Copy constructor
    UnitCell(const UnitCell& other) = default;
    UnitCell& operator=(const UnitCell& other) = default;
    /// Move constructor
    UnitCell(UnitCell&& other) = default;
    UnitCell& operator=(UnitCell&& other) = default;

    /// Construct an INFINITE unit cell
    UnitCell();
    /// Construct a cubic unit cell of side size `a`
    UnitCell(double a);
    /// Construct an ORTHOROMBIC unit cell of side size `a`, `b`, `c`
    UnitCell(double a, double b, double c);
    /// Construct a TRICLINIC unit cell of side size `a`, `b`, `c`, and cell
    /// angles `alpha`, `beta`, `gamma`
    UnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    /// Construct a cell of type `type`, with all lenghts set to 0 and all
    /// angles set to 90°
    UnitCell(CellShape shape);
    /// Construct a cell of type `type`, with all lenghts set to `a` and all
    /// angles set to 90°
    UnitCell(CellShape shape, double a);
    /// Construct a cell of type `type`, with lenghts set to `a` ,`b`, `d`,
    /// and all angles set to 90°
    UnitCell(CellShape shape, double a, double b, double c);

    /// Get a matricial representation of the cell.
    Matrix3D matricial() const {
        return h_;
    }
    /// Populate C-style matricial representation of the cell. The array should
    /// have a 3 x 3 size.
    void raw_matricial(double[3][3]) const;

    /// Get the cell shape
    CellShape shape() const { return shape_; }
    /// Set the cell shape to `shape`
    void shape(CellShape shape);

    /// Get the first lenght (a) of the cell
    double a() const { return a_; }
    /// Set the first lenght (a) of the cell
    void set_a(double val);
    /// Get the second lenght (b) of the cell
    double b() const { return b_; }
    /// Set the second lenght (b) of the cell
    void set_b(double val);
    /// Get the third lenght (c) of the cell
    double c() const { return c_; }
    /// Set the third lenght (c) of the cell
    void set_c(double val);

    /// Get the first angle (alpha) of the cell
    double alpha() const { return alpha_; }
    /// Set the first angle (alpha) of the cell if possible
    void set_alpha(double val);
    /// Get the second angle (beta) of the cell
    double beta() const { return beta_; }
    /// Set the second angle (beta) of the cell if possible
    void set_beta(double val);
    /// Get the third angle (gamma) of the cell
    double gamma() const { return gamma_; }
    /// Set the third angle (gamma) of the cell if possible
    void set_gamma(double val);

    /// Get the unit cell volume
    double volume() const;

    /// Wrap the vector `vect` in the unit cell, using periodic boundary
    /// conditions.
    Vector3D wrap(const Vector3D& vect) const;

private:
    /// Wrap a vector in orthorombic cell
    Vector3D wrap_orthorombic(const Vector3D& vect) const;
    /// Wrap a vector in triclinic cell
    Vector3D wrap_triclinic(const Vector3D& vect) const;
    /// Compute the cell matrix from the cell parameters
    void update_matrix();
    /// Caching the cell matrix
    Matrix3D h_;
    /// Caching the inverse of the cell matrix
    Matrix3D h_inv_;

    /// Cell lenghts
    double a_, b_, c_;
    /// Cell angles
    double alpha_, beta_, gamma_;
    /// Cell type
    CellShape shape_;
};

/// Exact comparison of unit cells.
///
/// This performs an exact comparison of the two unit cells, using floating
/// point equality. This means that the two cells have to be exactly identical,
/// not only very close.
CHFL_EXPORT bool operator==(const UnitCell& rhs, const UnitCell& lhs);
CHFL_EXPORT bool operator!=(const UnitCell& rhs, const UnitCell& lhs);

} // namespace chemfiles

#endif
