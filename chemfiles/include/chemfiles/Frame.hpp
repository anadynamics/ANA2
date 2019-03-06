// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FRAME_HPP
#define CHEMFILES_FRAME_HPP

#include "chemfiles/exports.hpp"
#include "chemfiles/optional.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/span.hpp"

#include "chemfiles/Topology.hpp"
#include "chemfiles/UnitCell.hpp"

namespace chemfiles {

/// A frame contains data from one simulation step The Frame class holds data
/// from one step of a simulation: the current topology, the positions, and the
/// velocities of the particles in the system.  If some information is missing
/// (topology or velocity or unit cell), the corresponding data is filled with a
/// default value. Specifically:
///
/// * `velocities` is the `nullopt` version of
///   `optional<std::vector<Vector3D>>`. Here, `optional<T>` refers to the
///   [std::optional] class as defined in C++17.
/// * `cell` is an infinite unit cell;
/// * `topology` is empty, and contains no data.
///
/// [std::optional]: http://en.cppreference.com/w/cpp/optional
class CHFL_EXPORT Frame {
public:
    /// Default constructor
    Frame();
    /// Constructor reserving some space for `natoms`
    explicit Frame(size_t natoms);
    /// Constructor reserving space for `topology.natoms()`, and using `cell`
    /// as unit cell. `cell` default to an `INFINITE` unit cell.
    explicit Frame(Topology topology, UnitCell cell = UnitCell());

    Frame(Frame&&) = default;
    Frame& operator=(Frame&&) = default;

    /// Get a clone (exact copy) of this frame.
    ///
    /// This replace the implicit copy constructor (which is disabled) to
    /// make an explicit copy of the frame.
    Frame clone() const {
        return *this;
    }

    /// Get a modifiable reference to the positions
    span<Vector3D> positions() { return positions_; }
    /// Get a const (non modifiable) reference to the positions
    const std::vector<Vector3D>& positions() const { return positions_; }

    /// Get an optional modifiable reference to the velocities
    optional<span<Vector3D>> velocities() {
        return velocities_ ? optional<span<Vector3D>>(as_span(*velocities_))
                           : optional<span<Vector3D>>(nullopt);
    }
    /// Get an optional const (non modifiable) reference to the velocities
    const optional<std::vector<Vector3D>>& velocities() const { return velocities_; }
    /// Add velocities to this frame. If velocities are already defined,
    /// this functions does nothing.
    void add_velocities();

    /// Get a modifiable reference to the internal topology
    Topology& topology() { return topology_; }
    /// Get a const (non-modifiable) reference to the internal topology
    const Topology& topology() const { return topology_; }
    /// Set the system topology
    void set_topology(const Topology& topology);

    /// Get a const (non-modifiable) reference to the unit cell of the system
    const UnitCell& cell() const { return cell_; }
    UnitCell& cell() { return cell_; }
    /// Set the unit cell fo the system
    void set_cell(const UnitCell& c) { cell_ = c; }

    /// Resize the frame to store data for `natoms` atoms. If the new size is
    /// bigger than the old one, missing data is initializd to 0. Pre-existing
    /// values are conserved.
    /// This function only resize the velocities if the data is present.
    void resize(size_t natoms);

    /// Reserve size in the frame to store data for `natoms` atoms.
    /// This function only reserve storage for the the velocities if the data
    /// is present.
    void reserve(size_t natoms);

    /// Add an `atom` at the given `position` and optionally with the given
    /// `velocity`. The `velocity` value will only be used if this frame
    /// contains velocity data.
    void add_atom(Atom atom, Vector3D position, Vector3D velocity = Vector3D());

    /// Get the number of atoms in the system
    size_t natoms() const;

    /// Remove the atom at index `i` in the system.
    ///
    /// @throws chemfiles::OutOfBounds if `i` is not in bounds
    void remove(size_t i);

    /// Get the current simulation step
    size_t step() const { return step_; }
    /// Set the current simulation step
    void set_step(size_t s) { step_ = s; }

    /// Guess the bonds, angles and dihedrals in the system. The bonds are
    /// guessed using a distance-based algorithm, and then angles and dihedrals
    /// are guessed from the bonds.
    void guess_topology();

    /// Get the distance between the atoms at indexes `i` and `j`, accounting
    /// for periodic boundary conditions. The distance is expressed in angstroms.
    ///
    /// @throws chemfiles::OutOfBounds if `i` or `j` are not in bounds
    double distance(size_t i, size_t j) const;

    /// Get the angle formed by the atoms at indexes `i`, `j` and `k`,
    /// accounting for periodic boundary conditions. The angle is expressed in
    /// radians.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j` or `k` are not in bounds
    double angle(size_t i, size_t j, size_t k) const;

    /// Get the dihedral angle formed by the atoms at indexes `i`, `j`, `k` and
    /// `m`, accounting for periodic boundary conditions. The angle is expressed
    /// in radians.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j`, `k` or `m` are not in bounds
    double dihedral(size_t i, size_t j, size_t k, size_t m) const;

    /// Get the out of plane distance formed by the atoms at indexes `i`, `j`,
    /// `k` and `m`, accounting for periodic boundary conditions. The distance
    /// is expressed in angstroms.
    ///
    /// This is the distance betweent the atom j and the ikm plane. The j atom
    /// is the center of the improper dihedral angle formed by i, j, k and m.
    ///
    /// @throws chemfiles::OutOfBounds if `i`, `j`, `k` or `m` are not in bounds
    double out_of_plane(size_t i, size_t j, size_t k, size_t m) const;

    /// Set an arbitrary property for this frame with the given `name` and
    /// `value`. If a property with this name already exist, it is replaced
    /// with the new value.
    void set(std::string name, Property value);

    /// Get the property with the given `name` for this frame if it exists.
    optional<const Property&> get(const std::string& name) const;

private:
    Frame(const Frame&) = default;
    Frame& operator=(const Frame&) = default;

    /// Current simulation step
    size_t step_;
    /// Positions of the particles
    std::vector<Vector3D> positions_;
    /// Velocities of the particles
    optional<std::vector<Vector3D>> velocities_;
    /// Topology of the described system
    Topology topology_;
    /// Unit cell of the system
    UnitCell cell_;
    /// Properties stored in this frame
    property_map properties_;
};

} // namespace chemfiles

#endif
