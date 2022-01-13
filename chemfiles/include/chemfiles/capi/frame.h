// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_CHFL_FRAME_H
#define CHEMFILES_CHFL_FRAME_H

#include "chemfiles/capi/types.h"
#ifdef __cplusplus
extern "C" {
#endif

/// Create a new empty frame.
///
/// The caller of this function should free the allocated memory using
/// `chfl_frame_free`.
///
/// @example{tests/capi/doc/chfl_frame/chfl_frame.c}
/// @return A pointer to the frame, or NULL in case of error. You can use
///         `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_FRAME* chfl_frame(void);

/// Get a copy of a `frame`.
///
/// The caller of this function should free the associated memory using
/// `chfl_frame_free`.
///
/// @example{tests/capi/doc/chfl_frame/copy.c}
/// @return A pointer to the new frame, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_FRAME* chfl_frame_copy(const CHFL_FRAME* const frame);

/// Get the current number of atoms in a `frame` in the integer pointed to by
/// `natoms`
///
/// @example{tests/capi/doc/chfl_frame/atoms_count.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_atoms_count(
    const CHFL_FRAME* const frame, uint64_t* natoms
);

/// Get a pointer to the positions array from a `frame`.
///
/// Positions are stored as a `natoms x 3` array, this function set the pointer
/// pointed to by `positions` to point to the first element of this array, and
/// give the number of atoms in the integer pointed to by `size`.
///
/// If the frame is resized (by writing to it, or calling `chfl_frame_resize`,
/// `chfl_frame_remove` or `chfl_frame_add_atom`), the pointer is invalidated.
/// If the frame is freed using `chfl_frame_free`, the pointer is freed too.
/// There is then no need to free the `*positions` pointer for the caller of
/// this function.
///
/// @example{tests/capi/doc/chfl_frame/positions.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_positions(
    CHFL_FRAME* const frame, chfl_vector3d** positions, uint64_t* size
);

/// Get a pointer to the velocities array from a `frame`.
///
/// Velocities are stored as a `natoms x 3` array, this function set the pointer
/// pointed to by `positions` to point to the first element of this array, and
/// give the number of atoms in the integer pointed to by `size`.
///
/// If the frame is resized (by writing to it, or calling `chfl_frame_resize`,
/// `chfl_frame_remove` or `chfl_frame_add_atom`), the pointer is invalidated.
/// If the frame is freed using `chfl_frame_free`, the pointer is freed too.
/// There is then no need to free the `*velocity` pointer for the caller of this
/// function.
///
/// If the frame does not have velocity, this will return an error. Use
/// `chfl_frame_add_velocities` to add velocities to a frame before calling
/// this function.
///
/// @example{tests/capi/doc/chfl_frame/velocities.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_velocities(
    CHFL_FRAME* const frame, chfl_vector3d** velocities, uint64_t* size
);

/// Add an `atom` and the corresponding `position` and `velocity` data to a
/// `frame`.
///
/// `velocity` can be `NULL` if no velocity is associated with the atom.
///
/// @example{tests/capi/doc/chfl_frame/add_atom.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_add_atom(
    CHFL_FRAME* const frame, const CHFL_ATOM* const atom,
    const chfl_vector3d position, const chfl_vector3d velocity
);

/// Remove the atom at index `i` in the `frame`.
///
/// This modify all the atoms indexes after `i`, and invalidate any pointer
/// obtained using `chfl_frame_positions` or `chfl_frame_velocities`.
///
/// @example{tests/capi/doc/chfl_frame/remove.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_remove(CHFL_FRAME* const frame, uint64_t i);

/// Resize the positions, velocities  and topology in the `frame`, to have space
/// for `natoms` atoms.
///
/// This function may invalidate any pointer to the positions or the velocities
/// if the new size is bigger than the old one. In all the cases, previous data
/// is conserved. This function conserve the presence or absence of velocities.
///
/// @example{tests/capi/doc/chfl_frame/resize.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_resize(
    CHFL_FRAME* const frame, uint64_t natoms
);

/// Add velocity data to this `frame`.
///
/// The velocities are initialized to `(chfl_vector3d){0, 0, 0}`. If the frame
/// already has velocities, this does nothing.
///
/// @example{tests/capi/doc/chfl_frame/add_velocities.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_add_velocities(CHFL_FRAME* const frame);

/// Check if this `frame` contains velocity data, and store the result in
/// `has_velocities`
///
/// @example{tests/capi/doc/chfl_frame/has_velocities.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_has_velocities(
    const CHFL_FRAME* const frame, bool* has_velocities
);

/// Set the unit cell of a `frame` to `cell`.
///
/// @example{tests/capi/doc/chfl_frame/set_cell.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_set_cell(
    CHFL_FRAME* const frame, const CHFL_CELL* const cell
);

/// Set the topology of a `frame` to `topology`.
///
/// Calling this function with a topology that does not contain the right
/// number of atom will return an error.
///
/// @example{tests/capi/doc/chfl_frame/set_topology.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_set_topology(
    CHFL_FRAME* const frame, const CHFL_TOPOLOGY* const topology
);

/// Get a `frame` step, *i.e.* the frame number in the trajectory in the integer
/// pointed to by `step`.
///
/// @example{tests/capi/doc/chfl_frame/step.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_step(
    const CHFL_FRAME* const frame, uint64_t* step
);

/// Set a `frame` step, *i.e.* the frame number in the trajectory to `step`.
///
/// @example{tests/capi/doc/chfl_frame/set_step.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_set_step(
    CHFL_FRAME* const frame, uint64_t step
);

/// Guess the bonds, angles and dihedrals in a `frame`.
///
/// The bonds are guessed using a distance-based algorithm, and then angles and
/// dihedrals are guessed from the bonds.
///
/// @example{tests/capi/doc/chfl_frame/guess_topology.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_guess_topology(CHFL_FRAME* const frame);

/// Get the distance between the atoms at indexes `i` and `j` in the `frame`,
/// accounting for periodic boundary conditions. The result is placed in
/// `distance`, and expressed in angstroms.
///
/// @example{tests/capi/doc/chfl_frame/distance.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_distance(
    const CHFL_FRAME* const frame, uint64_t i, uint64_t j, double* distance
);

/// Get the angle formed by the atoms at indexes `i`, `j` and `k` in the
/// `frame`, accounting for periodic boundary conditions. The result is placed
/// in `angle`, and expressed in radians.
///
/// @example{tests/capi/doc/chfl_frame/angle.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_angle(
    const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, double* angle
);

/// Get the angle formed by the atoms at indexes `i`, `j`, `k` and `m` in the
/// `frame`, accounting for periodic boundary conditions. The result is placed
/// in `dihedral`, and expressed in radians.
///
/// @example{tests/capi/doc/chfl_frame/dihedral.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_dihedral(
    const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, uint64_t m, double* dihedral
);

/// Get the out of plane distance formed by the atoms at indexes `i`, `j`, `k`
/// and `m` in the `frame`, accounting for periodic boundary conditions. The
/// result is placed in `distance` and expressed in angstroms.
///
/// This is the distance betweent the atom j and the ikm plane. The j atom
/// is the center of the improper dihedral angle formed by i, j, k and m.
///
/// @example{tests/capi/doc/chfl_frame/out_of_plane.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_out_of_plane(
    const CHFL_FRAME* const frame, uint64_t i, uint64_t j, uint64_t k, uint64_t m, double* distance
);

/// Add a new `property` with the given `name` to this `frame`.
///
/// If a property with the same name already exists, this function override the
/// existing property with the new one.
///
/// @example{tests/capi/doc/chfl_frame/property.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_frame_set_property(
    CHFL_FRAME* const frame, const char* name, const CHFL_PROPERTY* const property
);

/// Get a property with the given `name` in this `frame`.
///
/// This function returns `NULL` is no property exist with the given name.
///
/// The user of this function is responsible to deallocate memory using the
/// `chfl_property_free` function.
///
/// @example{tests/capi/doc/chfl_frame/property.c}
/// @return A pointer to the property, or NULL in case of error.
///         You can use `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_PROPERTY* chfl_frame_get_property(
    const CHFL_FRAME* const frame, const char* name
);

/// Free the memory associated with a `frame`.
///
/// @example{tests/capi/doc/chfl_frame/chfl_frame.c}
/// @return `CHFL_SUCCESS`
CHFL_EXPORT chfl_status chfl_frame_free(CHFL_FRAME* const frame);

#ifdef __cplusplus
}
#endif

#endif
