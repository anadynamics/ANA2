// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_CHFL_PROPERTY_H
#define CHEMFILES_CHFL_PROPERTY_H

#include "chemfiles/capi/types.h"
#ifdef __cplusplus
extern "C" {
#endif

/// Possible values holded by a CHFL_PROPERTY
typedef enum {
    /// Bool value
    CHFL_PROPERTY_BOOL = 0,
    /// Double value
    CHFL_PROPERTY_DOUBLE = 1,
    /// String value
    CHFL_PROPERTY_STRING = 2,
    /// chfl_vector3d value
    CHFL_PROPERTY_VECTOR3D = 3,
} chfl_property_kind;

/// Create a new property holding a boolean `value`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_property_free`.
///
/// @example{tests/capi/doc/chfl_property/bool.c}
/// @return A pointer to the property, or NULL in case of error. You can use
///         `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_PROPERTY* chfl_property_bool(bool value);

/// Create a new property holding a double `value`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_property_free`.
///
/// @example{tests/capi/doc/chfl_property/double.c}
/// @return A pointer to the property, or NULL in case of error. You can use
///         `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_PROPERTY* chfl_property_double(double value);

/// Create a new property holding a string `value`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_property_free`.
///
/// @example{tests/capi/doc/chfl_property/string.c}
/// @return A pointer to the property, or NULL in case of error. You can use
///         `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_PROPERTY* chfl_property_string(const char* value);

/// Create a new property holding a 3D vector `value`.
///
/// The caller of this function should free the allocated memory using
/// `chfl_property_free`.
///
/// @example{tests/capi/doc/chfl_property/vector3d.c}
/// @return A pointer to the property, or NULL in case of error. You can use
///         `chfl_last_error` to learn about the error.
CHFL_EXPORT CHFL_PROPERTY* chfl_property_vector3d(const chfl_vector3d value);

/// Get the type of value holded by this `property` in `kind`.
///
/// @example{tests/capi/doc/chfl_property/kind.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_property_get_kind(
    const CHFL_PROPERTY* const property, chfl_property_kind* kind
);

/// Get the boolean value holded by this `property` in the location pointed to
/// by `data`.
///
/// This function returns CHFL_PROPERTY_ERROR if the property is not a boolean
/// property.
///
/// @example{tests/capi/doc/chfl_property/bool.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_property_get_bool(
    const CHFL_PROPERTY* const property, bool* data
);

/// Get the double value holded by this `property` in the location pointed to
/// by `data`.
///
/// This function returns CHFL_PROPERTY_ERROR if the property is not a double
/// property.
///
/// @example{tests/capi/doc/chfl_property/double.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_property_get_double(
    const CHFL_PROPERTY* const property, double* data
);

/// Get the string value holded by this `property` in the given `buffer`.
///
/// This function returns CHFL_PROPERTY_ERROR if the property is not a sring
/// property.
///
/// The buffer size must be passed in `buffsize`. This function will truncate
/// the property to fit in the buffer.
///
/// @example{tests/capi/doc/chfl_property/string.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_property_get_string(
    const CHFL_PROPERTY* const property, char* const buffer, uint64_t buffsize
);

/// Get the 3D vector value holded by this `property` in the location pointed to
/// by `data`.
///
/// This function returns CHFL_PROPERTY_ERROR if the property is not a 3D vector
/// property.
///
/// @example{tests/capi/doc/chfl_property/vector3d.c}
/// @return The operation status code. You can use `chfl_last_error` to learn
///         about the error if the status code is not `CHFL_SUCCESS`.
CHFL_EXPORT chfl_status chfl_property_get_vector3d(
    const CHFL_PROPERTY* const property, chfl_vector3d data
);

/// Free the memory associated with a `property`.
///
/// @example{tests/capi/doc/chfl_property/bool.c}
/// @return `CHFL_SUCCESS`
CHFL_EXPORT chfl_status chfl_property_free(CHFL_PROPERTY* property);

#ifdef __cplusplus
}
#endif

#endif
