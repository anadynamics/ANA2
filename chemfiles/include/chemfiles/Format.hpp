// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_HPP
#define CHEMFILES_FORMAT_HPP

#include <memory>
#include <string>

#include "chemfiles/exports.hpp"

namespace chemfiles {
class Frame;

/// The `Format` class defines the interface to implement in order to add a new
/// format to chemfiles. It is possible to implement only one of `Format::read`;
/// `Format::read_step` or `Format::write`. In that case, only the corresponding
/// operations will be available from the corresponding `chemfiles::Trajectory`.
class CHFL_EXPORT Format {
public:
    Format() = default;
    virtual ~Format() noexcept = default;
    Format& operator=(const Format&) = delete;
    Format(const Format&) = delete;

    /// @brief Read a specific step from the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param step The step to read
    /// @param frame The frame to fill
    virtual void read_step(size_t step, Frame& frame);

    /// @brief Read a specific step from the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param frame The frame to fill
    virtual void read(Frame& frame);

    /// @brief Write a frame to the trajectory file.
    ///
    /// @throw FormatError if the file does not follow the format
    /// @throw FileError if their is an OS error while reading the file
    ///
    /// @param frame The frame to be writen
    virtual void write(const Frame& frame);

    /// @brief Get the number of frames in the associated file
    /// @return The number of frames
    virtual size_t nsteps() = 0;
};

} // namespace chemfiles

#endif
