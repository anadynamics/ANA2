// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_FACTORY_HPP
#define CHEMFILES_FORMAT_FACTORY_HPP

#include <memory>
#include <functional>
#include <unordered_map>

#include "chemfiles/exports.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/Error.hpp"

namespace chemfiles {

typedef std::function<std::unique_ptr<Format>(const std::string& path, File::Mode mode)> format_creator_t;

/// This class allow to register Format with names and file extensions
class CHFL_EXPORT FormatFactory {
private:
    FormatFactory();
    /// Files extensions to trajectory builder associations
    using trajectory_map_t = std::unordered_map<std::string, format_creator_t>;

public:
    /// Get the instance of the TrajectoryFactory
    static FormatFactory& get();

    /// @brief Get a `format_creator_t` from a format name.
    /// @param name the format name
    /// @return A `format_creator_t` corresponding to the format, if the format
    ///         name is found in the list of registered formats.
    ///
    /// @throws FormatError if the format can not be found
    format_creator_t name(const std::string& name);

    /// @brief Get a `format_creator_t` from a format extention.
    /// @param extension the format extention
    /// @return A `format_creator_t` corresponding to the format, if the format
    ///         extension is found in the list of registered extensions.
    ///
    /// @throws FormatError if the format can not be found
    format_creator_t extension(const std::string& extension);

    /// @brief Register a format `F` with the given `name`
    /// @param name the format name
    /// @throws FormatError if the name is already used by another format
    ///
    /// Usage example:
    /// ```
    /// FormatFactory::get().register_name<MyFormat>("my name");
    /// ```
    template <typename F>
    void register_name(const std::string& name) {
        if (formats_.find(name) != formats_.end()) {
            throw FormatError(
                "The name '" + name + "' is already associated with a format."
            );
        }
        formats_.emplace(name, [](const std::string& path, File::Mode mode) {
            return std::unique_ptr<Format>(new F(path, mode));
        });
    }

    /// @brief Register a format `F` with the given `extension`
    /// @param extension the format extension
    /// @throws FormatError if the extension is already used by another format
    ///
    /// Usage example:
    /// ```
    /// FormatFactory::get().register_extension<MyFormat>(".mft");
    /// ```
    template <typename F>
    void register_extension(const std::string& extension) {
        if (extensions_.find(extension) != extensions_.end()) {
            throw FormatError(
                "The extension '" + extension + "' is already associated with a format."
            );
        }
        extensions_.emplace(extension, [](const std::string& path, File::Mode mode) {
            return std::unique_ptr<Format>(new F(path, mode));
        });
    }

private:
    /// Trajectory map associating format descriptions and readers
    trajectory_map_t formats_;
    /// Trajectory map associating format descriptions and readers
    trajectory_map_t extensions_;
};

} // namespace chemfiles

#endif
