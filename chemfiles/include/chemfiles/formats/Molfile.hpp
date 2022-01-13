// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_MOLFILE_HPP
#define CHEMFILES_FORMAT_MOLFILE_HPP

extern "C" {
#include "vmdplugin.h"
#include "molfile_plugin.h"
}

#include "chemfiles/optional.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/Topology.hpp"

namespace chemfiles {

/// List all the VMD molfile plugins enabled. For more documentation about VMD
/// molfile plugins, please see:
/// http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/
enum MolfileFormat {
    DCD, ///< DCD binary file format
    GRO, ///< Gromacs .gro file format
    TRR, ///< Gromacs .trr file format
    XTC, ///< Gromacs .xtc file format
    TRJ, ///< Gromacs .trj file format
    LAMMPS, ///< Lammps trajectory files
};

/// A thin wrapper around the vmd plugin functions
template <MolfileFormat F>
struct MolfilePluginData {
    /// Function to initialize the plugin
    int init();
    /// Function to register the plugin.
    int registration(void* data, vmdplugin_register_cb callback);
    /// Function to unload the plugin
    int fini();

    /// Plugin format
    std::string format() const;
    /// Plugin name
    std::string plugin_name() const;
    /// Plugin reader, when a given plugin manages different formats
    std::string reader() const;
    /// Does this plugin have velocity data
    bool have_velocities() const;
};

/// Use of VMD Molfile plugins as format reader. This class is templated by a
/// value in the `MolfileFormat` enum.
template <MolfileFormat F>
class Molfile final: public Format {
public:
    Molfile(const std::string& path, File::Mode mode);
    ~Molfile() noexcept;

    void read(Frame& frame) override;
    size_t nsteps() override;
private:
    /// Convert a molfile timestep to a chemfiles frame
    void molfile_to_frame(const molfile_timestep_t& timestep, Frame& frame);
    /// Read topological information in the current file, if any.
    void read_topology();

    /// Path of the underlying file
    std::string path_;
    /// VMD plugin data
    MolfilePluginData<F> plugin_data_;
    /// VMD molfile plugin
    molfile_plugin_t* plugin_handle_;
    /// The file handler
    void* file_handle_;
    /// The number of atoms in the last trajectory read
    int natoms_;
    /// Store optional topological information
    optional<Topology> topology_;
};

} // namespace chemfiles

#endif
