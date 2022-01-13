// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_RESIDUE_HPP
#define CHEMFILES_RESIDUE_HPP

#include <string>
#include <algorithm>

#include "chemfiles/exports.hpp"
#include "chemfiles/optional.hpp"
#include "chemfiles/sorted_set.hpp"

namespace chemfiles {

/// A groupement of atoms in the same logical unit.
///
/// A `Residue` is a group of atoms belonging to the same logical unit. They
/// can be small molecules, amino-acids in a protein, monomers in polymers,
/// etc.
///
/// Iterating over a `Residue` gives the indexes of the atoms in the
/// associated `Topology`.
class CHFL_EXPORT Residue final {
public:
    /// Create a new residue with a given `name`.
    explicit Residue(std::string name);

    /// Create a new residue with a given `name` and residue id `resid`.
    Residue(std::string name, uint64_t resid);

    Residue(const Residue&) = default;
    Residue& operator=(const Residue&) = default;
    Residue(Residue&&) = default;
    Residue& operator=(Residue&&) = default;

    /// Get the name of the residue
    const std::string& name() const {
        return name_;
    }

    /// Get the residue identifier if it exists.
    optional<uint64_t> id() const {
        return id_;
    }
    /// Get the size of the residue, i.e. the number of atoms in this residue.
    size_t size() const {
        return atoms_.size();
    }

    /// Add an atom with index `i` to this residue
    void add_atom(size_t i);
    /// Check if the residue contains a given atom with index `i`
    bool contains(size_t i) const;

    using const_iterator = sorted_set<size_t>::const_iterator;
    // Iterators over the indexes of the atoms in the residue
    const_iterator begin() const {return atoms_.begin();}
    const_iterator end() const {return atoms_.end();}
    const_iterator cbegin() const {return atoms_.cbegin();}
    const_iterator cend() const {return atoms_.cend();}

private:
    /// Name of the residue
    std::string name_;
    /// Index of the residue in the initial topology file
    optional<uint64_t> id_;
    /// Indexes of the atoms in this residue. These indexes refers to the
    /// associated topology.
    sorted_set<size_t> atoms_;
};

inline bool operator==(const Residue& lhs, const Residue& rhs) {
    return lhs.id() == rhs.id() && lhs.name() == rhs.name() &&
           lhs.size() == rhs.size() && std::equal(
                lhs.begin(), lhs.end(), rhs.begin()
           );
}

inline bool operator!=(const Residue& lhs, const Residue& rhs) {
    return !(lhs == rhs);
}

} // namespace chemfiles

#endif
