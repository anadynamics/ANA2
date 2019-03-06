// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "chemfiles/Connectivity.hpp"
#include "chemfiles/ErrorFmt.hpp"

using namespace chemfiles;

Bond::Bond(size_t i, size_t j) {
    if (i == j) {
        throw error("can not have a bond between an atom and itself");
    }

    data_[0] = std::min(i, j);
    data_[1] = std::max(i, j);
}

size_t Bond::operator[](size_t i) const {
    if (i >= 2) {
        throw out_of_bounds("can not access atom n° {} in bond", i);
    }
    return data_[i];
}

Angle::Angle(size_t i, size_t j, size_t k) {
    if (i == j || i == k || j == k) {
        throw error("can not have the same atom twice in an angle");
    }

    data_[0] = std::min(i, k);
    data_[1] = j;
    data_[2] = std::max(i, k);
}

size_t Angle::operator[](size_t i) const {
    if (i >= 3) {
        throw out_of_bounds("can not access atom n° {} in angle", i);
    }
    return data_[i];
}

Dihedral::Dihedral(size_t i, size_t j, size_t k, size_t m) {
    if (i == j || j == k || k == m) {
        throw error("can not have an atom linked to itself in a dihedral angle");
    }

    if (i == k || j == m || i == m) {
        throw error("can not have an atom twice in a dihedral angle");
    }

    if (std::max(i, j) < std::max(k, m)) {
        data_[0] = i;
        data_[1] = j;
        data_[2] = k;
        data_[3] = m;
    } else {
        data_[0] = m;
        data_[1] = k;
        data_[2] = j;
        data_[3] = i;
    }
}

size_t Dihedral::operator[](size_t i) const {
    if (i >= 4) {
        throw out_of_bounds("can not access atom n° {} in dihedral", i);
    }
    return data_[i];
}

Improper::Improper(size_t i, size_t j, size_t k, size_t m) {
    if (j == i || j == k || j == m) {
        throw error("can not have an atom linked to itself in an improper dihedral angle");
    }

    if (i == k || i == m || k == m) {
        throw error("can not have an atom twice in an improper dihedral angle");
    }

    std::array<size_t, 3> others = {{i, k, m}};
    std::sort(others.begin(), others.end());
    data_[0] = others[0];
    data_[1] = j;
    data_[2] = others[1];
    data_[3] = others[2];
}

size_t Improper::operator[](size_t i) const {
    if (i >= 4) {
        throw out_of_bounds("can not access atom n° {} in improper", i);
    }
    return data_[i];
}

void Connectivity::recalculate() const {
    angles_.clear();
    dihedrals_.clear();
    impropers_.clear();

    // Generate the list of which atom is bonded to which one
    auto bonded_to = std::vector<std::vector<size_t>>(biggest_atom_ + 1);
    for (auto const& bond: bonds_) {
        assert(bond[0] < bonded_to.size());
        assert(bond[1] < bonded_to.size());
        bonded_to[bond[0]].push_back(bond[1]);
        bonded_to[bond[1]].push_back(bond[0]);
    }

    // Generate the list of angles
    for (auto const& bond: bonds_) {
        auto i = bond[0];
        auto j = bond[1];
        for (auto k: bonded_to[i]) {
            if (k != j) {
                angles_.insert(Angle(k, i, j));
            }
        }

        for (auto k: bonded_to[j]) {
            if (k != i) {
                angles_.insert(Angle(i, j, k));
            }
        }
    }

    // Generate the list of dihedrals
    for (auto const& angle: angles_) {
        auto i = angle[0];
        auto j = angle[1];
        auto k = angle[2];
        for (auto m: bonded_to[i]) {
            if (m != j && m != k) {
                dihedrals_.insert(Dihedral(m, i, j, k));
            }
        }

        for (auto m: bonded_to[k]) {
            if (m != i && m != j) {
                dihedrals_.insert(Dihedral(i, j, k, m));
            }
        }

        for (auto m: bonded_to[j]) {
            if (m != i && m != k) {
                impropers_.insert(Improper(i, j, k, m));
            }
        }
    }

    uptodate_ = true;
}

const sorted_set<Bond>& Connectivity::bonds() const {
    if (!uptodate_) {
        recalculate();
    }
    return bonds_;
}

const sorted_set<Angle>& Connectivity::angles() const {
    if (!uptodate_) {
        recalculate();
    }
    return angles_;
}

const sorted_set<Dihedral>& Connectivity::dihedrals() const {
    if (!uptodate_) {
        recalculate();
    }
    return dihedrals_;
}

const sorted_set<Improper>& Connectivity::impropers() const {
    if (!uptodate_) {
        recalculate();
    }
    return impropers_;
}

void Connectivity::add_bond(size_t i, size_t j) {
    uptodate_ = false;
    bonds_.emplace(i, j);
    if (i > biggest_atom_) {biggest_atom_ = i;}
    if (j > biggest_atom_) {biggest_atom_ = j;}
}

void Connectivity::remove_bond(size_t i, size_t j) {
    auto pos = bonds_.find(Bond(i, j));
    if (pos != bonds_.end()) {
        uptodate_ = false;
        bonds_.erase(pos);
    }
}
