// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <cstring>

#include "chemfiles/capi/atom.h"
#include "chemfiles/capi.hpp"

#include "chemfiles/ErrorFmt.hpp"
#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Topology.hpp"
using namespace chemfiles;

extern "C" CHFL_ATOM* chfl_atom(const char* name) {
    CHFL_ATOM* atom = nullptr;
    CHFL_ERROR_GOTO(
        atom = new Atom(name);
    )
    return atom;
error:
    delete atom;
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_copy(const CHFL_ATOM* const atom) {
    CHFL_ATOM* new_atom = nullptr;
    CHFL_ERROR_GOTO(
        new_atom = new Atom(*atom);
    )
    return new_atom;
error:
    delete new_atom;
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_from_frame(const CHFL_FRAME* const frame, uint64_t idx) {
    CHFL_ATOM* atom = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHFL_ERROR_GOTO(
        // Return NULL if the index is out of bounds
        if (idx >= frame->natoms()) {
            throw out_of_bounds(
                "out of bounds atomic index in `chfl_atom_from_frame`: we have {} atoms, but the index is {}",
                frame->natoms(), idx
            );
        }
        atom = new Atom(frame->topology()[checked_cast(idx)]);
    )
    return atom;
error:
    delete atom;
    return nullptr;
}

extern "C" CHFL_ATOM* chfl_atom_from_topology(const CHFL_TOPOLOGY* const topology, uint64_t idx) {
    CHFL_ATOM* atom = nullptr;
    CHECK_POINTER_GOTO(topology);
    CHFL_ERROR_GOTO(
        // Return NULL if the index is out of bounds
        if (idx >= topology->natoms()) {
            throw out_of_bounds(
                "out of bounds atomic index in `chfl_atom_from_topology`: we have {} atoms, but the index is {}",
                topology->natoms(), idx
            );
        }
        atom = new Atom((*topology)[checked_cast(idx)]);
    )
    return atom;
error:
    delete atom;
    return nullptr;
}

extern "C" chfl_status chfl_atom_mass(const CHFL_ATOM* const atom, double* mass) {
    CHECK_POINTER(atom);
    CHECK_POINTER(mass);
    CHFL_ERROR_CATCH(
        *mass = atom->mass();
    )
}

extern "C" chfl_status chfl_atom_set_mass(CHFL_ATOM* const atom, double mass) {
    CHECK_POINTER(atom);
    CHFL_ERROR_CATCH(
        atom->set_mass(mass);
    )
}

extern "C" chfl_status chfl_atom_charge(const CHFL_ATOM* const atom, double* charge) {
    CHECK_POINTER(atom);
    CHECK_POINTER(charge);
    CHFL_ERROR_CATCH(
        *charge = atom->charge();
    )
}

extern "C" chfl_status chfl_atom_set_charge(CHFL_ATOM* const atom, double charge) {
    CHECK_POINTER(atom);
    CHFL_ERROR_CATCH(
        atom->set_charge(charge);
    )
}

extern "C" chfl_status chfl_atom_type(const CHFL_ATOM* const atom, char* const type, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(type);
    CHFL_ERROR_CATCH(
        strncpy(type, atom->type().c_str(), checked_cast(buffsize) - 1);
        type[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_atom_set_type(CHFL_ATOM* const atom, const char* type) {
    CHECK_POINTER(atom);
    CHECK_POINTER(type);
    CHFL_ERROR_CATCH(
        atom->set_type(type);
    )
}

extern "C" chfl_status chfl_atom_name(const CHFL_ATOM* const atom, char* const name, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        strncpy(name, atom->name().c_str(), checked_cast(buffsize) - 1);
        name[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_atom_set_name(CHFL_ATOM* const atom, const char* name) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        atom->set_name(name);
    )
}

extern "C" chfl_status chfl_atom_full_name(const CHFL_ATOM* const atom, char* const name, uint64_t buffsize) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHFL_ERROR_CATCH(
        auto full_name = atom->full_name();
        if (full_name) {
            std::strncpy(name, full_name.value().c_str(), checked_cast(buffsize) - 1);
            name[buffsize - 1] = '\0';
        } else {
            std::memset(name, 0, checked_cast(buffsize));
        }
    )
}

extern "C" chfl_status chfl_atom_vdw_radius(const CHFL_ATOM* const atom, double* radius) {
    CHECK_POINTER(atom);
    CHECK_POINTER(radius);
    CHFL_ERROR_CATCH(
        *radius = atom->vdw_radius().value_or(0);
    )
}

extern "C" chfl_status chfl_atom_covalent_radius(const CHFL_ATOM* const atom, double* radius) {
    CHECK_POINTER(atom);
    CHECK_POINTER(radius);
    CHFL_ERROR_CATCH(
        *radius = atom->covalent_radius().value_or(0);
    )
}

extern "C" chfl_status chfl_atom_atomic_number(const CHFL_ATOM* const atom, uint64_t* number) {
    CHECK_POINTER(atom);
    CHECK_POINTER(number);
    CHFL_ERROR_CATCH(
        *number = atom->atomic_number().value_or(0ul);
    )
}

extern "C" chfl_status chfl_atom_set_property(CHFL_ATOM* const atom, const char* name, const CHFL_PROPERTY* const property) {
    CHECK_POINTER(atom);
    CHECK_POINTER(name);
    CHECK_POINTER(property);
    CHFL_ERROR_CATCH(
        atom->set(name, *property);
    )
}

extern "C" CHFL_PROPERTY* chfl_atom_get_property(const CHFL_ATOM* const atom, const char* name) {
    CHFL_PROPERTY* property = nullptr;
    CHECK_POINTER_GOTO(atom);
    CHECK_POINTER_GOTO(name);
    CHFL_ERROR_GOTO(
        auto atom_property = atom->get(name);
        if (atom_property) {
            property = new Property(*atom_property);
        } else {
            throw property_error("can not find a property named {} in this atom", name);
        }
    )
    return property;
error:
    delete property;
    return nullptr;
}

extern "C" chfl_status chfl_atom_free(CHFL_ATOM* const atom) {
    delete atom;
    return CHFL_SUCCESS;
}
