// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "chemfiles/capi/cell.h"
#include "chemfiles/capi.hpp"

#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Frame.hpp"
using namespace chemfiles;

static_assert(sizeof(chfl_cellshape) == sizeof(int), "Wrong size for chfl_cellshape enum");

extern "C" CHFL_CELL* chfl_cell(const chfl_vector3d lenghts) {
    CHFL_CELL* cell = nullptr;
    CHECK_POINTER_GOTO(lenghts);
    CHFL_ERROR_GOTO(
        cell = new UnitCell(lenghts[0], lenghts[1], lenghts[2]);
    )
    return cell;
error:
    delete cell;
    return nullptr;
}

extern "C" CHFL_CELL* chfl_cell_triclinic(const chfl_vector3d lenghts, const chfl_vector3d angles) {
    CHFL_CELL* cell = nullptr;
    CHECK_POINTER_GOTO(lenghts);
    CHECK_POINTER_GOTO(angles);
    CHFL_ERROR_GOTO(
        cell = new UnitCell(
            lenghts[0], lenghts[1], lenghts[2],
            angles[0], angles[1], angles[2]
        );
        // ensure that the unit cell shape is always TRICLINIC, even if the
        // three angles are 90°.
        cell->shape(UnitCell::TRICLINIC);
    )
    return cell;
error:
    delete cell;
    return nullptr;
}


extern "C" CHFL_CELL* chfl_cell_from_frame(const CHFL_FRAME* const frame) {
    CHFL_CELL* cell = nullptr;
    CHECK_POINTER_GOTO(frame);
    CHFL_ERROR_GOTO(
        cell = new UnitCell(frame->cell());
    )
    return cell;
error:
    delete cell;
    return nullptr;
}

extern "C" CHFL_CELL* chfl_cell_copy(const CHFL_CELL* const cell) {
    CHFL_CELL* new_cell = nullptr;
    CHFL_ERROR_GOTO(
        new_cell = new UnitCell(*cell);
    )
    return new_cell;
error:
    delete new_cell;
    return nullptr;
}

extern "C" chfl_status chfl_cell_volume(const CHFL_CELL* const cell, double* volume) {
    CHECK_POINTER(cell);
    CHECK_POINTER(volume);
    CHFL_ERROR_CATCH(
        *volume = cell->volume();
    )
}

extern "C" chfl_status chfl_cell_lengths(const CHFL_CELL* const cell, chfl_vector3d lenghts) {
    CHECK_POINTER(cell);
    CHECK_POINTER(lenghts);
    CHFL_ERROR_CATCH(
        lenghts[0] = cell->a();
        lenghts[1] = cell->b();
        lenghts[2] = cell->c();
    )
}

extern "C" chfl_status chfl_cell_set_lengths(CHFL_CELL* const cell, const chfl_vector3d lenghts) {
    CHECK_POINTER(cell);
    CHECK_POINTER(lenghts);
    CHFL_ERROR_CATCH(
        cell->set_a(lenghts[0]);
        cell->set_b(lenghts[1]);
        cell->set_c(lenghts[2]);
    )
}

extern "C" chfl_status chfl_cell_angles(const CHFL_CELL* const cell, chfl_vector3d angles) {
    CHECK_POINTER(cell);
    CHECK_POINTER(angles);
    CHFL_ERROR_CATCH(
        angles[0] = cell->alpha();
        angles[1] = cell->beta();
        angles[2] = cell->gamma();
    )
}

extern "C" chfl_status chfl_cell_set_angles(CHFL_CELL* const cell, const chfl_vector3d angles) {
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        cell->set_alpha(angles[0]);
        cell->set_beta(angles[1]);
        cell->set_gamma(angles[2]);
    )
}

extern "C" chfl_status chfl_cell_matrix(const CHFL_CELL* const cell, chfl_vector3d matrix[3]) {
    CHECK_POINTER(cell);
    CHECK_POINTER(matrix);
    CHFL_ERROR_CATCH(
        cell->raw_matricial(matrix);
    )
}

extern "C" chfl_status chfl_cell_shape(const CHFL_CELL* const cell, chfl_cellshape* const shape) {
    CHECK_POINTER(cell);
    CHECK_POINTER(shape);
    CHFL_ERROR_CATCH(
        *shape = static_cast<chfl_cellshape>(cell->shape());
    )
}

extern "C" chfl_status chfl_cell_set_shape(CHFL_CELL* const cell, chfl_cellshape shape) {
    CHECK_POINTER(cell);
    CHFL_ERROR_CATCH(
        cell->shape(static_cast<UnitCell::CellShape>(shape));
    )
}

extern "C" chfl_status chfl_cell_wrap(const CHFL_CELL* const cell, chfl_vector3d vector) {
    CHECK_POINTER(cell);
    CHECK_POINTER(vector);
    CHFL_ERROR_CATCH(
        auto result = cell->wrap(vector3d(vector));
        vector[0] = result[0];
        vector[1] = result[1];
        vector[2] = result[2];
    )
}

extern "C" chfl_status chfl_cell_free(CHFL_CELL* const cell) {
    delete cell;
    return CHFL_SUCCESS;
}
