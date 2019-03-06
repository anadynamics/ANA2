// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstring>

#include "chemfiles/capi/selection.h"
#include "chemfiles/capi.hpp"

#include "chemfiles/Selections.hpp"
using namespace chemfiles;

static_assert(
    CHFL_MAX_SELECTION_SIZE == Match::MAX_MATCH_SIZE,
    "CHFL_MAX_SELECTION_SIZE should match Match::MAX_MATCH_SIZE"
);

extern "C" struct CAPISelection {
    CAPISelection(Selection&& select): selection(std::move(select)), matches() {}
    Selection selection;
    std::vector<Match> matches;
};

extern "C" CHFL_SELECTION* chfl_selection(const char* selection) {
    CHFL_SELECTION* c_selection = nullptr;
    CHFL_ERROR_GOTO(
        c_selection = new CAPISelection(Selection(std::string(selection)));
    )
    return c_selection;
error:
    delete c_selection;
    return nullptr;
}

extern "C" CHFL_SELECTION* chfl_selection_copy(const CHFL_SELECTION* const selection) {
    CHFL_SELECTION* new_selection = nullptr;
    CHFL_ERROR_GOTO(
        new_selection = new CAPISelection(Selection(selection->selection.string()));
    )
    return new_selection;
error:
    delete new_selection;
    return nullptr;
}

extern "C" chfl_status chfl_selection_size(const CHFL_SELECTION* const selection, uint64_t* size) {
    CHECK_POINTER(selection);
    CHFL_ERROR_CATCH(
        *size = selection->selection.size();
    )
}

extern "C" chfl_status chfl_selection_string(const CHFL_SELECTION* const selection, char* const string, uint64_t buffsize) {
    CHECK_POINTER(selection);
    CHECK_POINTER(string);
    CHFL_ERROR_CATCH(
        strncpy(string, selection->selection.string().c_str(), checked_cast(buffsize) - 1);
        string[buffsize - 1] = '\0';
    )
}

extern "C" chfl_status chfl_selection_evaluate(CHFL_SELECTION* const selection, const CHFL_FRAME* const frame, uint64_t* n_matches) {
    CHECK_POINTER(selection);
    CHFL_ERROR_CATCH(
        selection->matches = selection->selection.evaluate(*frame);
        *n_matches = selection->matches.size();
    )
}

extern "C" chfl_status chfl_selection_matches(const CHFL_SELECTION* const selection, chfl_match* const matches, uint64_t n_matches) {
    CHECK_POINTER(selection);
    if (n_matches != selection->matches.size()) {
        CAPI_LAST_ERROR = "Wrong data size in function 'chfl_selection_matches'.";
        return CHFL_MEMORY_ERROR;
    }
    CHFL_ERROR_CATCH(
        auto size = selection->selection.size();
        for (size_t i=0; i<n_matches; i++) {
            matches[i].size = size;
            for (size_t j=0; j<size; j++) {
                matches[i].atoms[j] = selection->matches[i][j];
            }

            for (uint64_t j=size; j<CHFL_MAX_SELECTION_SIZE; j++) {
                matches[i].atoms[j] = static_cast<uint64_t>(-1);
            }
        }
    )
}

extern "C" chfl_status chfl_selection_free(CHFL_SELECTION* const selection) {
    delete selection;
    return CHFL_SUCCESS;
}
