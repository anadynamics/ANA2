#ifndef ANA_MODES_H
#define ANA_MODES_H
#include <ANA/Includes.hpp>

using std::size_t;

namespace ANA {
namespace NDD {

    class Modes {
    public:
        Modes(std::string const &modes_filename);

        void get_modes_from_raw(std::string_view const texto);

        void calpha_to_full_atom(std::string const &pdb_filename);

        void normalize_evectors() {
            if (_evectors.size() != 0) {
                normalize_matrix(_evectors);
            } else {
                throw std::runtime_error(
                    "_evectors empty. Can't normalize them.");
            }
        }

        void normalize_atm_evectors() {
            if (_atm_evectors.size() != 0) {
                normalize_matrix(_atm_evectors);
            } else {
                throw std::runtime_error(
                    "_atm_evectors empty. Can't normalize them.");
            }
        }

        void normalize_matrix(std::vector<std::vector<double>> &mtx);

        std::vector<std::vector<double>> _atm_evectors;
        std::vector<std::vector<double>> _evectors;
        std::vector<double> _evals;
        size_t _i, _j, _iatm = 0;
    };
}
}
#endif // _H
