#ifndef ANA_MODES_H
#define ANA_MODES_H
#include <ANA/Includes.hpp>
#include <ANA/TextUtils.hpp>

using std::size_t;

namespace ANA {
namespace NDD {

    class Modes {
    public:
        Modes(std::string const &modes_filename, bool const amber_modes) :
            _amber_modes(amber_modes) {

            std::unique_ptr<char[]> const buffer = slurp(modes_filename);
            if (_amber_modes) {
                get_amber_modes_from_raw(std::string_view(buffer.get()));
            } else {
                ;
            }
            return;
        }

        Modes(std::string const &modes_filename,
            std::string const &pdb_filename, bool const amber_modes) :
            Modes(modes_filename, amber_modes) {
            calpha_to_full_atom(pdb_filename);
        }

        void get_amber_modes_from_raw(std::string_view const texto);

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
        std::vector<double> _normas_atm;
        size_t _i, _j, _iatm = 0;
        bool const _amber_modes;
    };
}
}
#endif // _H
