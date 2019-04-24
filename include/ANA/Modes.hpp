#ifndef ANA_MODES_H
#define ANA_MODES_H
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>
#include <ANA/TextUtils.hpp>

using std::size_t;

namespace ANA {
namespace NDD {

    struct AmberTag {};
    struct ColumnTag {};
    struct RowTag {};

    class Modes {
    public:
        Modes() = default;

        // Create modes with eigenvalues.
        Modes(std::string const &modes_filename, AmberTag);

        // Create modes without eigenvalues.
        Modes(std::string const &modes_filename, RowTag);
        // Create modes without eigenvalues.
        Modes(std::string const &modes_filename, ColumnTag);

        // Create modes with eigenvalues.
        Modes(std::string const &modes_filename,
            std::string const &evals_filename, RowTag);
        // Create modes with eigenvalues.
        Modes(std::string const &modes_filename,
            std::string const &evals_filename, ColumnTag);

        // Create full atom modes.
        Modes(std::string const &modes_filename, AmberTag,
            std::string const &pdb_filename);
        // Create full atom modes.
        Modes(std::string const &modes_filename, RowTag,
            std::string const &pdb_filename);
        // Create full atom modes.
        Modes(std::string const &modes_filename, ColumnTag,
            std::string const &pdb_filename);

        void get_amber_modes_from_raw(std::string_view const texto);

        // TODO: get_rmajor_from_raw

        void get_cmajor_from_raw(std::string_view const texto);

        void get_evals_from_raw(std::string_view const texto);

        void read_evals(std::string const &evals_filename) {

            std::unique_ptr<char[]> const buffer_evals = slurp(evals_filename);
            get_evals_from_raw(std::string_view(buffer_evals.get()));
            return;
        }

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
    };

    // Call the proper Modes constructor.
    Modes create_modes(NDDOptions const &NDD_opts);
}
}
#endif // _H
