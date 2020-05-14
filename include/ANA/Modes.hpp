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
        Modes(NDDOptions const &NDD_opts, std::string const &pdb_filename,
            AmberTag);

        // Create modes reading from a Row Major file.
        Modes(NDDOptions const &NDD_opts, std::string const &pdb_filename,
            RowTag);
        // Create modes reading from a Column Major file.
        Modes(NDDOptions const &NDD_opts, std::string const &pdb_filename,
            ColumnTag);

        void get_amber_modes_from_raw(std::string_view const texto);

        void get_row_major_from_raw(std::string_view const texto);

        void get_col_major_from_raw(std::string_view const texto);

        void calpha_to_full_atom(std::string const &pdb_filename);

        void six_to_full_atom(std::string const &pdb_filename);

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

        auto six_to_full_atom_helper(int resi, std::string const &resn,
            std::vector<std::string> const &atm_names, size_t mode, int &idx,
            int &atm_cnt, std::vector<double> &atm_evector,
            std::vector<double>::iterator beg) -> std::vector<double>::iterator;

        auto get_res_info(chemfiles::Topology const &in_top)
            -> std::tuple<int, std::vector<std::string>, std::vector<int>,
                std::vector<std::string>>;

        std::vector<std::vector<double>> _atm_evectors;
        std::vector<std::vector<double>> _evectors;
        std::vector<double> _evalues;
        std::vector<double> _normas;
        // _i: length of original eigenvectors, _j: nbr of eigenvectors, _ii:
        // length of full atom eigenvectors. _ii = _natoms * 3.
        size_t _i, _j, _ii, _natoms = 0;
    };

    // Call the proper Modes constructor.
    Modes create_modes(
        NDDOptions const &NDD_opts, std::string const &pdb_filename);
}
}
#endif // _H
