#include <ANA/Modes.hpp>

namespace ANA {
namespace NDD {

    // Create modes with eigenvalues.
    Modes::Modes(
        NDDOptions const &NDD_opts, std::string const &pdb_filename, AmberTag) {

        auto[bufr_modes, fsz_modes] = slurp(NDD_opts._modes_ndd_filename);
        get_amber_modes_from_raw(std::string_view(bufr_modes.get(), fsz_modes));

        // Now, if these are coarse grain modes, turn them into full atom ones.
        switch (NDD_opts._particles_per_residue) {
        case NDDOptions::FullAtom: {
            // Nothing to be done.
            break;
        }
        case NDDOptions::AlphaCarbon: {
            calpha_to_full_atom(pdb_filename);
            break;
        }
        case NDDOptions::Six: {
            six_to_full_atom(pdb_filename);
            break;
        }
        default: {
            throw std::runtime_error("Modes constructor error. "
                                     "NDD_opts._particles_per_residue in "
                                     "invalid state. This shouldn't happen.");
        }
        }
        return;
    }

    // // Create modes without eigenvalues.
    // Modes::Modes(std::string const &modes_filename,
    //     std::string const &pdb_filename, RowTag) {

    //     auto[bufr_modes, fsz_modes] = slurp(modes_filename);
    //     get_row_major_from_raw(std::string_view(bufr_modes.get(),
    //     fsz_modes)); calpha_to_full_atom(pdb_filename); return;
    // }

    // // Create modes without eigenvalues.
    // Modes::Modes(std::string const &modes_filename,
    //     std::string const &pdb_filename, ColumnTag) {

    //     auto[bufr_modes, fsz_modes] = slurp(modes_filename);
    //     get_col_major_from_raw(std::string_view(bufr_modes.get(),
    //     fsz_modes)); calpha_to_full_atom(pdb_filename); return;
    // }

    // Create modes reading from a Row Major file.
    Modes::Modes(
        NDDOptions const &NDD_opts, std::string const &pdb_filename, RowTag) {
        // Get modes
        auto[bufr_modes, fsz_modes] = slurp(NDD_opts._modes_ndd_filename);
        get_row_major_from_raw(std::string_view(bufr_modes.get(), fsz_modes));

        // Get frequencies (eigenvalues) if available.
        if (NDD_opts._freqs_ndd_filename != "none") {
            auto[bufr_freqs, fsz_evals] = slurp(NDD_opts._freqs_ndd_filename);
            _evalues = get_values_from_raw(
                std::string_view(bufr_freqs.get(), fsz_evals));
            if (_evalues.size() != _j) {
                std::cerr << "Vector count: " << _j
                          << ". Scalar count: " << _evalues.size() << '\n';
                throw std::runtime_error(
                    "Frequencies don't match vectors. Aborting.");
            }
        }

        // Now, if these are coarse grain modes, turn them into full atom ones.
        switch (NDD_opts._particles_per_residue) {
        case NDDOptions::FullAtom: {
            // Nothing to be done.
            break;
        }
        case NDDOptions::AlphaCarbon: {
            calpha_to_full_atom(pdb_filename);
            break;
        }
        case NDDOptions::Six: {
            six_to_full_atom(pdb_filename);
            break;
        }
        default: {
            throw std::runtime_error("Modes constructor error. "
                                     "NDD_opts._particles_per_residue in "
                                     "invalid state. This shouldn't happen.");
        }
        }
        return;
    }

    // Create modes reading from a Column Major file.
    Modes::Modes(NDDOptions const &NDD_opts, std::string const &pdb_filename,
        ColumnTag) {
        // Get modes
        auto[bufr_modes, fsz_modes] = slurp(NDD_opts._modes_ndd_filename);
        get_col_major_from_raw(std::string_view(bufr_modes.get(), fsz_modes));

        // Get frequencies (eigenvalues) if available.
        if (NDD_opts._freqs_ndd_filename != "none") {
            auto[bufr_freqs, fsz_evals] = slurp(NDD_opts._freqs_ndd_filename);
            _evalues = get_values_from_raw(
                std::string_view(bufr_freqs.get(), fsz_evals));

            if (_evalues.size() != _j) {
                std::cerr << "Vector count: " << _j
                          << ". Scalar count: " << _evalues.size() << '\n';
                throw std::runtime_error(
                    "Frequencies don't match vectors. Aborting.");
            }
        }
        // Now, if these are coarse grain modes, turn them into full atom ones.
        switch (NDD_opts._particles_per_residue) {
        case NDDOptions::FullAtom: {
            // Nothing to be done.
            break;
        }
        case NDDOptions::AlphaCarbon: {
            calpha_to_full_atom(pdb_filename);
            break;
        }
        case NDDOptions::Six: {
            six_to_full_atom(pdb_filename);
            break;
        }
        default: {
            throw std::runtime_error("Modes constructor error. "
                                     "NDD_opts._particles_per_residue in "
                                     "invalid state. This shouldn't happen.");
        }
        }
        write_matrix(_evectors, _i, _j, "evectors");
        write_matrix(_atm_evectors, _ii, _j, "atm_evectors");
        std::cout << _i << " " << _j << "  " << _ii << "  " << _natoms << '\n';

        return;
    }

    // // Create full atom modes.
    // Modes::Modes(std::string const &modes_filename, AmberTag,
    //     std::string const &pdb_filename) :
    //     Modes(modes_filename, AmberTag()) {

    // }

    // // Create full atom modes.
    // Modes::Modes(std::string const &modes_filename, RowTag,
    //     std::string const &pdb_filename) :
    //     Modes(modes_filename, RowTag()) {

    //     calpha_to_full_atom(pdb_filename);
    // }

    // // Create full atom modes.
    // Modes::Modes(std::string const &modes_filename, ColumnTag,
    //     std::string const &pdb_filename) :
    //     Modes(modes_filename, ColumnTag()) {

    //     calpha_to_full_atom(pdb_filename);
    // }

    void Modes::get_amber_modes_from_raw(std::string_view const texto) {

        size_t beg = 57;
        if (isspace(texto[38])) {
            std::string_view st_j(texto.data() + 34, 3);
            _j = std::stoi(st_j.data());
            std::string_view st_i(texto.data() + 49, 3);
            _i = std::stoi(st_i.data());
        } else {
            std::string_view st_j(texto.data() + 34, 4);
            _j = std::stoi(st_j.data());
            std::string_view st_i(texto.data() + 48, 4);
            _i = std::stoi(st_i.data());
            ++beg;
        }
        _evectors.reserve(_j);
        _evalues.reserve(_j);
        // representa: " ****"
        size_t constexpr spacer = 6;
        size_t constexpr ch_elem = 11;
        size_t constexpr ch_num = 5;
        size_t constexpr ch_eval = 12;
        size_t constexpr ncolumns = 7;
        size_t const resto = _i % ncolumns;
        size_t const nlines = (resto) ? _i / ncolumns + 1 : _i / ncolumns;
        // número de caracteres q ocupa 1 modo en el formato de Amber.
        size_t const ch_mode = (nlines - 1) * ncolumns * ch_elem +
            resto * ch_elem + nlines + spacer;

        const char *cursor = texto.data() + beg + ch_mode;
        for (size_t K = 0; K < _j; ++K) {
            char const *start_eval = cursor + ch_num;
            char *end_eval = const_cast<char *>(start_eval) + ch_eval;
            _evalues.push_back(std::strtof(start_eval, &end_eval));

            auto idx = cursor + ch_num + ch_eval;
            size_t nl = 0;
            std::vector<double> temp;
            temp.reserve(_i);
            // nl is the new line counter. At the beginning, idx is pointing
            // to the new line character after the eigenvalue and before the
            // first element of the eigenvector, and since 0 % 7 evalueates
            // to false, nl starts is incremented in the first iteration.
            for (size_t k = 0; k < _i; ++k) {
                (k % 7) ? nl = 0 : nl = 1;
                char const *start_elem = idx;
                char *end_elem = const_cast<char *>(start_elem) + ch_elem;
                temp.push_back(std::strtof(start_elem, &end_elem));
                idx = idx + ch_elem + nl;
            }
            _evectors.push_back(temp);
            cursor = cursor + ch_mode + ch_num + ch_eval + 1;
        }

        return;
    }

    void Modes::get_row_major_from_raw(std::string_view const texto) {
        auto[ch_elem, ncols, nrows] = guess_format(texto);
        _j = nrows;
        _i = ncols;
        int line_length = ncols * ch_elem + 1;
        char const *cursor = texto.data();

        _evectors.resize(_j);
        for (size_t j = 0; j < _j; ++j) {
            std::vector<double> vector;
            vector.reserve(_i);

            for (size_t i = 0; i < _i; ++i) {
                char const *left{cursor + i * ch_elem};
                char *right{const_cast<char *>(left + ch_elem)};

                vector.push_back(std::strtof(left, &right));
            }

            _evectors.push_back(std::move(vector));
            cursor += line_length;
        }

        return;
    }

    void Modes::get_col_major_from_raw(std::string_view const texto) {
        auto[ch_elem, ncols, nrows] = guess_format(texto);
        _j = ncols;
        _i = nrows;
        int line_length = ncols * ch_elem + 1;
        char const *cursor = texto.data();

        _evectors.resize(_j);
        for (size_t i = 0; i < _i; ++i) {
            for (size_t j = 0; j < _j; ++j) {
                char const *left{cursor + j * ch_elem};
                char *right{const_cast<char *>(left + ch_elem)};

                _evectors[j].push_back(std::strtof(left, &right));
            }
            cursor += line_length;
        }

        return;
    }

    void Modes::normalize_matrix(std::vector<std::vector<double>> &mtx) {
        if (_normas.size() == _j) {
            // Vector norms are already calculated.
            for (auto &vector : mtx) {
                for (size_t ii = 0; ii != _i; ++ii) {
                    vector[ii] = vector[ii] / _normas[ii];
                }
            }
        } else {
            for (auto &vector : mtx) {
                double sum = 0.;
                for (auto const elmt : vector) {
                    sum += elmt * elmt;
                }
                double const norm = sqrt(sum);
                for (auto &elmt : vector) {
                    elmt = elmt / norm;
                }
            }
        }

        return;
    }

    void Modes::calpha_to_full_atom(std::string const &pdb_filename) {

        chemfiles::Trajectory in_trj(pdb_filename);
        chemfiles::Frame const in_frm = in_trj.read();
        auto const in_top = in_frm.topology();
        auto const res = in_top.residues();
        int const nres = res.size();
        assert(static_cast<size_t>(nres * 3) == _i &&
            static_cast<size_t>(nres * 3 - 6) == _j);

        // Get number of atoms per each residue.
        std::vector<int> atoms_per_res;
        atoms_per_res.reserve(nres);
        for (auto const &residue : res) {
            int const natoms = residue.size();
            atoms_per_res.push_back(natoms);
            _natoms += natoms;
        }

        // Length of full atom eigenvectors.
        _ii = _natoms * 3;

        FILE *natoms_file = std::fopen("natoms", "w");
        if (natoms_file) {
            for (auto const &residue : atoms_per_res) {
                fmt::print(natoms_file, "{}\n", residue);
            }
        }
        std::fclose(natoms_file);

        // Go eigenvector by eigenvector and repeat each of its elements
        // according to atoms_per_res. Store that into _atm_evectors.
        // Also store each vector's norm.
        _atm_evectors.reserve(_j);
        for (size_t kk = 0; kk < _j; ++kk) {
            std::vector<double> atm_evectors;
            atm_evectors.reserve(_ii);
            int res_number = 0;

            _normas.reserve(_j);
            double sum = 0.;

            for (int const natoms : atoms_per_res) {
                for (int k = 0; k < natoms; ++k) {
                    double const elmt_x = _evectors[kk][3 * res_number];
                    double const elmt_y = _evectors[kk][3 * res_number + 1];
                    double const elmt_z = _evectors[kk][3 * res_number + 2];

                    atm_evectors.push_back(elmt_x);
                    atm_evectors.push_back(elmt_y);
                    atm_evectors.push_back(elmt_z);
                    sum += elmt_x * elmt_x + elmt_y * elmt_y + elmt_z * elmt_z;
                }
                ++res_number;
            }

            _atm_evectors.push_back(std::move(atm_evectors));
            _normas.push_back(sqrt(sum));
        }
        // normalize_matrix(_atm_evectors);
        return;
    }

    void Modes::six_to_full_atom(std::string const &pdb_filename) {

        chemfiles::Trajectory in_trj(pdb_filename);
        chemfiles::Frame const in_frm = in_trj.read();
        auto const in_top = in_frm.topology();
        auto const res = in_top.residues();
        int const nres = res.size();
        assert(static_cast<size_t>(nres * 3) == _i &&
            static_cast<size_t>(nres * 3 - 6) == _j);

        // Get number of atoms per each residue.
        std::vector<int> atoms_per_res;
        atoms_per_res.reserve(nres);
        for (auto const &residue : res) {
            int const natoms = residue.size();
            atoms_per_res.push_back(natoms);
            _natoms += natoms;
        }

        // Length of full atom eigenvectors.
        _ii = _natoms * 3;

        // Go eigenvector by eigenvector and repeat each of its elements
        // according to atoms_per_res. Store that into _atm_evectors.
        // Also store each vector's norm.
        _atm_evectors.reserve(_j);
        for (size_t kk = 0; kk < _j; ++kk) {
            std::vector<double> atm_evectors;
            atm_evectors.reserve(_ii);
            int res_number = 0;

            _normas.reserve(_j);
            double sum = 0.;

            for (int const natoms : atoms_per_res) {
                for (int k = 0; k < natoms; ++k) {
                    double const elmt = _evectors[kk][res_number];
                    atm_evectors.push_back(elmt);
                    sum += elmt * elmt;
                }
                ++res_number;
            }

            _atm_evectors.push_back(std::move(atm_evectors));
            _normas.push_back(sqrt(sum));
        }
        normalize_matrix(_atm_evectors);
        return;
    }

    // Call the proper Modes constructor.
    Modes create_modes(
        NDDOptions const &NDD_opts, std::string const &pdb_filename) {

        if (NDD_opts._modes_format == "amber") {
            return {NDD_opts, pdb_filename, AmberTag()};

        } else if (NDD_opts._modes_format == "row") {
            return {NDD_opts, pdb_filename, RowTag()};

        } else if (NDD_opts._modes_format == "column") {
            return {NDD_opts, pdb_filename, ColumnTag()};
        }
        throw std::invalid_argument(
            "No Modes format input could be parsed. This shouldn't happen.");
        return {};
    }
}
}