#include <ANA/Modes.hpp>

namespace ANA {
namespace NDD {

    Modes::Modes(std::string const &filename) {

        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            ifs.seekg(0, ifs.end);
            size_t const fsz = static_cast<size_t>(ifs.tellg());
            std::unique_ptr<char[]> const buffer =
                std::make_unique<char[]>(fsz);
            ifs.seekg(0);
            ifs.read(buffer.get(), fsz);
            ifs.close();

            get_modes_from_raw(std::string_view(buffer.get(), fsz));
        } else
            std::cout << "Could not get vectors from: " << filename << '\n';
    }

    void Modes::get_modes_from_raw(std::string_view const texto) {

        size_t beg = 57;
        if (isspace(texto[38])) {
            std::string_view st_j(texto.data() + 34, 3);
            j = std::stoi(st_j.data());
            std::string_view st_i(texto.data() + 49, 3);
            i = std::stoi(st_i.data());
        } else {
            std::string_view st_j(texto.data() + 34, 4);
            j = std::stoi(st_j.data());
            std::string_view st_i(texto.data() + 48, 4);
            i = std::stoi(st_i.data());
            ++beg;
        }
        evectors.reserve(j);
        evals.reserve(j);
        size_t constexpr ch_elem = 11;
        size_t constexpr ch_num = 5;
        size_t constexpr ch_eval = 12;
        size_t constexpr ncolumns = 7;
        size_t const resto = i % ncolumns;
        size_t const nlines = (resto) ? i / ncolumns + 1 : i / ncolumns;
        // representa: " ****"
        size_t constexpr spacer = 6;
        // número de caracteres q ocupa 1 modo en el formato de Amber.
        size_t const ch_mode = (nlines - 1) * ncolumns * ch_elem +
            resto * ch_elem + nlines + spacer;

        const char *cursor = texto.data() + beg + ch_mode;
        for (size_t K = 0; K < j; ++K) {
            char const *start_eval = cursor + ch_num;
            char *end_eval = const_cast<char *>(start_eval) + ch_eval;
            evals.push_back(std::strtof(start_eval, &end_eval));

            auto idx = cursor + ch_num + ch_eval;
            size_t nl = 0;
            std::vector<double> temp;
            temp.reserve(i);
            // nl is the new line counter. At the beginning, idx is pointing
            // to the new line character after the eigenvalue and before the
            // first element of the eigenvector, and since 0 % 7 evalueates
            // to false, nl starts is incremented in the first iteration.
            for (size_t k = 0; k < i; ++k) {
                (k % 7) ? nl = 0 : nl = 1;
                char const *start_elem = idx;
                char *end_elem = const_cast<char *>(start_elem) + ch_elem;
                temp.push_back(std::strtof(start_elem, &end_elem));
                idx = idx + ch_elem + nl;
            }
            evectors.push_back(std::move(temp));
            cursor = cursor + ch_mode + ch_num + ch_eval + 1;
        }

        return;
    }
}
}