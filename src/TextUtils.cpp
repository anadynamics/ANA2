#include <ANA/TextUtils.hpp>
namespace ANA {
namespace NDD {
    void write_vector(std::vector<double> vec, std::string const &filename) {
        FILE *out_file = std::fopen(filename.c_str(), "w");
        if (out_file) {
            for (auto const &each : vec) {
                fmt::print(out_file, "{:8.3f}\n", each);
            }
        } else {
            printf("Could not write NDD output to: %s.\n", filename.c_str());
        }
        return;
    }
}

void write_vector(std::vector<double> vec, std::string const &filename) {
    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        for (size_t i = 0; i < vec.size(); ++i) {
            fmt::print(out_file, "{: <5}{: >5}{:8.3f}\n", "Frame", i, vec[i]);
        }
    } else {
        printf("Could not write NDD output to: %s.\n", filename.c_str());
    }
    return;
}

void write_matrix(std::vector<std::vector<double>> mtx, size_t nrows,
    size_t ncols, std::string const &filename) {
    FILE *out_file = std::fopen(filename.c_str(), "w");

    if (mtx.size() < ncols || mtx[0].size() < nrows) {
        throw std::runtime_error("write_matrix(): Matrix dimensions are "
                                 "smaller than nrows and ncols specified.\n");
    }

    if (out_file) {
        for (size_t i = 0; i != nrows; ++i) {
            for (size_t j = 0; j != ncols; ++j) {
                fmt::print(out_file, " {:9.6f}", mtx[j][i]);
            }
            fmt::print(out_file, "\n");
        }
    } else {
        printf("Could not write NDD output to: %s.\n", filename.c_str());
    }
    return;
}

// Read whole file into memory and return unique pointer to it and its size.
auto slurp(std::string const &filename)
    -> std::tuple<std::unique_ptr<char[]>, size_t> {

    std::ifstream ifs(filename);
    if (ifs.is_open()) {
        ifs.seekg(0, ifs.end);
        size_t const fsz = static_cast<size_t>(ifs.tellg());
        std::unique_ptr<char[]> buffer = std::make_unique<char[]>(fsz);
        ifs.seekg(0);
        ifs.read(buffer.get(), fsz);
        ifs.close();
        return {std::move(buffer), fsz};
    } else {
        throw std::runtime_error("Could not read Amber vectors.");
    }
}

// Returns number of chars in each element (including whitespaces), nbr of
// elements in line, nbr of lines.
auto guess_format(std::string_view texto)
    -> std::tuple<size_t, size_t, size_t> {

    size_t ch_elem = 0, line_length = 0;
    size_t const fsz = texto.size();
    bool reading_element = false;

    for (size_t c = 0; c < fsz; ++c) {
        if (iscntrl(texto[c])) {
            // A new line was found. Done counting characters per line.
            line_length = c;
            if (ch_elem == 0) {
                // newline char read before a space char. There is only 1
                // element per line.
                ch_elem = line_length;
            }
            break;
        } else if (isspace(texto[c])) {
            if (!reading_element) {
                // Leading whitespace.
                continue;
            } else if (ch_elem == 0) {
                // A whitespace was found after reading an element. Done
                // counting characters per element.
                ch_elem = c;
            } else {
                // Done reading another element. Go on until the whole line is
                // read.
                continue;
            }
        } else {
            // This character is neither whitespace nor a newline.
            reading_element = true;
        }
    }

    bool const malformed_element = (line_length % ch_elem != 0);
    bool const malformed_line = (fsz - 1) % (line_length + 1) != 0;
    bool const malformed_line_last_newline = fsz % (line_length + 1) != 0;

    if (malformed_element || (malformed_line && malformed_line_last_newline)) {
        std::cerr << "File size: " << fsz << " Line length: " << line_length
                  << " Element length: " << ch_elem << '\n';
        throw std::runtime_error(
            "guess_format(): corrupted input. Line length is not a "
            "multiple of element length.\nRemember to use fixed-width "
            "format.");
    }

    size_t elems_per_line{line_length / ch_elem};
    size_t lines_per_file{(fsz - 1) / (line_length + 1)};
    if (malformed_line) {
        // There's a final newline at the end of the file.
        lines_per_file = fsz / (line_length + 1);
    }

    return {ch_elem, elems_per_line, lines_per_file};
}

auto get_values_from_raw(std::string_view const texto) -> std::vector<double> {
    auto[ch_elem, ncols, nrows] = guess_format(texto);
    int line_length = ncols * ch_elem + 1;
    char const *cursor = texto.data();
    std::vector<double> data;

    data.reserve(nrows);
    for (size_t j = 0; j < nrows; ++j) {
        char const *left{cursor};
        char *right{const_cast<char *>(left + ch_elem)};
        data.push_back(std::strtof(left, &right));
        cursor += line_length;
    }

    return data;
}

} // namespace ANA