#ifndef TEXT_UTILS_H
#define TEXT_UTILS_H
#include <ANA/Includes.hpp>

namespace ANA {
namespace NDD {
    void write_vector(std::vector<double> vec, std::string const &filename);
}

void write_vector(std::vector<double> vec, std::string const &filename);

void write_matrix(std::vector<std::vector<double>> mtx, size_t nrows,
    size_t ncols, std::string const &filename);

// Read whole file into memory and return unique pointer to it and its size.
auto slurp(std::string const &filename)
    -> std::tuple<std::unique_ptr<char[]>, size_t>;

// Returns number of chars in each element, nbr of elements in line, nbr of
// lines.
auto guess_format(std::string_view texto) -> std::tuple<size_t, size_t, size_t>;

auto get_values_from_raw(std::string_view const texto) -> std::vector<double>;
}
#endif // _H