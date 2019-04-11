#ifndef TEXT_UTILS_H
#define TEXT_UTILS_H
#include <ANA/Includes.hpp>

namespace ANA {
namespace NDD {
    void write_vector(std::vector<double> vec, std::string const &filename);
}

void write_vector(std::vector<double> vec, std::string const &filename);

// Read whole file into memory and return unique pointer to it and its size.
auto slurp(std::string const &filename) -> std::unique_ptr<char[]>;

//
auto guess_format(std::string_view texto) -> std::tuple<size_t, size_t, size_t>;
}
#endif // _H