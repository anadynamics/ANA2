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

} // namespace ANA