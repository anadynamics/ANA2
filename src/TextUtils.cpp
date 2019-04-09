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

// Read whole file into memory and return unique pointer to it and its size.
auto slurp(std::string const &filename) -> std::unique_ptr<char[]> {

    std::ifstream ifs(filename);
    if (ifs.is_open()) {
        ifs.seekg(0, ifs.end);
        size_t const fsz = static_cast<size_t>(ifs.tellg());
        std::unique_ptr<char[]> buffer = std::make_unique<char[]>(fsz);
        ifs.seekg(0);
        ifs.read(buffer.get(), fsz);
        ifs.close();
        return std::move(buffer);
    } else {
        throw std::runtime_error("Could not read Amber vectors.");
    }
}

} // namespace ANA