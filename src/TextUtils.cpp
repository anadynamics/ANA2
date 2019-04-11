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

//
auto guess_format(std::string_view texto)
    -> std::tuple<size_t, size_t, size_t> {

    size_t ch_elem = 0, line_length = 0;
    size_t const fsz = texto.size();
    bool done_ch_elem = false, pre_blank = false;

    for (size_t c = 0; c < fsz; ++c) {
        if (iscntrl(texto[c])) {
            line_length = c;
            if (!done_ch_elem) {
                // newline char read before a space char. There is only 1
                // element per line.
                ch_elem = line_length;
            }
            break;
        }
        if (!done_ch_elem) {
            if (isspace(texto[c])) {
                // In case there are leading whitespaces.
                if (c == 0 || pre_blank) {
                    pre_blank = true;
                    continue;
                }
                ch_elem = c;
                done_ch_elem = true;
            } else
                pre_blank = false;
        }
    }

    if (line_length == 0 || line_length % ch_elem != 0 ||
        fsz % (line_length + 1) != 0) {
        throw std::runtime_error(
            "guess_format(): corrupted input. Could not read them.\nRemember "
            "to use fixed-width format.");
    }

    return {ch_elem, line_length / ch_elem, fsz / (line_length + 1)};
}

} // namespace ANA