#ifndef ANA_MODES_H
#define ANA_MODES_H
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using std::size_t;

namespace ANA {
namespace NDD {

    class Modes {
    public:
        Modes(std::string const &filename);

        void get_modes_from_raw(std::string_view const texto);

        std::vector<std::vector<double>> atm_evectors;

        std::vector<std::vector<double>> evectors;

        std::vector<double> evals;

        size_t i, j;
    };

}
}
#endif // _H
