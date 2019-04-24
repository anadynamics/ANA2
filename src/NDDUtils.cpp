#include <ANA/NDDUtils.hpp>

namespace ANA::NDD {

// On-site NDD.
void ndd(
    Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts) {

    std::vector<double> output_volumes;

    Modes const modos = create_modes(NDD_opts);
    std::vector<double> scaling_factors =
        initialize_scaling_factors(NDD_opts, modos._j);

    std::vector<double> pos_vols_ndd, neg_vols_ndd, der_vols_ndd;
    pos_vols_ndd.reserve(modos._j);
    neg_vols_ndd.reserve(modos._j);
    der_vols_ndd.reserve(modos._j);

    if (CH._dynamic) {
        for (size_t j = 0; j < modos._j; ++j) {
            double const mul = NDD_opts._step / scaling_factors[j];

            // In the positive direction.
            ConvexHull CH_pos(CH, modos._evectors[j], mul);
            Cavity hueco_pos(hueco, modos._evectors[j], mul);
            carve_CH_into_cavity(hueco_pos, CH_pos);
            double const pos_vol = hueco_pos._volume + hueco_pos._outer_volume;

            // In the negative direction.
            ConvexHull CH_neg(CH, modos._evectors[j], -mul);
            Cavity hueco_neg(hueco, modos._evectors[j], -mul);
            carve_CH_into_cavity(hueco_neg, CH_neg);
            double const neg_vol = hueco_neg._volume + hueco_neg._outer_volume;

            // Numerical derivative.
            double const der_vol {(pos_vol - neg_vol) / mul};

            der_vols_ndd.push_back(der_vol);
            pos_vols_ndd.push_back(pos_vol);
            neg_vols_ndd.push_back(neg_vol);
        }

        if (NDD_opts._derivative) {
            std::string const filename = std::to_string(NDD_opts._step) + "_" +
                NDD_opts._out_ndd_filename;
            write_vector(der_vols_ndd, filename);
        } else {
            std::string const filename_pos = std::to_string(NDD_opts._step) +
                "_pos_" + NDD_opts._out_ndd_filename;
            write_vector(pos_vols_ndd, filename_pos);

            std::string const filename_neg = std::to_string(NDD_opts._step) +
                "_neg_" + NDD_opts._out_ndd_filename;
            write_vector(neg_vols_ndd, filename_neg);
        }

    } else {
        ;
    }

    return;
}

// If scaling values were provided, read them. If not, set a uniform distro
// of 1.
auto initialize_scaling_factors(NDDOptions const &NDD_opts, std::size_t j)
    -> std::vector<double> {
    if (NDD_opts._scaling_ndd_filename != "none") {

        std::unique_ptr<char[]> const buffer_evals =
            slurp(NDD_opts._scaling_ndd_filename);
        std::vector<double> scaling_factors =
            get_values_from_raw(std::string_view(buffer_evals.get()));

        if (scaling_factors.size() != j) {
            std::cerr << "Vector count: " << j
                      << ". Scalar count: " << scaling_factors.size() << '\n';
            throw std::runtime_error(
                "Frequencies don't match vectors. Aborting.");
        }

        return scaling_factors;
    }

    std::vector<double> scaling_factors;
    scaling_factors.reserve(j);
    for (std::size_t i = 0; i < j; ++i) {
        scaling_factors.push_back(1.);
    }

    return scaling_factors;
}

} // namespace ANA::NDD