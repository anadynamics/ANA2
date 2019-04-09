#include <ANA/NDDUtils.hpp>

namespace ANA::NDD {

// On-site NDD.
void ndd(
    Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts) {

    std::vector<double> output_volumes;

    Modes const modos(NDD_opts._modes_ndd_filename, NDD_opts._amber_modes);
    std::vector<double> pos_vols_ndd, neg_vols_ndd, der_vols_ndd;
    pos_vols_ndd.reserve(modos._j);
    neg_vols_ndd.reserve(modos._j);
    der_vols_ndd.reserve(modos._j);

    if (CH._dynamic) {
        for (size_t j = 0; j < modos._j; ++j) {
            double const mul = NDD_opts._step / modos._evals[j];
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
            double const der_vol =
                ((pos_vol - neg_vol) * modos._evals[j]) / NDD_opts._step;
            der_vols_ndd.push_back(der_vol);
            pos_vols_ndd.push_back(pos_vol);
            neg_vols_ndd.push_back(neg_vol);
        }

        std::string const filename =
            std::to_string(NDD_opts._step) + "_" + NDD_opts._out_ndd_filename;
        write_vector(der_vols_ndd, filename);
    } else {
        ;
    }

    return;
}

} // namespace ANA::NDD