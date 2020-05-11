#include <ANA/NDDUtils.hpp>

namespace ANA::NDD {

// On-site NDD.
void ndd(Cavity const &hueco, ConvexHull const &CH, NDDOptions const &NDD_opts,
    std::string const &pdb_filename) {

    std::vector<double> output_volumes;

    Modes const modos = create_modes(NDD_opts, pdb_filename);

    // write_matrix(modos._atm_evectors, modos._ii, modos._j, "atm_evectors");

    std::vector<double> scaling_factors =
        initialize_scaling_factors(modos, NDD_opts);

    std::vector<double> neg_vols_ndd, pos_vols_ndd, vgv;
    neg_vols_ndd.reserve(modos._j);
    pos_vols_ndd.reserve(modos._j);
    vgv.reserve(modos._j);

    if (CH._dynamic) {
        // update the convex hull (included area) if needed.
        for (size_t j = 0; j < modos._j; ++j) {
            double const mul = NDD_opts._size / scaling_factors[j];

            // In the negative direction.
            ConvexHull CH_neg(CH, modos._atm_evectors[j], -mul);
            Cavity hueco_neg(hueco, modos._atm_evectors[j], -mul);
            carve_CH_into_cavity(hueco_neg, CH_neg);
            double const neg_vol = hueco_neg._volume + hueco_neg._outer_volume;

            // In the positive direction.
            ConvexHull CH_pos(CH, modos._atm_evectors[j], mul);
            Cavity hueco_pos(hueco, modos._atm_evectors[j], mul);
            carve_CH_into_cavity(hueco_pos, CH_pos);
            double const pos_vol = hueco_pos._volume + hueco_pos._outer_volume;

            // 1st step
            pos_vols_ndd.push_back(pos_vol);
            neg_vols_ndd.push_back(neg_vol);

            // 2nd step: numerical derivative.
            double const der_vol{(pos_vol - neg_vol) / mul};
            vgv.push_back(der_vol);
        }

    } else {
        ;
    }

    write_result(NDD_opts, modos, neg_vols_ndd, pos_vols_ndd, vgv);

    return;
}

// If scaling values were provided, read them. If not, set a linear
// distribution.
auto initialize_scaling_factors(Modes const &modos, NDDOptions const &NDD_opts)
    -> std::vector<double> {

    // _scaling_ndd_filename overrides _scale_w_freqs if both are set.
    if (NDD_opts._scaling_ndd_filename != "none") {

        auto[bufr_scaling_ftor, fsz] = slurp(NDD_opts._scaling_ndd_filename);
        std::vector<double> scaling_factors =
            get_values_from_raw(std::string_view(bufr_scaling_ftor.get(), fsz));

        if (scaling_factors.size() != modos._j) {
            std::cerr << "Vector count: " << modos._j
                      << ". Scalar count: " << scaling_factors.size() << '\n';
            throw std::runtime_error(
                "Frequencies don't match vectors. Aborting.");
        }

        return scaling_factors;
    } else if (NDD_opts._scale_w_freqs) {
        return modos._evalues;
    } else {

        std::vector<double> scaling_factors;
        scaling_factors.reserve(modos._j);
        for (std::size_t i = 0; i < modos._j; ++i) {
            scaling_factors.push_back(static_cast<double>(i + 5 * 0.5));
        }
        return scaling_factors;
    }
}

void write_result(NDDOptions const &NDD_opts, Modes const &modos,
    std::vector<double> const &neg_vols_ndd,
    std::vector<double> const &pos_vols_ndd, std::vector<double> const &vgv) {

    switch (NDD_opts._step) {
    case NDDOptions::Volumes: {
        std::string const filename_neg = std::to_string(NDD_opts._size) +
            "_neg_" + NDD_opts._out_ndd_filename;
        std::string const filename_pos = std::to_string(NDD_opts._size) +
            "_pos_" + NDD_opts._out_ndd_filename;

        write_vector(neg_vols_ndd, filename_neg);
        write_vector(pos_vols_ndd, filename_pos);
        break;
    }
    case NDDOptions::Gradient: {
        std::string const filename =
            std::to_string(NDD_opts._size) + "_" + NDD_opts._out_ndd_filename;
        write_vector(vgv, filename);
        break;
    }
    case NDDOptions::Index: {
        barletta_index(modos, vgv);
        break;
    }
    }
    return;
}

void barletta_index(Modes const &modos, std::vector<double> const &vgv) {
    // VGV needs to be normalized before calculating the flexibility index.
    double squared_norm = 0;
    for (std::size_t j = 0; j < modos._j; ++j) {
        squared_norm += vgv[j] * vgv[j];
    }
    // Instead of dividing each VGV element by the vector's norm, just
    // divide each squared element by the squared norm.
    double sum = 0;
    for (std::size_t j = 0; j < modos._j; ++j) {
        sum += modos._evalues[j] * modos._evalues[j] * vgv[j] * vgv[j] /
            squared_norm;
    }

    double const barletta_index = 1 / (cte * sum);
    printf("Flexibility:  %.10f\n", barletta_index);

    return;
}

} // namespace ANA::NDD