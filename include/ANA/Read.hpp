#ifndef ANA_READ_H
#define ANA_READ_H
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/PrimitivesUtils.hpp>

namespace ANA {
// Refine the provided list of amino acids. If its not present, then return an
// array of two '0' elements.
template <class T>
bool adapt_AA_list(
    std::string &aa_list_proto, std::vector<T> &aa_list, int top = 0) {
    int aa;

    if (aa_list_proto == "none") {
        aa_list.push_back(0);
        aa_list.push_back(0);
        return false;

    } else {
        // borro el caracter ' pq me caga todo
        std::replace(aa_list_proto.begin(), aa_list_proto.end(), '\'', ' ');
        std::stringstream stream_aa(aa_list_proto);
        std::string temp_aa;

        while (!stream_aa.eof()) {
            stream_aa >> temp_aa;
            try {
                aa = std::stoi(temp_aa);
            } catch (std::invalid_argument const &ia) {
                // some character present. Doesn't matter, skip it and keep
                // reading.
                continue;
            } catch (std::out_of_range const &oor) {
                // int is too large to be represented by int
                std::cerr << "Invalid residue number: " << temp_aa << '\n';
                std::exit(1);
            } catch (...) {
                // some other exception. Doesn't matter, skip it and keep
                // reading.
                std::cerr << "Invalid residue number: " << temp_aa
                          << ". Will keep reading." << '\n';
                continue;
            }
            aa_list.push_back(aa);
        }
        // sort list of included amino acids
        std::sort(aa_list.begin(), aa_list.end());

        if (top != 0) {
            if (aa_list[aa_list.size() - 1] > top) {
                std::cerr
                    << "Atom / residue list goes out of bounds. Check this "
                       "input list and your input PDB atom count. Quiting now."
                    << '\n';
                exit(0);
            }
        }

        return true;
    }
}

// Read coordinates in pdb format using chemfiles.
bool read_static(std::string const &filename,
    bool const triangulate_only_included_aas, bool const atom_only,
    std::string &aa_list_proto, std::string &exclude_ca_for_ASA_proto,
    std::string &include_CH_aa_proto, std::string &include_CH_atom_proto,
    std::string &sphere_proto, std::string &cylinder_proto,
    std::string &prism_proto, std::string const &include_CH_filename,
    ANA_molecule &molecule_points, CPoint &cm, std::vector<int> &aa_list,
    std::vector<int> &CA_indices, std::vector<CPoint> &CAs_Points,
    std::vector<int> &include_CH_aa, Triang_Vector &CH_triangs,
    std::vector<int> &hetatm_atoms);

// Read coordinates in netcdf format.
void read_MD(const chemfiles::Frame &in_frame, bool const requested_CH,
    std::string const &sphere_proto, std::string const &cylinder_proto,
    std::string const &prism_proto, const std::vector<int> &hetatm_atoms,
    std::vector<int> &include_CH_atoms, std::string const &include_CH_filename,
    Triang_Vector &CH_triangs, std::string const &ASA_method,
    const std::vector<int> &CA_indices, std::vector<CPoint> &CAs_points,
    ANA_molecule &molecule_points);

// Get the center of mass from a chemfiles molecule.
CPoint getCM(
    const chemfiles::span<chemfiles::Vector3D> &in_xyz, int const natoms);

// Read PDB to draw included area for MD
void read_included_area(
    std::string const &filename, std::vector<CPoint> &area_points);
}
#endif // _H
