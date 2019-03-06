#ifndef ANA_UTILS_H
#define ANA_UTILS_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Cavity.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Molecule.hpp>
#include <ANA/NDDUtils.hpp>
#include <ANA/Options.hpp>
#include <ANA/Primitives.hpp>
#include <ANA/Read.hpp>
#include <ANA/Write.hpp>

namespace ANA {

// Substract the volume filled with the 4 atoms from the total volume of
// the corresponding cell.
double refine_cell_volume(
    double const entire_cell_vol, Finite_cells_iterator const cell_iterator);

// Cluster neighbouring cells.
void cluster_cells_cgal(NA_Vector const &input_cells, NA_Matrix &output_cells,
    int const min_cells_cluster);

// Given a cell, get all neighbouring cells that haven't been discovered yet
void get_neighbors(NA_Vector const &input_cells,
    Finite_cells_iterator const query_cell, std::vector<int> &except,
    NA_Vector &output_cells);

// Cluster neighbouring cells. Iso-oriented boxes method
void cluster_cells_boxes(NA_Vector const &input_cells, NA_Matrix &output_cells);

// Keep cells that correspond to the included amino acids
void keep_included_aa_cells(NA_Vector const &input_cells,
    const std::vector<int> &aa_list, int const nbr_of_vertices_to_include,
    NA_Vector &output_cells);

// Discard exposed cells
void discard_ASA_dot_pdt_cm(CPoint const &cm,
    std::vector<CPoint> const &Calpha_xyz, double const min_dot,
    double const max_length, std::string const only_side_ASA,
    NA_Vector const &input_cells, NA_Vector &output_cells);

// Discard exposed cells. Calpha convex hull method
void discard_ASA_CACH(std::vector<CPoint> const &Calpha_xyz,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells);

// Discard exposed cells. Alpha shape method
void discard_ASA_dot_pdt_axes(std::vector<CPoint> const &Calpha_xyz,
    double const min_dot, double const max_length,
    std::string const only_side_ASA, NA_Vector const &input_cells,
    NA_Vector &output_cells);

// Triangulate only the specified amino acids or the whole molecule.
inline Delaunay triangulate(ANA_molecule const &molecule_points) {
    Delaunay T;
    T.insert(molecule_points.begin(), molecule_points.end());
    return T;
}
// Calc volume of the input cells.
inline double get_void_volume(NA_Vector const &input_cells) {
    double volumen = 0;

    double current_cell_vol;
    for (Finite_cells_iterator const fc_ite : input_cells) {
        current_cell_vol = volume(fc_ite);

        current_cell_vol = refine_cell_volume(current_cell_vol, fc_ite);

        volumen = volumen + current_cell_vol;
    }
    return volumen;
}
// Disregard lone cells.
inline void disregard_lone_cells(
    NA_Vector const &input_cells, NA_Vector &output_cells) {

    int i, j, cell_cnt = input_cells.size();

    Finite_cells_iterator cell_ite1, cell_ite2;

    for (i = 0; i != cell_cnt; ++i) {
        cell_ite1 = input_cells[i];
        for (j = 1; j != cell_cnt; ++j) {
            cell_ite2 = input_cells[j];

            if (cell_ite1->has_neighbor(cell_ite2)) {
                output_cells.push_back(cell_ite1);
                break;
            }
        }
    }
}
// Calculate area of the cell's facets
inline void cell_facets_areas(Finite_cells_iterator const cell_iterator,
    double &facet_0_area, double &facet_1_area, double &facet_2_area,
    double &facet_3_area) {

    facet_0_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(1)->point(), cell_iterator->vertex(2)->point(),
        cell_iterator->vertex(3)->point()));

    facet_1_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(2)->point(), cell_iterator->vertex(3)->point(),
        cell_iterator->vertex(0)->point()));

    facet_2_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(3)->point()));

    facet_3_area = CGAL::to_double(CGAL::squared_area(
        cell_iterator->vertex(0)->point(), cell_iterator->vertex(1)->point(),
        cell_iterator->vertex(2)->point()));

    return;
}

// Fill the 2 input vectors with iterators for the outer and inner cells
// respectively.
void partition_triangulation(
    Delaunay const &T, NA_Vector &outer_cells, NA_Vector &inner_cells);

// Calc volume and get the proper cells.
double get_all_voids(Delaunay const &T, NA_Vector &big_cells,
    CellFilteringOptions const cell_opts);

// Discard cells without a vertex inside the specified convex hull. Lo
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells);

// Discard cells without a vertex inside the specified convex hull. Hi
// precision.
void discard_CH_0(NA_Vector const &in_cells, Triang_Vector const &CH_triangs,
    NA_Vector &out_cells, NA_Vector &out_intersecting_cells,
    std::vector<std::array<bool, 4>> &intersecting_bool,
    std::vector<int> &intersecting_total);

// Discard parts of cells outside the specified triangulation using
// intersecitons.
double discard_CH_1(NA_Vector const &in_intersecting_cells,
    Triang_Vector const &CH_triangs,
    const std::vector<std::array<bool, 4>> &intersecting_bool,
    const std::vector<int> &intersecting_total, Poly_Vector &border_poly,
    std::vector<std::array<double, 3>> &in_vtces_radii, int &atom_cnt_poly);

// Discard cells without a vertex inside the specified convex hull
void discard_CH(
    NA_Vector const &in_cells, Polyhedron &CH, NA_Vector &out_cells);

// Extract vertices coordinates from the cells and store them in an "md_vector".
void na_vector_into_ndd_vector(
    NA_Vector const &in_cells, NDD_Vector &out_cells);

// Tool for reading PDB to draw included area.
void tool_PDB_to_CH(
    std::string const &in_filename, std::string const &out_filename);

// Tool for normalizing PDB, by renumbering its atoms and residues.
void tool_PDB_norm(
    std::string const &in_filename, std::string const &tool_pdb_norm);

// Helper function for inserting elements in ordered vectors.
template <class CVector, class T>
void insert_into_ord_vtor(CVector &v, const T &to_insert) {
    typename CVector::iterator i =
        std::lower_bound(v.begin(), v.end(), to_insert);
    if (i == v.end() || to_insert < *i) {
        v.insert(i, to_insert);
    }
    return;
}

// Helper function for taking the "i" number in the 'in_vec'' that doesn't
// match the query.
template <class T>
T get_i_not_equal(const std::vector<T> &in_vec, const T &query, int const i) {
    if (static_cast<int>(in_vec.size()) < i) {
        throw std::invalid_argument("get_i_not_equal(): specified \"i\" "
                                    "position larger than input vector");
    }

    int cont = 0;
    for (auto const &each : in_vec) {
        if (each != query) {
            ++cont;
            if (cont >= i) {
                return each;
            }
        }
    }

    // Fail
    return i;
}

// Helper function for taking the "i" number in the 'in_vec'' that doesn't
// match the query vector.
template <class T>
T get_i_not_equal(const std::vector<T> &in_vec, const std::vector<T> &query_vec,
    int const i) {
    if (static_cast<int>(in_vec.size()) < i) {
        throw std::invalid_argument("get_i_not_equal(): specified \"i\" "
                                    "position larger than input vector");
    }

    int cont = 0;
    for (auto const &each : in_vec) {
        bool each_bool = true;
        for (auto const &query : query_vec) {
            each_bool = each_bool && (each != query);
        }
        if (each_bool) {
            ++cont;
            if (cont >= i) {
                return each;
            }
        }
    }

    // Fail
    return i;
}

// Helper function for getting the indices that sort a vector.
template <typename T>
std::vector<int> sort_indices(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indices based on comparing values in v
    sort(
        idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] < v[i2]; });

    return idx;
}

// Helper function to use binary search to find the lowest bound of a query in
// a sorted vector in ascending order. It returns true if a match is found,
// and
// stores the index of the element that satisfies the lower bound condition in
// the variable "first".
template <typename T, typename K>
bool lb(const std::vector<T> &v1, const T q1, K &first) {

    static_assert(std::is_integral<K>::value);
    K count = v1.size(), step, current;
    first = 0;

    while (count > 0) {
        step = count / 2;
        current = first;
        current += step;

        if (v1[current] < q1) {
            first = ++current;
            count -= (step + 1);
        } else
            count = step;
    }

    // Did the query match?
    if (first == count) {
        return false;
    } else
        return true;
}

// Helper function to use binary search to find the lowest bound of a query in
// an unsorted vector and vector of indices that sorts it in ascending order.
// It returns true if a match is found, and stores the index of the
// element that satisfies the lower bound condition in the variable "first".
// This variable also serves as a starting point in the search, to start
// searching in an arbitrary position and forward.
template <typename T, typename K>
bool lb_with_indices(const std::vector<T> &v1, const std::vector<int> &indices,
    const T q1, K &first) {

    static_assert(std::is_integral<K>::value);
    K count = v1.size(), step, current;
    first = 0;

    while (count > 0) {
        step = count / 2;
        current = first;
        current += step;

        if (v1[indices[current]] < q1) {
            first = ++current;
            count -= (step + 1);
        } else
            count = step;
    }

    // Is the query larger?
    if (first == count) {
        return false;
    } else
        return true;
}

} // namespace ANA
#endif // _H
