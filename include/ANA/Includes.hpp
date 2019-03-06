#ifndef ANA_INCLUDES_H
#define ANA_INCLUDES_H

#include "chemfiles.hpp"
#include <algorithm>
#include <array>
#include <assert.h>
#include <boost/program_options.hpp>
#include <cmath>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Origin.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/convex_hull_3.h>

struct VertexInfo {
public:
    VertexInfo() = default;

    VertexInfo(int const index, double const radius, int const resn,
        std::string_view const resi) :
        _index(index),
        _radius(radius), _resn(resn), _resi(resi) {}

    int _index;
    double _radius;
    // atom's residue number
    int _resn;
    // atom's residue name in 3 letter format
    std::string _resi;
};

// clang-format off
// Basic definitions
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, EPIC>;
using Tds = CGAL::Triangulation_data_structure_3<Vb>;
using Delaunay = CGAL::Delaunay_triangulation_3<EPIC, Tds>;

// Definitions for Delaunay triangulation
using CVector = CGAL::Vector_3<EPIC>;
using Segment = CGAL::Segment_3<EPIC>;
using CPoint = Delaunay::Point;
using Vertex_handle = Delaunay::Vertex_handle;
using Cell_handle = Delaunay::Cell_handle;
using Vertex_iterator = Delaunay::Vertex_iterator;
using Finite_vertices_iterator = Delaunay::Finite_vertices_iterator;
using All_vertices_iterator = Delaunay::All_vertices_iterator;
using Edge_iterator = Delaunay::Edge_iterator;
using Facet_iterator = Delaunay::Facet_iterator;
using Finite_facets_iterator = Delaunay::Finite_facets_iterator;
using Cell_iterator = Delaunay::Cell_iterator;
using Finite_cells_iterator = Delaunay::Finite_cells_iterator;
using All_cells_iterator = Delaunay::All_cells_iterator;
using Cell_circulator = Delaunay::Cell_circulator;

// Definitions for convex hull
using Polyhedron = CGAL::Polyhedron_3<EPIC>;
using P_Facet_iterator = Polyhedron::Facet_iterator;
using P_Facet_const_iterator = Polyhedron::Facet_const_iterator;
using P_Edge_iterator = Polyhedron::Edge_iterator;
using P_Edge_const_iterator = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator =
    Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator =
    Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator = Polyhedron::Vertex_const_iterator;

// Miscellaneous definitions
using Object = CGAL::Object;
using CTriangle = CGAL::Triangle_3<EPIC>;
using Triang_Vector = std::vector<CTriangle>;
using Segment = EPIC::Segment_3;
using CTetrahedron = CGAL::Tetrahedron_3<EPIC>;
using Tetra_Vector = std::vector<CTetrahedron>;
using Box =
    CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Finite_cells_iterator>;

// ANA definitions
using MD_Element = std::array<CPoint, 4>;
// Cell
using MD_Vector = std::vector<MD_Element>;
// Pocket
using MD_Matrix = std::vector<MD_Vector>;
// All voids
using NDD_Element = std::array<std::pair<CPoint, double>, 4>;
// Cell
using NDD_Vector = std::vector<NDD_Element>;
// Pocket
using NDD_Matrix = std::vector<NDD_Vector>;
// All voids
using NDD_IElement = std::array<int, 4>;
// Cell indices
using NDD_IVector = std::vector<NDD_IElement>;
// Pocket indices
using NDD_IMatrix = std::vector<NDD_IVector>;
// All voids indices
using NA_Vector = std::vector<Finite_cells_iterator>;
// Pocket
using NA_Matrix = std::vector<NA_Vector>;
// All voids
using Poly_Vector = std::vector<Polyhedron>;
// Pocket border cells
using Poly_Matrix = std::vector<Poly_Vector>;
// All null areas border cells
using ANA_molecule = std::vector<std::pair<CPoint, VertexInfo>>;

#endif
