#pragma once

#include <CGAL/Real_timer.h>

// Correct and slow
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Lazy_exact_nt.h>
//typedef CGAL::Lazy_exact_nt<double> NT;
//typedef CGAL::Simple_cartesian<NT> K;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>


typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_with_history_2<Traits_2> Arrangement_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2_NoHistory;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef CGAL::Arr_segment_2<K> Curve_2;

typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Face_handle Face_handle;
typedef Arrangement_2::Vertex_const_handle Vertex_const_handle;
typedef Arrangement_2::Vertex_iterator Vertex_iterator;

// Const iterators and handles
typedef Arrangement_2::Face_const_iterator Face_const_iterator;

//typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle Face_const_handle;

// Used to find which face a point is in 
//#include <CGAL/Arr_trapezoid_ric_point_location.h>
//typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> Point_location;

//#include <CGAL/Arr_naive_point_location.h>
//typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

#include <CGAL/Arr_walk_along_line_point_location.h>
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;

#include <CGAL/Intersections_2/Segment_2_Segment_2.h>  // Correct header for CGAL 6




//
// Used for computing Barycentric coordinates for fibers
//
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Kernel/global_functions_2.h>  // For bounded_side_2
#include <CGAL/Polygon_2.h>

typedef CGAL::Simple_cartesian<double> CartesianKernel;
typedef CartesianKernel::Point_2 CartesianPoint;
typedef CartesianKernel::Segment_2 CartesianSegment;
typedef CartesianKernel::Line_2 CartesianLine;
typedef CGAL::Polygon_2<CartesianKernel> CartesianPolygon_2;

#include <CGAL/centroid.h>
