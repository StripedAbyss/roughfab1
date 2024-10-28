#pragma once
#include "liblgp.hpp"
#include "libhgp.h"
#include "RI.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Do_intersect_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/Ray_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::FT FT;  // 浮点数类型
typedef K::Segment_2 Segment;  // 二维线段类型
typedef K::Vector_2 Vector_2;

typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;  // Alpha Shape 顶点类型
typedef CGAL::Alpha_shape_face_base_2<K> Fb;    // Alpha Shape 面类型
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point_2 Point;

typedef CGAL::Partition_traits_2<K> Traits;
typedef Traits::Polygon_2 Polygon2;
typedef std::list<Polygon2> Polygon_list;
typedef CGAL::Ray_2<K> Ray;
typedef K::Triangle_2 Triangle;


using namespace libhgp;
using namespace liblgp;

namespace utils {
    Polygon_2 Convert_Vector2d1_to_Polygon_2(Vector2d1 v2d);//类型转换Vector2d1转化为Polygon
    Vector2d1 Convert_Polygon_2_to_Vector2d1(Polygon_2 p2);//类型转换Polygon转化为Vector2d1
    double pointToPolygonDist(const Point_2& p, const Polygon_2& polygon);
}