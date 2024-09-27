#pragma once
#include "tree.h"
#include "liblgp.hpp"
#include "libhgp.h"
#include "RI.hpp"
//#include "tinyxml2.hpp"
#include <opencv2/opencv.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <opencv2/opencv.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Do_intersect_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
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
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;  // 三角剖分数据结构类型
typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation_2;  // Delaunay 三角剖分类型
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;  // Alpha Shape 类型
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;  // Alpha Shape 边迭代器
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point_2 Point;

typedef CGAL::Partition_traits_2<K> Traits;
typedef Traits::Polygon_2 Polygon2;
typedef std::list<Polygon2> Polygon_list;
typedef CGAL::Ray_2<K> Ray;
typedef K::Triangle_2 Triangle;

using namespace std;
using namespace libhgp;
using namespace liblgp;

struct Candidate {
    int CandidateId = 0;
    std::vector<Vector2d1> polygons;  // 解决方案中的多边形集合
    double score;                    // 解决方案的得分
    std::vector<int> typenum;          //解决方案中多边形序号
    int previous = 1;
    // 构造函数
    Candidate(int id, const std::vector<Vector2d1>& polys, double sc, const std::vector<int>& tn, int p) : CandidateId(id), polygons(polys), score(sc), typenum(tn), previous(p) {}
    bool operator<(const Candidate& other) const {
        return score < other.score; //排序
    }
    bool operator>(const Candidate& other) const {
        return score > other.score; //排序
    }
};

class Beamsearch {
public:
    string image_path = "image_outputs";
    Tree gml_tree;
    vector<Vector2d1> origin_polygons;//未粗加工化多边形集合
    vector<Vector2d1> polygons;//多边形集合
    vector<vector<Vector2d1>> process_solutions;
    vector<double> score;
    //void generateNodeGML(TreeNode* node, std::ofstream& outfile);
    //void generateEdgeGML(TreeNode* node, std::ofstream& outfile);
    //void Convert_GmlTree_To_GML();
    pair<double, int> calculateScore(const std::vector<Vector2d1>& polygons, int previous);
    bool doPolygonsCollide2(const Vector2d1& poly1, const vector<Vector2d1>& poly2);
    Vector2d1 translatePolygon(const Vector2d1& polygon, double dx, double dy);
    std::vector<Vector2d1> beamSearch(const std::vector<Vector2d1>& inputPolygons, int beamWidth, const Vector2d1& boundingRect);
    void work();
    void test();
    void get_points_to_polygon();
    std::vector<Vector2d1> perior_geometry_put();
    void geometry_layer_output(vector<Vector2d1> a);
    void geometry_layer_output2(vector<Vector2d1> a, vector<Vector2d1> b);

    void geometry_layer_save(vector<Vector2d1> a, int num, double score);
    void geometry_layer_save1(vector<Vector2d1> a, vector<Vector2d1> b);
    vector<double> GetScore();
    void PolygonModification();
    void PolygonModification1();
    void PolygonModification2();
};

