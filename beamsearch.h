#pragma once
#include "tree.h"
#include "polygons.h"
#include "output.h"
//#include "tinyxml2.hpp"
#include <opencv2/opencv.hpp>

using namespace std;
using namespace utils;

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

struct SegmentComparator {
    bool operator()(const Segment_2& a, const Segment_2& b) const {
        // 实现比较逻辑，返回 true 如果 seg1 应该排在 seg2 之前
        if (a.source().x() != b.source().x()) {
            return a.source().x() < b.source().x();
        }
        else if (a.source().y() != b.source().y()) {
            return a.source().y() < b.source().y();
        }
        else if (a.target().x() != b.target().x()) {
            return a.target().x() < b.target().x();
        }
        else {
            return a.target().y() < b.target().y();
        }

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

