#include "polygons.h"

namespace utils {


Polygon_2 Convert_Vector2d1_to_Polygon_2(Vector2d1 v2d)//����ת��Vector2d1ת��ΪPolygon
{
    vector<Point_2> pts;
    for (auto itt = v2d.begin(); itt != v2d.end(); itt++)
    {
        Point_2 temp((*itt).x, (*itt).y);
        pts.push_back(temp);
    }
    Polygon_2 target(pts.begin(), pts.end());
    return target;
}

Vector2d1 Convert_Polygon_2_to_Vector2d1(Polygon_2 p2)//����ת��Polygonת��ΪVector2d1
{
    Vector2d1 target;
    for (auto itt = p2.begin(); itt != p2.end(); itt++)
    {
        Vector2d temp((*itt).x(), (*itt).y());
        target.push_back(temp);
    }
    return target;
}

double pointToPolygonDist(const Point_2& p, const Polygon_2& polygon) {//�㵽����εľ���;������ɣ�û����
    int count = 0;
    double minDist = INFINITY;

    for (auto e = polygon.edges_begin(); e != polygon.edges_end(); ++e) {
        const Point_2& a = e->source();
        const Point_2& b = e->target();

        // �жϵ� p ���߶� ab �Ƿ���ͬһˮƽ���ϣ����� p �� ab �����
        if ((a.y() > p.y() != b.y() > p.y()) &&
            (p.x() < (b.x() - a.x()) * (p.y() - a.y()) / (b.y() - a.y()) + a.x())) {
            count++;
        }

        // ����� p ���߶� ab ����̾��룬��������С����
        Segment_2 s(a, b);
        minDist = std::min(minDist, squared_distance(p, s));
    }

    if (count % 2 == 0) {
        return std::sqrt(abs(minDist));

    }
    return -std::sqrt(abs(minDist));
}

bool doPolygonsCollide2(const Vector2d1& poly1, const vector<Vector2d1>& poly2) {//��ײ��⣬�������
    for (const Vector2d1& one_polygon : poly2) {
        if (PL().HGP_2D_Two_Polygons_Intersection_C(poly1, one_polygon) > 0) {
            return true; // ������ײ
        }
    }
    return false; // δ������ײ
}

Vector2d1 translatePolygon(const Vector2d1& polygon, double dx, double dy) {
    std::vector<Vector2d> translatedVertices;
    // �������ж��㣬��ÿ���������ƽ�Ʋ���������ӵ��µĶ����б���
    for (auto it = polygon.begin(); it != polygon.end(); ++it) {

        Vector2d translatedPoint((*it).x + dx, (*it).y + dy);
        translatedVertices.push_back(translatedPoint);
    }
    // ʹ���µĶ����б���һ���µ�Vector2d1���󲢷���
    return Vector2d1(translatedVertices.begin(), translatedVertices.end());
}

pair<int, int> rayExtension(int x1, int y1, int x2, int y2) {
    int deltaY = y2 - y1;
    int deltaX = x2 - x1;
    y2 = y1 + deltaY * 100;
    x2 = x1 + deltaX * 100;
    return pair<int, int>(x2, y2);
}

void GetMarkEdges(const Polygon_2 &poly, vector<pair<Segment_2, int>> &MarkEdges) {

    for (Polygon_2::Edge_const_iterator itt = poly.edges_begin(); itt != poly.edges_end(); ++itt) {
        Point_2 start = itt->source();
        Point_2 end = itt->target();
        Ray ray1(start, end);
        Ray ray2(end, start);
        int nums1 = 0, nums2 = 0;
        for (Polygon_2::Edge_const_iterator w = poly.edges_begin(); w != poly.edges_end(); ++w) {
            if ((w->source() == start && w->target() == end) || (w->source() == end && w->target() == start))continue;
            if (CGAL::do_intersect(ray1, *w))nums1++;
            if (CGAL::do_intersect(ray2, *w))nums2++;

        }
        pair<Segment_2, int> MarkEdge;
        MarkEdge.first = (*itt);
        MarkEdge.second = 0;
        assert(nums1 >= 2);
        assert(nums2 >= 2);
        if (nums1 > 2) {
            MarkEdge.second |= 1;
        }
        if (nums2 > 2) {
            MarkEdge.second |= 2;
        }
        MarkEdges.push_back(MarkEdge);
    }
}

void GetCutDir(const Polygon_2 &poly, vector<pair<int, int>> &cutIndex) {
    int n = poly.size();
    if (n < 3)return;
    for (int i = 0; i < n; ++i) {
        int j = i + 1;
        if (j == n) {
            j = 0;
        }
        Point_2 start = poly[i];
        Point_2 end = poly[j];
        Ray ray1(start, end);
        Ray ray2(end, start);
        int nums1 = 0, nums2 = 0;
        for (Polygon_2::Edge_const_iterator w = poly.edges_begin(); w != poly.edges_end(); ++w) {
            if ((w->source() == start && w->target() == end) || (w->source() == end && w->target() == start))continue;
            if (CGAL::do_intersect(ray1, *w))nums1++;
            if (CGAL::do_intersect(ray2, *w))nums2++;
        }
        if (nums1 == 2) {
            cutIndex.push_back({ i,j });
        }
        else if (nums2 == 2) {
            cutIndex.push_back({ j,i });
        }
     
    }
}

}