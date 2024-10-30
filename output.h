#pragma once
#include "polygons.h"
#include <opencv2/opencv.hpp>

namespace utils {
    struct PolygonWithCutDir {
        Vector2d1 polygon;  // 多边形各点坐标
        std::vector<pair<int, int>> cutIndex;   //切割射线起点与方向点的index
        // 构造函数
        PolygonWithCutDir() {}
        PolygonWithCutDir(const Vector2d1 &polygon) : polygon(polygon){}
        PolygonWithCutDir(const Vector2d1 &polygon, const std::vector<pair<int, int>> & cutIndex) : polygon(polygon), cutIndex(cutIndex){}
        inline void addAntiClockWiseCut() {
            int n = polygon.size();
            if (n < 3)return;
            for (int i = 0; i < n - 1; ++i) {
                addCut(i, i + 1);
            }
            addCut(n - 1, 0);
        }
        inline void addCut(int x, int y) {
            cutIndex.push_back({ x,y });
        }
    };

    class outputLayer {
    private:
        double boxx, boxy;
        std::vector<PolygonWithCutDir> input;
    public:
        outputLayer(double boxx, double boxy) :boxx(boxx), boxy(boxy) {}
        inline void setX(double x) { boxx = x; }
        inline void setY(double y) { boxy = y; }
        inline void setXY(double x, double y) { setX(x); setY(y); }
        inline void setInput(const std::vector<PolygonWithCutDir>& vec) { input = vec; }
        inline void addPolygon(const PolygonWithCutDir& polygon) {
            input.push_back(polygon);
        }
        inline void addPolygon(const Vector2d1& polygon) {
            PolygonWithCutDir temp(polygon);
            input.push_back(temp);
        }
        inline void addPolygon(const Vector2d1& polygon, const std::vector<pair<int, int>>& cut) {
            PolygonWithCutDir temp(polygon, cut);
            input.push_back(temp);
        }
        void geometry_layer_output();
    };
}