#include "output.h"

namespace utils {
	void outputLayer::geometry_layer_output() {
        // ����ͼ��ĳߴ�

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_64FC3����ʼֵΪ��ɫ
        cv::Mat image(boxy, boxx, CV_64FC3, cv::Scalar(0, 0, 0));

        // ���ƶ����
        for (const auto& a : input) {
            std::vector<cv::Point> points;
            for (const auto& vertex : a.polygon) {
                // ����������ת��ΪOpenCVͼ������ϵ�е�����
                int x = vertex.x;
                int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
                cv::Point point(x, y);
                points.push_back(point);
            }
            const cv::Point* pts = points.data();
            int num_points = points.size();
            // ���ƶ��������
            cv::polylines(image, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);

            for (const auto& pp : a.cutIndex) {
                if (pp.first > a.polygon.size() || pp.second > a.polygon.size()) {
                    std::cout << "cutIndex OOR" << std::endl;
                    continue;
                }
                int x1 = a.polygon[pp.first].x;
                int y1 = boxy - a.polygon[pp.first].y;
                int x2 = a.polygon[pp.second].x;
                int y2 = boxy - a.polygon[pp.second].y;

                auto exPair = rayExtension(x1, y1, x2, y2);
                x2 = exPair.first;
                y2 = exPair.second;

                cv::Point point1(x1, y1);
                cv::Point point2(x2, y2);
                cv::line(image, point1, point2, cv::Scalar(0, 255, 255), 1);
            }

        }

        cv::imshow("Polygons", image);
        cv::waitKey(0);
        return;
	}

}
