#include "beamsearch.h"
double boxx, boxy;
// �������ڶ��������ıȽϺ���
bool beamssort(const Vector2d1& poly1, const Vector2d1& poly2)
{
    // ��������洢����ε�����ֵ
    double xmax_score1, xmin_score1, xlength_score1, ylength_score1, area_score1;
    double xmax_score2, xmin_score2, xlength_score2, ylength_score2, area_score2;
    Vector2d max1, min1, max2, min2;
    // �����һ������εı߽������ֵ
    PL().HGP_2d_Polygon_Boundingbox_C(poly1, min1, max1);
    xmax_score1 = max1.x / boxx;
    xmin_score1 = min1.x / boxx;
    xlength_score1 = (max1.x - min1.x) / boxx;
    ylength_score1 = (max1.y - min1.y) / boxy;
    // ����ڶ�������εı߽������ֵ
    PL().HGP_2d_Polygon_Boundingbox_C(poly2, min2, max2);
    xmax_score2 = max2.x / boxx;
    xmin_score2 = min2.x / boxx;
    xlength_score2 = (max2.x - min2.x) / boxx;
    ylength_score2 = (max2.y - min2.y) / boxy;
    // �������ε��������ֵ
    area_score1 = PL().HGP_2D_Polygon_Area_C(poly1) / (boxx * (max1.y - min1.y));
    area_score2 = PL().HGP_2D_Polygon_Area_C(poly2) / (boxx * (max2.y - min2.y));
    // �������ε��ܵ÷�
    double score1 = xmax_score1 + xmin_score1 + xlength_score1 + ylength_score1 + area_score1;
    double score2 = xmax_score2 + xmin_score2 + xlength_score2 + ylength_score2 + area_score2;
    // �Ƚ���������ε��ܵ÷�
    return score1 > score2;
}


bool comPoints(const Point_2& p1, const Point_2& p2) {
    if (p1.x() != p2.x()) {
        return p1.x() < p2.x();
    }
    return p1.y() < p2.y();
}

pair<double, int> Beamsearch::calculateScore(const std::vector<Vector2d1>& polygons, int previous)
{
    double offset = 15;
    double areas = 0;
    double score = 0;
    boxx /= 2;
    int now = 0;
    vector<Polygon_2> pys;
    vector<Polygon_2> ans;
    vector<Vector2d1> output;
    for (auto it = polygons.begin(); it != polygons.end(); it++)
    {
        pys.push_back(Convert_Vector2d1_to_Polygon_2(*it));
    }
    vector<Point_2> getit;
    vector<Vector2d1> outerFrame{ { Vector2d(-offset, -offset - 1),Vector2d(boxx + offset, -offset - 1),Vector2d(boxx + offset, -offset),Vector2d(-offset, -offset) },
                     { Vector2d(boxx + offset + 1, -offset),Vector2d(boxx + offset + 1, boxy + offset),Vector2d(boxx + offset, boxy + offset),Vector2d(boxx + offset, -offset) },
                     { Vector2d(boxx + offset, boxy + offset + 1),Vector2d(-offset, boxy + offset + 1),Vector2d(-offset, boxy + offset),Vector2d(boxx + offset, boxy + offset) },
                     { Vector2d(-offset - 1, boxy + offset),Vector2d(-offset - 1, -offset),Vector2d(-offset, -offset),Vector2d(-offset, boxy + offset) }
    };
    boxx *= 2;
    vector<Vector2d1> polygonsWithFrame{polygons};
    polygonsWithFrame.insert(polygonsWithFrame.end(), outerFrame.begin(), outerFrame.end());
    PL().HGP_2D_Polygons_One_Offsets_C(polygonsWithFrame, -offset, output);

    for (auto it = output.begin(); it != output.end(); it++) {
        ans.push_back(Convert_Vector2d1_to_Polygon_2(*it));
    }

    geometry_layer_output2(output, polygons);
    for (auto it = ans.begin(); it != ans.end(); it++)
    {
        if (abs(it->area()) < 500) {//��С�ľ������迼��
            continue;
        }
        else {//������
            //���￼�ǵ��Ǿ�������������Χ��������ܳ����ܳ�
            areas += it->bbox().x_span() * it->bbox().y_span();
            double length = 0;
            now++;
            for (auto itt = it->edges_begin(); itt != it->edges_end(); itt++)
            {
                length += sqrt(itt->squared_length());
            }
            score += it->bbox().x_span() * it->bbox().y_span() * (it->bbox().x_span() + it->bbox().y_span()) * 2 / (length);
        }
    }
    if (areas == 0)return make_pair(0, now);
    score /= areas;
    score = score * 0.9 + 0.1 * min(previous / now, 1);//��ʽ�������ɵ���
    return make_pair(score, now);//���ص÷��뾧������
}

std::vector<Vector2d1> Beamsearch::beamSearch(const std::vector<Vector2d1>& inputPolygons, int beamWidth, const Vector2d1& boundingRect) {
    // ��ʼ����ѡ���������id
    int id = 0;
    // ����һ�����ȶ��У����ڴ洢��ѡ������������յ÷ִӸߵ�������
    std::priority_queue <Candidate, std::vector<Candidate>, less<Candidate>> candidates;
    // �������ڵ�
    Candidate root(id++, {}, 0.0, {}, 1);//�յĽڵ㣬û�м������Σ�����Ҳ��0
    // �����ڵ�����ѡ�����������
    candidates.push(root);
    // ����һ������ڵ㣬��־���ڵ�
    gml_tree.insert(-1, 0);//��϶������gml_tree���Ի󣬼���������ĵ�
    std::vector<Vector2d1> sortedPolygons = inputPolygons;
    // ����ԭʼ�����
    std::vector<Vector2d1> ori_Polygons = inputPolygons;
    std::sort(sortedPolygons.begin(), sortedPolygons.end(), beamssort);// ������Ķ���ν������򣬰���һ���Ĺ������ﰴ��beamsort����õ��Ƕ��������

    // ����ÿ�������õĶ����
    for (int times = 0; times < sortedPolygons.size(); times++) {
        // ���ڴ洢��һ�ִεĺ�ѡ������������ȶ���
        std::priority_queue < Candidate, std::vector<Candidate>, less<Candidate>> nextCandidates;

        // ����ǰ�ִε�ÿ����ѡ�������
        while (!candidates.empty()) {
            // ��ȡ��ǰ���ŵĺ�ѡ�������
            Candidate candidate = candidates.top();
            candidates.pop();
            // ���Ʒ��ô����ı���
            int tab = 0;

            // �����ڵ�ǰλ�÷��ò�ͬ�Ķ����
            for (int i = 0; i < sortedPolygons.size(); i++) {
                // ���ö�����Ƿ��Ѿ������ڽ��������
                bool type_tab = 0;
                for (auto types : candidate.typenum) {
                    if (types == i) {
                        type_tab = 1;
                        break;
                    }
                }
                // ���������Ѿ������ڽ�������У�������
                if (type_tab != 0) {
                    continue;
                }
                // ���Ʒ��ô�������ೢ��3��
                if (tab < 3) tab++;
                else break;

                // ������η��������������
                Vector2d bomin, bomax, somin, somax;
                PL().HGP_2d_Polygon_Boundingbox_C(boundingRect, bomin, bomax);
                PL().HGP_2d_Polygon_Boundingbox_C(sortedPolygons[i], somin, somax);
                double dx = 0.0;
                double dy = bomax.y - somax.y;
                Vector2d1 finalPolygon = translatePolygon(sortedPolygons[i], dx, dy);
                // ������ú�����ײ�򳬳��߽磬������
                if (doPolygonsCollide2(finalPolygon, candidate.polygons) || somax.y > boxy) {
                    tab--;
                    continue;
                }
                // ʹ�ö��ַ�����ƽ�ƣ�ֱ��������ײ
                PL().HGP_2d_Polygon_Boundingbox_C(finalPolygon, bomin, bomax);
                double bottom_distance = bomin.y;
                bool judge = 1;
                double pymin = bomin.y;
                while (bottom_distance > 10) {
                    pymin -= bottom_distance;
                    if (pymin < 0) {
                        pymin += bottom_distance;
                        bottom_distance /= 2.0;
                        continue;
                    }
                    Vector2d1 translatedPolygon = translatePolygon(finalPolygon, 0.0, -bottom_distance);
                    judge = doPolygonsCollide2(translatedPolygon, candidate.polygons);
                    if (judge == true) {
                        pymin += bottom_distance;
                        bottom_distance /= 2.0;
                    }
                    else {
                        finalPolygon = translatedPolygon;
                    }
                }
                // �����µĺ�ѡ�������
                std::vector<Vector2d1> newPolygons = candidate.polygons;
                std::vector<int> temp = candidate.typenum;
                temp.push_back(i);//ѹ���¶�������
                newPolygons.push_back(finalPolygon);

                // �����µĽ�������ĵ÷�
                pair<double, int> sc_pv = calculateScore(newPolygons, candidate.previous);
                double newScore = sc_pv.first;

                // �������
                cout << "����id" << id << ":" << newScore << endl;
                // ����ͼ�񲢼�¼�÷�
                geometry_layer_save(newPolygons, id, newScore);
                score.push_back(newScore);
                process_solutions.push_back(newPolygons);
                // ���µĽ�����������ѡ����
                gml_tree.insert(candidate.CandidateId, id);//gml���Ĳ���
                Candidate son(id++, newPolygons, newScore, temp, sc_pv.second);
                nextCandidates.push(son);
                // ���ֺ�ѡ���еĴ�С�����������
                while (nextCandidates.size() > beamWidth) {
                    nextCandidates.pop();
                }
            }
        }
        // ���º�ѡ�����������
        candidates = nextCandidates;
    }

    // ��ȡ��ѵĺ�ѡ�������
    while (candidates.size() > 1) {
        candidates.pop();
    }
    // ����������
    cout << "���score" << candidates.top().score << endl;
    vector<Vector2d1> a = origin_polygons;
    // ������ѽ��������ͼ��
    if (!candidates.top().polygons.empty()) {
        // ��ȡ��ѽ�������Ķ���κ�������ԭʼ�����е�����˳��
        vector<Vector2d1> final_plan = candidates.top().polygons;
        vector<int> final_nums = candidates.top().typenum;
        int i = 0;
        // ������ѽ�������е�ÿ�������
        for (auto it = final_plan.begin(); it != final_plan.end(); it++) {
            // ��ȡ��ǰ�������ԭʼ�����е�����
            int sort_index = final_nums[i];
            i++;
            // ���㵱ǰ������� x �� y ���ϵ�λ����
            double delta_x = ((*it).begin())->x - (sortedPolygons[sort_index].begin())->x;
            double delta_y = ((*it).begin())->y - (sortedPolygons[sort_index].begin())->y;
            cout << "deltas" << delta_x << " " << delta_y << endl;
            int j = 0;
            // ����ԭʼ����Ķ����
            for (auto ooo : ori_Polygons) {
                // �����ǰ������뵱ǰ������ԭʼ�������ͬһ�������
                if (ooo == sortedPolygons[sort_index]) {
                    // �Ըö���ε�ÿ���������λ�ƣ�ʹ���뵱ǰ����ε�λ�ö���
                    for (auto it = a[j].begin(); it != a[j].end(); it++) {
                        (*it).x += delta_x;
                        (*it).y += delta_y;
                    }
                }
                else {
                    j++;
                }
            }
        }
        // ������������ѽ��������ԭʼ�������ε�ͼ��
        geometry_layer_save1(final_plan, a);
    }
    // ������ѽ������
    return candidates.top().polygons;
}

void create_folder(string a) {//���ߺ����������ļ���
    string folderPath = "./" + a;
    CreateDirectory(folderPath.c_str(), NULL);
    return;
}

void Beamsearch::work() {
    string output_filename = image_path;//��Ŀǰ����£�������Ҫ��packing�仯���̵�ÿһ��ͼ����������image_path��һ���ļ��е�ַ��������Beamsearch����
    create_folder(output_filename);//�����������ͼ����ļ���
    SYSTEMTIME st;//��ȡʱ�䣬�ԶԱ����packing����ͼ����и���
    GetSystemTime(&st);
    string time_path = image_path + "/" + to_string(st.wYear) + "_" + to_string(st.wMonth) + "_" + to_string(st.wDay) + "_" + to_string(st.wHour) + "_" + to_string(st.wMinute) + "_" + to_string(st.wSecond);//�˴�packing���̵õ����ļ�����
    create_folder(time_path);//һ��packing��һ���ļ���
    this->image_path = "./" + time_path;
    get_points_to_polygon();//�����ļ����ڵ�Ԫ���ļ�
    //PolygonModification();//��Ԫ�����д��ϻ�
    Vector2d1 boundingRect;//Բ������2ά����ľ���
    int beamWidth = 10;//�����beamsearch�㷨������
    boundingRect.push_back(Vector2d(0, 0));
    boundingRect.push_back(Vector2d(boxx / 2, 0));
    boundingRect.push_back(Vector2d(boxx / 2, boxy));
    boundingRect.push_back(Vector2d(0, boxy));
    std::vector<Vector2d1> wtf = perior_geometry_put();//���ǵ�test.txt�е�Ԫ��������ͬ���Ǹ����ϣ����Ҫ������ͬ��Ԫ������ʹ�øú����ظ������ӦԪ��
    std::vector<Vector2d1> bestSolution = beamSearch(wtf, beamWidth, boundingRect);//�����㷨��beamsearch�㷨

    // �����ѽ������
    std::cout << "��ѽ��������" << std::endl;
    for (const Vector2d1& polygon : bestSolution) {
        // �������ε�����
        for (const Vector2d& point : polygon) {
            std::cout << "(" << point.x << ", " << point.y << ") ";
        }
        std::cout << std::endl;
    }
    if (bestSolution.empty())cout << "�޷����ɽ��������" << endl;//�������Ϊ�գ�������޳ɹ���������˵���⼸��ԭ������ô���ö��ᷢ����ײ��ͻ
    else geometry_layer_output(bestSolution);//�����������
}

void Beamsearch::test() {//���Ժ��������ڲ��Եľ���PolygonModification2()����������кܴ�����⣬������Ŀ�ͽ�չ������
    get_points_to_polygon();
    vector<Vector2d1> a = polygons;
    PolygonModification2();
    geometry_layer_save1(a, polygons);
}

void Beamsearch::get_points_to_polygon() {
    boxx = 700;
    boxy = 1000;
    string address = "test.txt";
    ifstream infile;
    infile.open(address);
    if (!infile.is_open()) {
        std::cout << "�ļ���ʧ��" << endl;
        return;
    }
    string line;
    while (getline(infile, line)) {//ÿ�δ��ļ���ȡһ��
        istringstream iss(line);
        Vector2d1 points;
        int n;
        iss >> n;
        double x, y;
        for (int i = 0; i < n; i++)
        {
            iss >> x >> y;
            points.push_back(Vector2d(x, y));
        }
        if (PL().HGP_2D_Polygon_Is_Clockwise_Oriented_C(points))//��ֹ���˳��ߵ����������ɵĶ�����Ǹ��ģ�����α���������Ϊ������ģ��������ǵÿ���˳��ʱ�����⣩
        {
            std::reverse(points.begin(), points.end());
        }
        polygons.push_back(points);
    }
    origin_polygons = polygons;//������ڴ����е��ʼ�Ķ������
}

std::vector<Vector2d1> Beamsearch::perior_geometry_put()//��������ظ���Ԫ��
{
    std::vector<Vector2d1> ans = polygons;
    std::cout << "ͼ����������Ƿ��ظ�" << endl;
    bool ques;
    cin >> ques;
    if (ques) {
        std::cout << "�ظ����м���" << endl;
        int n; cin >> n;
        while (n--) {
            int type;
            cin >> type;
            if (type > polygons.size()) {
                std::cout << "û�и����͵ļ��νṹŶ������������" << endl;
                n++;
                continue;
            }
            Vector2d1 temp = polygons[type - 1];
            ans.push_back(temp);
        }
    }
    return ans;
}

void Beamsearch::geometry_layer_output(vector<Vector2d1> a) {
    // ����ͼ��ĳߴ�

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_64FC3����ʼֵΪ��ɫ
    cv::Mat rightimage(boxy, boxx / 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // ���ƶ����
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // ���ƶ��������
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }
    /*
    for (const auto& polygon : a) {
        int cnt = polygon.size();
        for (int i = 0; i < cnt; ++i) {
            auto& vertex1 = polygon[i];
            auto& vertex2 = polygon[(i + 1) % cnt];

            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x1 = vertex1.x;
            int y1 = boxy - vertex1.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
    }
    */
    // ���ҷ�תͼ��
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // ����һ����ֱ��
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // ��ʾͼ��
    cv::imshow("Polygons", symmetric_image);
    cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_output2(vector<Vector2d1> a, vector<Vector2d1> b) {
    // ����ͼ��ĳߴ�
    static int count = 0;

    string path = this->image_path;
    path = path + "output" + to_string(count) + ".jpg";
    cout << path << endl;
    ++count;
    
    int delta = 0;
    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_64FC3����ʼֵΪ��ɫ
    cv::Mat rightimage(boxy + delta * 2, boxx / 2 + delta * 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // ���ƶ����
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x + delta * 1;
            int y = boxy - vertex.y + delta; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // ���ƶ��������
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }

    // ���ƶ����
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x + delta * 1;
            int y = boxy - vertex.y + delta; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // ���ƶ��������
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 0, 0), 2);
    }

    // ���ҷ�תͼ��
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // ����һ����ֱ��
    cv::Point point1(boxx / 2 + delta * 2, boxy + delta * 2);
    cv::Point point2(boxx / 2 + delta * 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // ��ʾͼ��
    //cv::imshow("Polygons", symmetric_image);
    cv::imwrite(path, symmetric_image);
    //cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_save(vector<Vector2d1> a, int num, double score) {
    // ����ͼ��ĳߴ�
    string path = this->image_path;
    path = path + "/�ڵ�" + to_string(num) + "���֣�" + to_string(score) + ".jpg";
    cout << path << endl;

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_64FC3����ʼֵΪ��ɫ
    cv::Mat rightimage(boxy, boxx / 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // ���ƶ����
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // ���ƶ��������
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }

    // ���ҷ�תͼ��
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // ����һ����ֱ��
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // ����ͼ��
    cv::imwrite(path, symmetric_image);
    cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_save1(vector<Vector2d1> a, vector<Vector2d1> b) {
    // ����ͼ��ĳߴ�
    string path = this->image_path;
    string path1 = path + "/Roughing.jpg";
    string path2 = path + "/Finishing.jpg";
    string path3 = path + "/BothOfThem.jpg";

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_8UC3����ʼֵΪ��ɫ
    cv::Mat rightimage1(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));
    std::vector<std::vector<cv::Point>> pts;

    // ���ƶ���Σ�ʹ�ú�ɫ���
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts.push_back(points);
    }
    // �������
    cv::fillPoly(rightimage1, pts, cv::Scalar(0, 0, 0));

    // ���ҷ�תͼ��
    cv::Mat leftimage1;
    cv::flip(rightimage1, leftimage1, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image1;
    cv::hconcat(leftimage1, rightimage1, symmetric_image1);

    // ����һ����ֱ��
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image1, point1, point2, cv::Scalar(0, 0, 255), 1);

    // ����ͼ��
    cv::imwrite(path1, symmetric_image1);

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_8UC3����ʼֵΪ��ɫ
    cv::Mat rightimage2(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));
    std::vector<std::vector<cv::Point>> pts1;

    // ���ƶ���Σ�ʹ����ɫ���
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts1.push_back(points);
    }
    // �������
    cv::fillPoly(rightimage2, pts1, cv::Scalar(255, 0, 0));

    // ���ҷ�תͼ��
    cv::Mat leftimage2;
    cv::flip(rightimage2, leftimage2, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image2;
    cv::hconcat(leftimage2, rightimage2, symmetric_image2);

    // ����һ����ֱ��
    cv::Point point3(boxx / 2, boxy);
    cv::Point point4(boxx / 2, 0);
    cv::line(symmetric_image2, point3, point4, cv::Scalar(0, 0, 255), 1);

    // ����ͼ��
    cv::imwrite(path2, symmetric_image2);

    // ����һ����ɫ��ͼ�񣬳ߴ�Ϊ(boxy, boxx / 2)����������ΪCV_8UC3����ʼֵΪ��ɫ
    cv::Mat rightimage3(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));

    std::vector<std::vector<cv::Point>> pts2;

    // ���ƶ���Σ�ʹ�ú�ɫ���
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts2.push_back(points);
    }
    // ������Σ�ʹ�ú�ɫ���
    cv::fillPoly(rightimage3, pts2, cv::Scalar(0, 0, 0));

    std::vector<std::vector<cv::Point>> pts3;

    // ���ƶ���Σ�ʹ����ɫ���
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // ����������ת��ΪOpenCVͼ������ϵ�е�����
            int x = vertex.x;
            int y = boxy - vertex.y; // ��OpenCV�У�ͼ���ԭ��λ�����Ͻǣ�������Ҫ��תy��
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts3.push_back(points);
    }
    // ������Σ�ʹ����ɫ���
    cv::fillPoly(rightimage3, pts3, cv::Scalar(0, 0, 255));

    // ���ҷ�תͼ��
    cv::Mat leftimage3;
    cv::flip(rightimage3, leftimage3, 1);

    // ƴ������ͼ�񣬵õ��Գ�ͼ��
    cv::Mat symmetric_image3;
    cv::hconcat(leftimage3, rightimage3, symmetric_image3);

    // ����һ����ֱ��
    cv::Point point5(boxx / 2, boxy);
    cv::Point point6(boxx / 2, 0);
    cv::line(symmetric_image3, point5, point6, cv::Scalar(0, 0, 255), 1);

    // ����ͼ��
    cv::imwrite(path3, symmetric_image3);
    return;
}

vector<double> Beamsearch::GetScore()
{
    return score;
}

void Beamsearch::PolygonModification() {

}

void Beamsearch::PolygonModification1()
{
}

Polygon_2 convex_output(std::vector<Point_2> points)
{
    Polygon_2 convex_hull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull));
    return convex_hull;
}

void findMaxTriangleArea(std::vector<Point>& points, std::vector<Point>& all_points) {
    //ͼ�δ��ϻ�������Ҫ�������ߺ�����ͼ�εĵ���Ϣ����˳������������ɢȡ�㣬�ҵ����Ž⣬��û�кϷ��㣬�ӳ����ߣ����֣�ĩ�����Χ���ཻ���õ��µ����߽��м���
    double lenth = 0;
    Point f, t;
    auto num2 = points.begin();
    auto num3 = points.end();
    num2++;
    num2++;
    num3--;
    //��¼��β����Ϣ
    f = *points.begin();
    t = *(num3);
    Point final_point;//���յ�
    double final_area = -1;//�������������
    for (auto i = num2; i != num3; i++) {//�����ܳ�
        Point point1, point2;
        point1 = *i;
        i--;
        point2 = *i;
        lenth += sqrt(CGAL::squared_distance(point1, point2));
        i++;
    }
    lenth /= 100;//�ܳ���1/100��Ϊ����
    int f_num = 0, t_num = 0;//��¼��β�����߲��Ϸ��Ĵ���������ѡ���ٵ�һ���ӳ�
    for (auto i = num2; i != num3; i++) {
        double x = (*i).x();
        double y = (*i).y();
        Point point1, point2;
        point1 = *i;
        i--;
        point2 = *i;
        Segment_2 segment(point2, point1);//�жϸñ��ϵĵ��Ƿ�Ϸ�
        x -= (*i).x();
        y -= (*i).y();
        Point en = *i;
        i++;
        double div = sqrt(x * x + y * y);
        x /= div;
        y /= div;
        x *= lenth;
        y *= lenth;
        while (CGAL::squared_distance(segment, en) == 0) {//���ڱ���
            //�����Ƿ�û����ײ�����Ƿ�Ϸ�
            Segment ray1(f, en);
            Segment ray2(t, en);
            bool is = true;
            for (auto j = num2; j != num3; j++) {
                Point pointj1, pointj2;
                pointj1 = *j;
                j--;
                pointj2 = *j;
                j++;
                Segment_2 seg(pointj2, pointj1);
                if (seg == segment)continue;
                if (CGAL::do_intersect(ray1, seg)) {
                    f_num++;
                    is = false;
                    continue;
                }
            }
            for (auto j = num2; j != num3; j++) {
                Point pointj1, pointj2;
                pointj1 = *j;
                j--;
                pointj2 = *j;
                j++;
                Segment_2 seg(pointj2, pointj1);
                if (seg == segment)continue;
                if (CGAL::do_intersect(ray2, seg)) {
                    t_num++;
                    is = false;
                    continue;
                }
            }
            if (is) {
                Triangle triangle(f, t, en);
                double area = abs(triangle.area());
                if (area > final_area) {
                    final_point = en;
                    final_area = area;
                }
            }
            en = Point(en.x() + x, en.y() + y);
        }
    }
    if (final_area != -1) {//�ҵ����Ž�
        vector<Point> ans;
        ans.push_back(f);
        ans.push_back(final_point);
        ans.push_back(t);
        points = ans;
        return;
    }
    else {//û�ҵ����Ž⣬��Ҫ�ӳ�����
        Point f_next = *(num2);
        auto tt = points.end();
        tt -= 2;
        Point t_next = *(tt);

        Ray ray1(f_next, f);
        Ray ray2(t_next, t);
        Polygon_2 the_convex_polygon = convex_output(all_points);//��Χ��
        std::vector<Point_2> intersectionPointsRay1;
        Point f_inter_point;
        Point t_inter_point;
        for (Polygon_2::Edge_const_iterator itt = the_convex_polygon.edges_begin(); itt != the_convex_polygon.edges_end(); ++itt)
        {//�������Χ�н���
            Segment seg = *itt;
            auto result1 = CGAL::intersection(seg, ray1);
            auto result2 = CGAL::intersection(seg, ray2);
            if (result1) {
                if (const Point_2* ipoint = boost::get<Point_2>(&*result1)) {
                    f_inter_point = *ipoint;
                }
            }
            if (result2) {
                if (const Point_2* ipoint = boost::get<Point_2>(&*result2)) {
                    t_inter_point = *ipoint;
                }
            }
        }
        if (f_num > t_num) {
            double fx = (f_inter_point).x();
            double fy = (f_inter_point).y();
            fx -= f.x();
            fy -= f.y();
            fx /= 7;//��Ϊ�߷ݣ�����ֱ�ӵ����㣬����һ������������
            fy /= 7;
            Point f_new(f.x() + fx, f.y() + fy);
            while (CGAL::squared_distance(Segment(f, f_inter_point), f_new) == 0) {//���������ٴ�Ѱ��
                for (auto i = num2; i != num3; i++) {
                    double x = (*i).x();
                    double y = (*i).y();
                    Point point1, point2;
                    point1 = *i;
                    i--;
                    point2 = *i;
                    Segment_2 segment(point2, point1);
                    x -= (*i).x();
                    y -= (*i).y();
                    Point en = *i;
                    i++;
                    double div = sqrt(x * x + y * y);
                    x /= div;
                    y /= div;
                    x *= lenth;
                    y *= lenth;
                    while (CGAL::squared_distance(segment, en) == 0) {
                        en = Point(en.x() + x, en.y() + y);
                        Segment ray1(f_new, en);
                        Segment ray2(t, en);
                        bool is = true;
                        for (auto j = num2; j != num3; j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray1, seg)) {
                                f_num++;
                                is = false;
                                break;
                            }
                        }
                        for (auto j = num2; j != num3; j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray2, seg)) {
                                t_num++;
                                is = false;
                                break;
                            }
                        }
                        if (is) {
                            Triangle triangle(f, t, en);
                            double area = abs(triangle.area());
                            if (area > final_area) {
                                final_point = en;
                                final_area = area;
                            }
                        }
                    }
                }
                if (final_area != -1) {
                    vector<Point> ans;
                    ans.push_back(f_new);
                    ans.push_back(final_point);
                    ans.push_back(t);
                    points = ans;
                    return;
                }
                f_new = Point(f_new.x() + fx, f_new.y() + fy);
            }
            f = f_inter_point;
            double tx = (t_inter_point).x();
            double ty = (t_inter_point).y();
            tx -= t.x();
            ty -= t.y();
            tx /= 7;
            ty /= 7;
            Point t_new(t.x() + tx, t.y() + ty);
            while (CGAL::squared_distance(Segment(t, t_inter_point), t_new) == 0) {
                for (auto i = num2; i != num3; i++) {
                    double x = (*i).x();
                    double y = (*i).y();
                    Point point1, point2;
                    point1 = *i;
                    i--;
                    point2 = *i;
                    Segment_2 segment(point2, point1);
                    x -= (*i).x();
                    y -= (*i).y();
                    Point en = *i;
                    i++;
                    double div = sqrt(x * x * +y * y);
                    x /= div;
                    y /= div;
                    x *= lenth;
                    y *= lenth;
                    while (CGAL::squared_distance(segment, en) == 0) {
                        en = Point(en.x() + x, en.y() + y);
                        Segment ray1(f, en);
                        Segment ray2(t_new, en);
                        bool is = true;
                        for (auto j = num2; j != points.end(); j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray1, seg)) {
                                f_num++;
                                is = false;
                                break;
                            }
                        }
                        for (auto j = num2; j != points.end(); j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray2, seg)) {
                                t_num++;
                                is = false;
                                break;
                            }
                        }
                        if (is) {
                            Triangle triangle(f, t, en);
                            double area = abs(triangle.area());
                            if (area > final_area) {
                                final_point = en;
                                final_area = area;
                            }
                        }
                    }
                }
                if (final_area != -1) {
                    vector<Point> ans;
                    ans.push_back(f);
                    ans.push_back(final_point);
                    ans.push_back(t_new);
                    points = ans;
                    return;
                }
                t_new = Point(t_new.x() + tx, t_new.y() + ty);
            }
        }
        else {//�յ����
            double tx = (t_inter_point).x();
            double ty = (t_inter_point).y();
            tx -= t.x();
            ty -= t.y();
            tx /= 7;
            ty /= 7;
            Point t_new(t.x() + tx, t.y() + ty);
            while (CGAL::squared_distance(Segment(t, t_inter_point), t_new) == 0) {
                for (auto i = num2; i != num3; i++) {
                    double x = (*i).x();
                    double y = (*i).y();
                    Segment_2 segment(*i, *(--i));
                    x -= (*i).x();
                    y -= (*i).y();
                    Point en = *i;
                    i++;
                    double div = sqrt(x * x * +y * y);
                    x /= div;
                    y /= div;
                    x *= lenth;
                    y *= lenth;
                    while (CGAL::squared_distance(segment, en) == 0) {
                        en = Point(en.x() + x, en.y() + y);
                        Segment ray1(f, en);
                        Segment ray2(t_new, en);
                        bool is = true;
                        for (auto j = num2; j != num3; j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray1, seg)) {
                                f_num++;
                                is = false;
                                break;
                            }
                        }
                        for (auto j = num2; j != num3; j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray2, seg)) {
                                t_num++;
                                is = false;
                                break;
                            }
                        }
                        if (is) {
                            Triangle triangle(f, t, en);
                            double area = abs(triangle.area());
                            if (area > final_area) {
                                final_point = en;
                                final_area = area;
                            }
                        }
                    }
                }
                if (final_area != -1) {
                    vector<Point> ans;
                    ans.push_back(f);
                    ans.push_back(final_point);
                    ans.push_back(t_new);
                    points = ans;
                    return;
                }
                t_new = Point(t_new.x() + tx, t_new.y() + ty);
            }
            t = t_inter_point;
            double fx = (f_inter_point).x();
            double fy = (f_inter_point).y();
            fx -= f.x();
            fy -= f.y();
            fx /= 7;
            fy /= 7;
            Point f_new(f.x() + fx, f.y() + fy);
            while (CGAL::squared_distance(Segment(f, f_inter_point), f_new) == 0) {
                for (auto i = num2; i != num3; i++) {
                    double x = (*i).x();
                    double y = (*i).y();
                    Segment_2 segment(*i, *(--i));
                    x -= (*i).x();
                    y -= (*i).y();
                    Point en = *i;
                    i++;
                    double div = sqrt(x * x * +y * y);
                    x /= div;
                    y /= div;
                    x *= lenth;
                    y *= lenth;
                    while (CGAL::squared_distance(segment, en) == 0) {
                        en = Point(en.x() + x, en.y() + y);
                        Segment ray1(f_new, en);
                        Segment ray2(t, en);
                        bool is = true;
                        for (auto j = num2; j != points.end(); j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray1, seg)) {
                                f_num++;
                                is = false;
                                break;
                            }
                        }
                        for (auto j = num2; j != points.end(); j++) {
                            Point pointj1, pointj2;
                            pointj1 = *j;
                            j--;
                            pointj2 = *j;
                            j++;
                            Segment_2 seg(pointj2, pointj1);
                            if (seg == segment)continue;
                            if (CGAL::do_intersect(ray2, seg)) {
                                t_num++;
                                is = false;
                                break;
                            }
                        }
                        if (is) {
                            Triangle triangle(f, t, en);
                            double area = abs(triangle.area());
                            if (area > final_area) {
                                final_point = en;
                                final_area = area;
                            }
                        }
                    }
                }
                if (final_area != -1) {
                    vector<Point> ans;
                    ans.push_back(f_new);
                    ans.push_back(final_point);
                    ans.push_back(t);
                    points = ans;
                    return;
                }
                f_new = Point(f_new.x() + fx, f_new.y() + fy);
            }

        }
    }
    if (final_area == -1) {
        cout << "����";
    }
}

void Beamsearch::PolygonModification2() {

    //ʹ��findMaxTriangleArea����±���Ƕ±�
    vector<Vector2d1> newpolygons;
    for (auto it = polygons.begin(); it != polygons.end(); it++)
    {
        Polygon_2 wewant = Convert_Vector2d1_to_Polygon_2(*it);
        vector<Point> poly_points;
        for (Polygon_2::Edge_const_iterator itt = wewant.edges_begin(); itt != wewant.edges_end(); ++itt) {
            poly_points.push_back(itt->source());
        }
        bool stopit = 0;
        while (1)
        {
            vector<Point> all_points;
            for (Polygon_2::Edge_const_iterator itt = wewant.edges_begin(); itt != wewant.edges_end(); ++itt) {
                all_points.push_back(itt->source());
            }
            vector<pair<Segment_2, int>> MarkEdges;
            for (Polygon_2::Edge_const_iterator itt = wewant.edges_begin(); itt != wewant.edges_end(); ++itt) {
                Point_2 start = itt->source();
                Point_2 end = itt->target();
                Ray ray1(start, end);
                Ray ray2(end, start);
                int nums1 = 0, nums2 = 0;
                for (Polygon_2::Edge_const_iterator w = wewant.edges_begin(); w != wewant.edges_end(); ++w) {
                    if ((w->source() == start && w->target() == end) || (w->source() == end && w->target() == start))continue;
                    if (CGAL::do_intersect(ray1, *w))nums1++;
                    if (CGAL::do_intersect(ray2, *w))nums2++;

                }
                pair<Segment_2, int> MarkEdge;
                MarkEdge.first = (*itt);
                MarkEdge.second = 0;
                if ((nums1 >= 3) && (nums2 >= 3)) {
                    MarkEdge.second = 2;
                }
                else if (nums1 == 2 && nums2 == 2)
                {
                    MarkEdge.second = 0;
                }
                else {
                    MarkEdge.second = 1;
                }
                MarkEdges.push_back(MarkEdge);
            }
            stopit = 1;
            for (int i = 0; i < MarkEdges.size(); i++)
            {
                if (MarkEdges[i].second == 2)
                {
                    stopit = 0;
                    break;
                }
            }
            if (stopit == 1)break;
            for (int i = 0; i < MarkEdges.size(); i++)
            {
                int size = MarkEdges.size();
                int tab = MarkEdges[i].second;
                if (tab == 2)
                {
                    int first = (i - 1 + size) % size;
                    while (true)
                    {
                        assert(MarkEdges[first].second != 0);
                        if (MarkEdges[first].second == 1)break;
                        first = (first - 1 + size) % size;
                    }
                    Segment first_edge = MarkEdges[first].first;
                    Point point1 = first_edge.source();

                    int last = (i + 1) % size;
                    while (true)
                    {
                        assert(MarkEdges[last].second != 0);
                        if (MarkEdges[last].second == 1)break;
                        last = (last + 1) % size;
                    }
                    Segment last_edge = MarkEdges[last].first;
                    Point point2 = last_edge.target();
                    vector<Point> points;//your target point set
                    points.push_back(point1);
                    int start = first + 1;
                    while (1) {
                        Point temp_point;
                        Segment newedge = MarkEdges[start].first;
                        temp_point = newedge.source();
                        points.push_back(temp_point);
                        start++;
                        if (newedge.target() == point2)
                        {
                            points.push_back(point2);
                            break;
                        }
                    }

                    vector<Point> ans;
                    vector<Point> origin_points = points;
                    cout << "����һ���±߼���" << endl;
                    for (auto pp = points.begin(); pp != points.end(); pp++)cout << (*pp) << "||||";
                    cout << endl;
                    cout << "��ǰȫ�㼯��" << endl;
                    for (auto pp = all_points.begin(); pp != all_points.end(); pp++)cout << (*pp) << "||||";
                    cout << endl;
                    Polygon_2 pf(all_points.begin(), all_points.end());
                    vector<Vector2d1> pfv;
                    pfv.push_back(Convert_Polygon_2_to_Vector2d1(pf));

                    outputLayer outL(boxx, boxy);
                    PolygonWithCutDir poly;
                    poly.polygon = Convert_Polygon_2_to_Vector2d1(pf);
                    poly.addAntiClockWiseCut();
                    outL.addPolygon(poly);
                    
                    outL.geometry_layer_output();

                    //geometry_layer_output(pfv);


                    findMaxTriangleArea(points, all_points);
                    cout << "�����ǰ�Ĺ���Ķ±߼�" << endl;
                    for (auto pp = points.begin(); pp != points.end(); pp++)cout << (*pp) << "||||";
                    cout << endl;


                    for (int gg = 0; gg < all_points.size(); gg++)
                    {
                        bool is = 0, isfirst = 0;
                        for (int hh = 0; hh < origin_points.size(); hh++)
                        {
                            if (all_points[gg] == origin_points[hh])is = 1;
                        }
                        if (is == 0)ans.push_back(all_points[gg]);
                        if (all_points[gg] == point1) {
                            for (int hh = 0; hh < points.size(); hh++)
                            {
                                ans.push_back(points[hh]);
                            }
                        }
                    }
                    cout << "�����ǰ�Ĺ����ȫ�㼯" << endl;
                    for (auto pp = ans.begin(); pp != ans.end(); pp++)cout << (*pp) << "||||";
                    cout << endl;

                    Polygon_2 newans(ans.begin(), ans.end());
                    vector<Vector2d1> plv;
                    plv.push_back(Convert_Polygon_2_to_Vector2d1(newans));
                    //geometry_layer_output2(plv, pfv);
                    wewant = newans;
                    break;
                }
            }
        }
        newpolygons.push_back(Convert_Polygon_2_to_Vector2d1(wewant));
        Polygon_2 pf(poly_points.begin(), poly_points.end());
        vector<Vector2d1> pfv;
        pfv.push_back(Convert_Polygon_2_to_Vector2d1(pf));
        geometry_layer_output2(pfv, newpolygons);
    }
    polygons = newpolygons;
}
