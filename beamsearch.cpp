#include "beamsearch.h"
double boxx, boxy;
// 定义用于多边形排序的比较函数
bool beamssort(const Vector2d1& poly1, const Vector2d1& poly2)
{
    // 定义变量存储多边形的特征值
    double xmax_score1, xmin_score1, xlength_score1, ylength_score1, area_score1;
    double xmax_score2, xmin_score2, xlength_score2, ylength_score2, area_score2;
    Vector2d max1, min1, max2, min2;
    // 计算第一个多边形的边界框特征值
    PL().HGP_2d_Polygon_Boundingbox_C(poly1, min1, max1);
    xmax_score1 = max1.x / boxx;
    xmin_score1 = min1.x / boxx;
    xlength_score1 = (max1.x - min1.x) / boxx;
    ylength_score1 = (max1.y - min1.y) / boxy;
    // 计算第二个多边形的边界框特征值
    PL().HGP_2d_Polygon_Boundingbox_C(poly2, min2, max2);
    xmax_score2 = max2.x / boxx;
    xmin_score2 = min2.x / boxx;
    xlength_score2 = (max2.x - min2.x) / boxx;
    ylength_score2 = (max2.y - min2.y) / boxy;
    // 计算多边形的面积特征值
    area_score1 = PL().HGP_2D_Polygon_Area_C(poly1) / (boxx * (max1.y - min1.y));
    area_score2 = PL().HGP_2D_Polygon_Area_C(poly2) / (boxx * (max2.y - min2.y));
    // 计算多边形的总得分
    double score1 = xmax_score1 + xmin_score1 + xlength_score1 + ylength_score1 + area_score1;
    double score2 = xmax_score2 + xmin_score2 + xlength_score2 + ylength_score2 + area_score2;
    // 比较两个多边形的总得分
    return score1 > score2;
}

double pointToPolygonDist(const Point_2& p, const Polygon_2& polygon) {//点到多边形的距离;检验完成，没报错
    int count = 0;
    double minDist = INFINITY;

    for (auto e = polygon.edges_begin(); e != polygon.edges_end(); ++e) {
        const Point_2& a = e->source();
        const Point_2& b = e->target();

        // 判断点 p 与线段 ab 是否在同一水平线上，并且 p 在 ab 的左侧
        if ((a.y() > p.y() != b.y() > p.y()) &&
            (p.x() < (b.x() - a.x()) * (p.y() - a.y()) / (b.y() - a.y()) + a.x())) {
            count++;
        }

        // 计算点 p 到线段 ab 的最短距离，并更新最小距离
        Segment_2 s(a, b);
        minDist = std::min(minDist, squared_distance(p, s));
    }

    if (count % 2 == 0) {
        return std::sqrt(abs(minDist));

    }
    return -std::sqrt(abs(minDist));
}

void CGAL_2D_Polygon_Dart_Sampling_a(const vector<Polygon_2>& py, const double& d, vector<Point_2>& sampling_points)
{//离散取点方式：均匀取点
    Functs::MAssert(d > 0 && d < 1.0, "CGAL_2D_Polygon_Dart_Sampling if (!(d > 0 && d < 1.0))");

    double xmin = 0;
    double ymin = 0;
    double xmax = boxx;
    double ymax = boxy;
    double minimal_d = 15;
    int total_iter = 1 / d;
    double dx = xmax / total_iter;
    double dy = ymax / total_iter;
    vector<Point_2> insert_points;
    for (int i = 0; i < total_iter; i++)
    {
        double x = xmin + dx * i;
        for (int j = 0; j < total_iter; j++) {

            double y = ymin + dy * j;
            Point_2 point(x, y);
            double distance = CGAL_IA_MAX_DOUBLE;
            for (int num = 0; num < py.size(); num++) {
                Polygon_2 poly = py[num];
                distance = min(distance, pointToPolygonDist(point, poly));

            }
            if (distance > minimal_d)
            {
                insert_points.push_back(Point_2(x, y));
            }
        }
    }
    sampling_points = insert_points;
}

void findPolygons(vector<Segment_2> edges, vector<vector<Segment_2>>& polys) {
    //将无序的边连成图形，可返回多个图形
    while (!edges.empty()) {
        vector<Segment_2> poly;
        Segment_2 currentEdge = edges[0];
        poly.push_back(currentEdge);
        edges.erase(edges.begin());
        while (true) {
            Segment_2 nextEdge;
            bool found = false;
            for (auto it = edges.begin(); it != edges.end(); it++) {
                if (it->source().x() == currentEdge.target().x() && it->source().y() == currentEdge.target().y())
                {
                    nextEdge = *it;
                    edges.erase(it);
                    found = true;
                    break;
                }
                else if (it->target().x() == currentEdge.target().x() && it->target().y() == currentEdge.target().y())
                {
                    Segment_2 wewant(Point_2(it->target().x(), it->target().y()), Point_2(it->source().x(), it->source().y()));
                    nextEdge = wewant;
                    edges.erase(it);
                    found = true;
                    break;
                }
            }
            if (found) {
                poly.push_back(nextEdge);
                currentEdge = nextEdge;
            }
            else {
                break;
            }
        }
        polys.push_back(poly);
    }
}

Polygon_2 Convert_Vector2d1_to_Polygon_2(Vector2d1 v2d)//类型转换Vector2d1转化为Polygon
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

Vector2d1 Convert_Polygon_2_to_Vector2d1(Polygon_2 p2)//类型转换Polygon转化为Vector2d1
{
    Vector2d1 target;
    for (auto itt = p2.begin(); itt != p2.end(); itt++)
    {
        Vector2d temp((*itt).x(), (*itt).y());
        target.push_back(temp);
    }
    return target;
}
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

bool comPoints(const Point_2& p1, const Point_2& p2) {
    if (p1.x() != p2.x()) {
        return p1.x() < p2.x();
    }
    return p1.y() < p2.y();
}

vector<Polygon_2> get_triangulation_net(vector<Point_2> point_set, vector<Polygon_2> py)
{//使用Delaunay生成网格图，去除过长的边，找到连续的边缘点，即可返回图形
    Delaunay triangulation;
    triangulation.insert(point_set.begin(), point_set.end());//生成网格图
    double maxlength = 400;
    vector<Segment_2> one_face_edges;
    map<Segment_2, int, SegmentComparator> edge_map;
    for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); ++fit) {
        // 遍历三角形的三条边.
        bool tab = false;
        Delaunay::Vertex_handle v1 = fit->vertex(0);
        Delaunay::Vertex_handle v2 = fit->vertex(1);
        Delaunay::Vertex_handle v3 = fit->vertex(2);
        Point_2 wt1 = v1->point();
        Point_2 wt2 = v2->point();
        Point_2 wt3 = v3->point();
        Segment_2 edge1(wt2, wt1);
        Segment_2 edge2(wt3, wt1);
        Segment_2 edge3(wt3, wt2);
        if (comPoints(wt1, wt2)) {
            Segment_2 edge(wt1, wt2);
            edge1 = edge;
        }
        else {
            Segment_2 edge(wt2, wt1);
            edge1 = edge;
        }

        if (comPoints(wt2, wt3)) {
            Segment_2 edge(wt2, wt3);
            edge2 = edge;
        }
        else {
            Segment_2 edge(wt3, wt2);
            edge2 = edge;
        }

        if (comPoints(wt1, wt3)) {
            Segment_2 edge(wt1, wt3);
            edge3 = edge;
        }
        else {
            Segment_2 edge(wt3, wt1);
            edge3 = edge;
        }
        double edge1_length = CGAL::to_double(edge1.squared_length());
        double edge2_length = CGAL::to_double(edge2.squared_length());
        double edge3_length = CGAL::to_double(edge3.squared_length());
        double maxedge = max(edge1_length, max(edge2_length, edge3_length));
        if (maxedge < maxlength) {
            edge_map[edge1]++;
            edge_map[edge2]++;
            edge_map[edge3]++;
        }
    }
    for (auto it = edge_map.begin(); it != edge_map.end(); it++) {
        if (it->second != 2) {
            one_face_edges.push_back(it->first);
        }
    }
    vector<vector<Segment_2>> ans;
    findPolygons(one_face_edges, ans);
    vector<Polygon_2> true_ans;
    for (auto itt = ans.begin(); itt != ans.end(); itt++)
    {
        vector<Point_2> pots;
        //cv::Mat image(boxy + 200, boxx + 200, CV_8UC3, cv::Scalar(255, 255, 255));
        for (auto it = (*itt).begin(); it != (*itt).end(); it++)
        {
            Segment_2 edge = (*it);
            pots.push_back(Point_2(edge.source().x(), edge.source().y()));
            //cv::line(image, cv::Point(edge[0].x() + 100, edge[0].y() + 100),
                //cv::Point(edge[1].x() + 100, edge[1].y() + 100), cv::Scalar(0, 0, 0), 1);
        }
        Polygon_2 wewant(pots.begin(), pots.end());
        true_ans.push_back(wewant);
        // 显示图像
        //cv::imshow("Image with Point", image);
        //cv::waitKey(0);
    }
    return true_ans;
}

pair<double, int> Beamsearch::calculateScore(const std::vector<Vector2d1>& polygons, int previous)
{
    double areas = 0;
    double score = 0;
    int now = 0;
    vector<Polygon_2> pys;
    vector<Polygon_2> ans;
    vector<Vector2d1> output;
    for (auto it = polygons.begin(); it != polygons.end(); it++)
    {
        pys.push_back(Convert_Vector2d1_to_Polygon_2(*it));
    }
    vector<Point_2> getit;
    CGAL_2D_Polygon_Dart_Sampling_a(pys, 0.002, getit);//离散取点，判断是否在图形外侧，返回点集
    //CGAL_2D_Polygon_Dart_Sampling_b(pys, 0.5, getit, 100);//离散取点，判断是否在图形外侧，返回点集
    ans = get_triangulation_net(getit, pys);//根据点集生成晶胞
    for (auto it = ans.begin(); it != ans.end(); it++) {
        output.push_back(Convert_Polygon_2_to_Vector2d1(*it));
    }
    geometry_layer_output2(output, polygons);
    for (auto it = ans.begin(); it != ans.end(); it++)
    {
        if (abs(it->area()) < 500) {//较小的晶胞不予考虑
            continue;
        }
        else {//处理晶胞
            //这里考虑的是晶胞的数量，包围盒面积与周长，周长
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
    score = score * 0.9 + 0.1 * min(previous / now, 1);//公式，参数可调整
    return make_pair(score, now);//返回得分与晶胞数量
}

bool Beamsearch::doPolygonsCollide2(const Vector2d1& poly1, const vector<Vector2d1>& poly2) {//碰撞检测，多边形求交
    for (const Vector2d1& one_polygon : poly2) {
        if (PL().HGP_2D_Two_Polygons_Intersection_C(poly1, one_polygon) > 0) {
            return true; // 发生碰撞
        }
    }
    return false; // 未发生碰撞
}

Vector2d1 Beamsearch::translatePolygon(const Vector2d1& polygon, double dx, double dy) {
    std::vector<Vector2d> translatedVertices;
    // 遍历所有顶点，对每个顶点进行平移操作，并添加到新的顶点列表中
    for (auto it = polygon.begin(); it != polygon.end(); ++it) {

        Vector2d translatedPoint((*it).x + dx, (*it).y + dy);
        translatedVertices.push_back(translatedPoint);
    }
    // 使用新的顶点列表构造一个新的Vector2d1对象并返回
    return Vector2d1(translatedVertices.begin(), translatedVertices.end());
}

std::vector<Vector2d1> Beamsearch::beamSearch(const std::vector<Vector2d1>& inputPolygons, int beamWidth, const Vector2d1& boundingRect) {
    // 初始化候选解决方案的id
    int id = 0;
    // 创建一个优先队列，用于存储候选解决方案，按照得分从高到低排序
    std::priority_queue <Candidate, std::vector<Candidate>, less<Candidate>> candidates;
    // 创建根节点
    Candidate root(id++, {}, 0.0, {}, 1);//空的节点，没有加入多边形，评分也是0
    // 将根节点加入候选解决方案队列
    candidates.push(root);
    // 插入一个虚拟节点，标志根节点
    gml_tree.insert(-1, 0);//你肯定对这个gml_tree很迷惑，见代码详解文档
    std::vector<Vector2d1> sortedPolygons = inputPolygons;
    // 备份原始多边形
    std::vector<Vector2d1> ori_Polygons = inputPolygons;
    std::sort(sortedPolygons.begin(), sortedPolygons.end(), beamssort);// 对输入的多边形进行排序，按照一定的规则，这里按照beamsort，获得的是多边形排序

    // 处理每个待放置的多边形
    for (int times = 0; times < sortedPolygons.size(); times++) {
        // 用于存储下一轮次的候选解决方案的优先队列
        std::priority_queue < Candidate, std::vector<Candidate>, less<Candidate>> nextCandidates;

        // 处理当前轮次的每个候选解决方案
        while (!candidates.empty()) {
            // 获取当前最优的候选解决方案
            Candidate candidate = candidates.top();
            candidates.pop();
            // 控制放置次数的变量
            int tab = 0;

            // 尝试在当前位置放置不同的多边形
            for (int i = 0; i < sortedPolygons.size(); i++) {
                // 检查该多边形是否已经放置在解决方案中
                bool type_tab = 0;
                for (auto types : candidate.typenum) {
                    if (types == i) {
                        type_tab = 1;
                        break;
                    }
                }
                // 如果多边形已经放置在解决方案中，则跳过
                if (type_tab != 0) {
                    continue;
                }
                // 控制放置次数，最多尝试3次
                if (tab < 3) tab++;
                else break;

                // 将多边形放置在容器的最顶部
                Vector2d bomin, bomax, somin, somax;
                PL().HGP_2d_Polygon_Boundingbox_C(boundingRect, bomin, bomax);
                PL().HGP_2d_Polygon_Boundingbox_C(sortedPolygons[i], somin, somax);
                double dx = 0.0;
                double dy = bomax.y - somax.y;
                Vector2d1 finalPolygon = translatePolygon(sortedPolygons[i], dx, dy);
                // 如果放置后发生碰撞或超出边界，则跳过
                if (doPolygonsCollide2(finalPolygon, candidate.polygons) || somax.y > boxy) {
                    tab--;
                    continue;
                }
                // 使用二分法进行平移，直到发生碰撞
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
                // 生成新的候选解决方案
                std::vector<Vector2d1> newPolygons = candidate.polygons;
                std::vector<int> temp = candidate.typenum;
                temp.push_back(i);//压入新多边形序号
                newPolygons.push_back(finalPolygon);

                // 计算新的解决方案的得分
                pair<double, int> sc_pv = calculateScore(newPolygons, candidate.previous);
                double newScore = sc_pv.first;

                // 输出评分
                cout << "方案id" << id << ":" << newScore << endl;
                // 保存图像并记录得分
                geometry_layer_save(newPolygons, id, newScore);
                score.push_back(newScore);
                process_solutions.push_back(newPolygons);
                // 将新的解决方案加入候选队列
                gml_tree.insert(candidate.CandidateId, id);//gml树的插入
                Candidate son(id++, newPolygons, newScore, temp, sc_pv.second);
                nextCandidates.push(son);
                // 保持候选队列的大小不超过束宽度
                while (nextCandidates.size() > beamWidth) {
                    nextCandidates.pop();
                }
            }
        }
        // 更新候选解决方案队列
        candidates = nextCandidates;
    }

    // 获取最佳的候选解决方案
    while (candidates.size() > 1) {
        candidates.pop();
    }
    // 输出最佳评分
    cout << "最佳score" << candidates.top().score << endl;
    vector<Vector2d1> a = origin_polygons;
    // 保存最佳解决方案的图像
    if (!candidates.top().polygons.empty()) {
        // 获取最佳解决方案的多边形和它们在原始输入中的索引顺序
        vector<Vector2d1> final_plan = candidates.top().polygons;
        vector<int> final_nums = candidates.top().typenum;
        int i = 0;
        // 遍历最佳解决方案中的每个多边形
        for (auto it = final_plan.begin(); it != final_plan.end(); it++) {
            // 获取当前多边形在原始输入中的索引
            int sort_index = final_nums[i];
            i++;
            // 计算当前多边形在 x 和 y 轴上的位移量
            double delta_x = ((*it).begin())->x - (sortedPolygons[sort_index].begin())->x;
            double delta_y = ((*it).begin())->y - (sortedPolygons[sort_index].begin())->y;
            cout << "deltas" << delta_x << " " << delta_y << endl;
            int j = 0;
            // 遍历原始输入的多边形
            for (auto ooo : ori_Polygons) {
                // 如果当前多边形与当前遍历的原始多边形是同一个多边形
                if (ooo == sortedPolygons[sort_index]) {
                    // 对该多边形的每个顶点进行位移，使其与当前多边形的位置对齐
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
        // 保存调整后的最佳解决方案和原始输入多边形的图像
        geometry_layer_save1(final_plan, a);
    }
    // 返回最佳解决方案
    return candidates.top().polygons;
}

void create_folder(string a) {//工具函数，创建文件夹
    string folderPath = "./" + a;
    CreateDirectory(folderPath.c_str(), NULL);
    return;
}

void Beamsearch::work() {
    string output_filename = image_path;//在目前情况下，我们需要将packing变化过程的每一张图保存下来，image_path是一个文件夹地址，定义在Beamsearch类里
    create_folder(output_filename);//创建保存输出图像的文件夹
    SYSTEMTIME st;//获取时间，以对保存的packing过程图像进行赋名
    GetSystemTime(&st);
    string time_path = image_path + "/" + to_string(st.wYear) + "_" + to_string(st.wMonth) + "_" + to_string(st.wDay) + "_" + to_string(st.wHour) + "_" + to_string(st.wMinute) + "_" + to_string(st.wSecond);//此次packing过程得到的文件夹名
    create_folder(time_path);//一次packing，一个文件夹
    this->image_path = "./" + time_path;
    get_points_to_polygon();//导入文件夹内的元件文件
    //PolygonModification();//将元件进行粗料化
    Vector2d1 boundingRect;//圆柱材料2维截面的矩形
    int beamWidth = 10;//这个即beamsearch算法的束宽
    boundingRect.push_back(Vector2d(0, 0));
    boundingRect.push_back(Vector2d(boxx, 0));
    boundingRect.push_back(Vector2d(boxx, boxy));
    boundingRect.push_back(Vector2d(0, boxy));
    std::vector<Vector2d1> wtf = perior_geometry_put();//我们的test.txt中的元件各不相同，是个集合，如果要加入相同的元件，就使用该函数重复加入对应元件
    std::vector<Vector2d1> bestSolution = beamSearch(wtf, beamWidth, boundingRect);//核心算法，beamsearch算法

    // 输出最佳解决方案
    std::cout << "最佳解决方案：" << std::endl;
    for (const Vector2d1& polygon : bestSolution) {
        // 输出多边形的坐标
        for (const Vector2d& point : polygon) {
            std::cout << "(" << point.x << ", " << point.y << ") ";
        }
        std::cout << std::endl;
    }
    if (bestSolution.empty())cout << "无法生成解决方案！" << endl;//如果返回为空，则代表无成功方案，这说明这几个原件再怎么放置都会发生碰撞冲突
    else geometry_layer_output(bestSolution);//绘制输出函数
}

void Beamsearch::test() {//测试函数，现在测试的就是PolygonModification2()，这个函数有很大的问题，我们项目就进展到这了
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
        std::cout << "文件打开失败" << endl;
        return;
    }
    string line;
    while (getline(infile, line)) {//每次从文件读取一行
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
        if (PL().HGP_2D_Polygon_Is_Clockwise_Oriented_C(points))//防止点的顺序颠倒而导致生成的多边形是负的（多边形边你可以理解为是有向的，所以我们得考虑顺逆时针问题）
        {
            std::reverse(points.begin(), points.end());
        }
        polygons.push_back(points);
    }
    origin_polygons = polygons;//获得现在待排列的最开始的多边形们
}

std::vector<Vector2d1> Beamsearch::perior_geometry_put()//处理出现重复的元件
{
    std::vector<Vector2d1> ans = polygons;
    std::cout << "图形种类加入是否重复" << endl;
    bool ques;
    cin >> ques;
    if (ques) {
        std::cout << "重复的有几个" << endl;
        int n; cin >> n;
        while (n--) {
            int type;
            cin >> type;
            if (type > polygons.size()) {
                std::cout << "没有该类型的几何结构哦，请重新输入" << endl;
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
    // 计算图像的尺寸

    // 创建一个黑色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_64FC3，初始值为黑色
    cv::Mat rightimage(boxy, boxx / 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // 绘制多边形
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // 绘制多边形线条
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }

    // 左右翻转图像
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // 绘制一条垂直线
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // 显示图像
    cv::imshow("Polygons", symmetric_image);
    cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_output2(vector<Vector2d1> a, vector<Vector2d1> b) {
    // 计算图像的尺寸
    static int count = 0;

    string path = this->image_path;
    path = path + "output" + to_string(count) + ".jpg";
    cout << path << endl;
    ++count;

    // 创建一个黑色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_64FC3，初始值为黑色
    cv::Mat rightimage(boxy, boxx / 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // 绘制多边形
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // 绘制多边形线条
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }

    // 绘制多边形
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // 绘制多边形线条
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 0, 0), 2);
    }

    // 左右翻转图像
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // 绘制一条垂直线
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // 显示图像
    //cv::imshow("Polygons", symmetric_image);
    cv::imwrite(path, symmetric_image);
    //cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_save(vector<Vector2d1> a, int num, double score) {
    // 计算图像的尺寸
    string path = this->image_path;
    path = path + "/节点" + to_string(num) + "评分：" + to_string(score) + ".jpg";
    cout << path << endl;

    // 创建一个黑色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_64FC3，初始值为黑色
    cv::Mat rightimage(boxy, boxx / 2, CV_64FC3, cv::Scalar(0, 0, 0));

    // 绘制多边形
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        const cv::Point* pts = points.data();
        int num_points = points.size();
        // 绘制多边形线条
        cv::polylines(rightimage, &pts, &num_points, 1, true, cv::Scalar(255, 255, 255), 2);
    }

    // 左右翻转图像
    cv::Mat leftimage;
    cv::flip(rightimage, leftimage, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image;
    cv::hconcat(leftimage, rightimage, symmetric_image);

    // 绘制一条垂直线
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image, point1, point2, cv::Scalar(0, 0, 255), 1);

    // 保存图像
    cv::imwrite(path, symmetric_image);
    cv::waitKey(0);
    return;
}

void Beamsearch::geometry_layer_save1(vector<Vector2d1> a, vector<Vector2d1> b) {
    // 计算图像的尺寸
    string path = this->image_path;
    string path1 = path + "/Roughing.jpg";
    string path2 = path + "/Finishing.jpg";
    string path3 = path + "/BothOfThem.jpg";

    // 创建一个白色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_8UC3，初始值为白色
    cv::Mat rightimage1(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));
    std::vector<std::vector<cv::Point>> pts;

    // 绘制多边形，使用黑色填充
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts.push_back(points);
    }
    // 填充多边形
    cv::fillPoly(rightimage1, pts, cv::Scalar(0, 0, 0));

    // 左右翻转图像
    cv::Mat leftimage1;
    cv::flip(rightimage1, leftimage1, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image1;
    cv::hconcat(leftimage1, rightimage1, symmetric_image1);

    // 绘制一条垂直线
    cv::Point point1(boxx / 2, boxy);
    cv::Point point2(boxx / 2, 0);
    cv::line(symmetric_image1, point1, point2, cv::Scalar(0, 0, 255), 1);

    // 保存图像
    cv::imwrite(path1, symmetric_image1);

    // 创建一个白色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_8UC3，初始值为白色
    cv::Mat rightimage2(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));
    std::vector<std::vector<cv::Point>> pts1;

    // 绘制多边形，使用蓝色填充
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts1.push_back(points);
    }
    // 填充多边形
    cv::fillPoly(rightimage2, pts1, cv::Scalar(255, 0, 0));

    // 左右翻转图像
    cv::Mat leftimage2;
    cv::flip(rightimage2, leftimage2, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image2;
    cv::hconcat(leftimage2, rightimage2, symmetric_image2);

    // 绘制一条垂直线
    cv::Point point3(boxx / 2, boxy);
    cv::Point point4(boxx / 2, 0);
    cv::line(symmetric_image2, point3, point4, cv::Scalar(0, 0, 255), 1);

    // 保存图像
    cv::imwrite(path2, symmetric_image2);

    // 创建一个白色的图像，尺寸为(boxy, boxx / 2)，数据类型为CV_8UC3，初始值为白色
    cv::Mat rightimage3(boxy, boxx / 2, CV_8UC3, cv::Scalar(255, 255, 255));

    std::vector<std::vector<cv::Point>> pts2;

    // 绘制多边形，使用黑色填充
    for (const auto& polygon : a) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts2.push_back(points);
    }
    // 填充多边形，使用黑色填充
    cv::fillPoly(rightimage3, pts2, cv::Scalar(0, 0, 0));

    std::vector<std::vector<cv::Point>> pts3;

    // 绘制多边形，使用蓝色填充
    for (const auto& polygon : b) {
        std::vector<cv::Point> points;
        for (const auto& vertex : polygon) {
            // 将顶点坐标转换为OpenCV图像坐标系中的坐标
            int x = vertex.x;
            int y = boxy - vertex.y; // 在OpenCV中，图像的原点位于左上角，所以需要翻转y轴
            cv::Point point(x, y);
            points.push_back(point);
        }
        pts3.push_back(points);
    }
    // 填充多边形，使用蓝色填充
    cv::fillPoly(rightimage3, pts3, cv::Scalar(0, 0, 255));

    // 左右翻转图像
    cv::Mat leftimage3;
    cv::flip(rightimage3, leftimage3, 1);

    // 拼接左右图像，得到对称图像
    cv::Mat symmetric_image3;
    cv::hconcat(leftimage3, rightimage3, symmetric_image3);

    // 绘制一条垂直线
    cv::Point point5(boxx / 2, boxy);
    cv::Point point6(boxx / 2, 0);
    cv::line(symmetric_image3, point5, point6, cv::Scalar(0, 0, 255), 1);

    // 保存图像
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

void Beamsearch::PolygonModification2()
{
}
