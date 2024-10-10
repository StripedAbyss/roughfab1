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
    // 直接生成偏移后的多边形
	for (const auto& poly : pys) {
    Polygon_2 offset_poly = poly.outline().offset(offset_value); // 设定的偏移值
    if (!offset_poly.is_empty()) {
        ans.push_back(offset_poly);
    }
	}
	// 处理生成的多边形
    for (const auto& it : ans) {
        output.push_back(Convert_Polygon_2_to_Vector2d1(it));
        if (abs(it.area()) < 500) {
            continue; // 跳过较小的晶胞
        }

        // 计算包围盒面积和周长
        areas += it.bbox().x_span() * it.bbox().y_span();
        double length = 0;
        
        // 计算边长总和
        for (auto itt = it.edges_begin(); itt != it.edges_end(); itt++) {
            length += sqrt(itt->squared_length());
        }

        // 计算得分
        score += it.bbox().x_span() * it.bbox().y_span() * (it.bbox().x_span() + it.bbox().y_span()) * 2 / length;
        now++;
    }

    if (areas == 0) return make_pair(0, now);
    score /= areas;
    score = score * 0.9 + 0.1 * min(previous / now, 1); // 调整得分公式

    return make_pair(score, now); // 返回得分与晶胞数量
}
 /*   vector<Point_2> getit;
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
}*/ 
