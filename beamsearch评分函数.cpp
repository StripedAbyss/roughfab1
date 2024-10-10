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
    // ֱ������ƫ�ƺ�Ķ����
	for (const auto& poly : pys) {
    Polygon_2 offset_poly = poly.outline().offset(offset_value); // �趨��ƫ��ֵ
    if (!offset_poly.is_empty()) {
        ans.push_back(offset_poly);
    }
	}
	// �������ɵĶ����
    for (const auto& it : ans) {
        output.push_back(Convert_Polygon_2_to_Vector2d1(it));
        if (abs(it.area()) < 500) {
            continue; // ������С�ľ���
        }

        // �����Χ��������ܳ�
        areas += it.bbox().x_span() * it.bbox().y_span();
        double length = 0;
        
        // ����߳��ܺ�
        for (auto itt = it.edges_begin(); itt != it.edges_end(); itt++) {
            length += sqrt(itt->squared_length());
        }

        // ����÷�
        score += it.bbox().x_span() * it.bbox().y_span() * (it.bbox().x_span() + it.bbox().y_span()) * 2 / length;
        now++;
    }

    if (areas == 0) return make_pair(0, now);
    score /= areas;
    score = score * 0.9 + 0.1 * min(previous / now, 1); // �����÷ֹ�ʽ

    return make_pair(score, now); // ���ص÷��뾧������
}
 /*   vector<Point_2> getit;
    CGAL_2D_Polygon_Dart_Sampling_a(pys, 0.002, getit);//��ɢȡ�㣬�ж��Ƿ���ͼ����࣬���ص㼯
    //CGAL_2D_Polygon_Dart_Sampling_b(pys, 0.5, getit, 100);//��ɢȡ�㣬�ж��Ƿ���ͼ����࣬���ص㼯
    ans = get_triangulation_net(getit, pys);//���ݵ㼯���ɾ���
    for (auto it = ans.begin(); it != ans.end(); it++) {
        output.push_back(Convert_Polygon_2_to_Vector2d1(*it));
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
}*/ 
