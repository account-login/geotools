#include "edgetree.h"

#include <cmath>
#include <cassert>


namespace geotools {

    struct EdgeInsertCtx {
        GeoLine     line;
        GeoBox      box;
        uint32_t    depth;
        uint32_t    max_depth;

        EdgeInsertCtx(const GeoLine &line, const GeoBox &box, uint32_t depth, uint32_t max_depth)
            : line(line), box(box), depth(depth), max_depth(max_depth)
        {}
    };

    static EdgeNode *make_node(EdgeNode *node) {
        if (!node) {
            node = new EdgeNode(EDGENODE_INNER);
        }
        return node;
    }

    // 计算 line 跟原点组成的平面，line 不能跟经线平行
    // a * x1 + b * y1 = z1
    // a * x2 + b * y2 = z2

    // 直角坐标转经纬度
    // A 是经度，B 是纬度
    // Z = sin(B),
    // X = cos(B) * cos(A),
    // Y = cos(B) * sin(A),
    static void calc_ab(const GeoLine &line, double &a, double &b) {
        const GeoLonLat &p1 = line.src;
        const GeoLonLat &p2 = line.dst;
        assert(p1.lon != p2.lon);

        double a1 = deg2rad(p1.lon);
        double a2 = deg2rad(p2.lon);
        double sin_a1 = sin(a1);
        double sin_a2 = sin(a2);
        double cos_a1 = cos(a1);
        double cos_a2 = cos(a2);
        double tan_b1 = tan(deg2rad(p1.lat));
        double tan_b2 = tan(deg2rad(p2.lat));

        // sin(a1 - a2)
        double sin_a1_subs_a2 = sin_a1 * cos_a2 - sin_a2 * cos_a1;
        //     sin(a1) * tan(b2) - sin(a2) * tan(b1)
        // a = -------------------------------------
        //                sin(a1 - a2)
        a = (sin_a1 * tan_b2 - sin_a2 * tan_b1) / sin_a1_subs_a2;
        //     cos(a2) * tan(b1) - cos(a1) * tan(b2)
        // b = -------------------------------------
        //                sin(a1 - a2)
        b = (cos_a2 * tan_b1 - cos_a1 * tan_b2) / sin_a1_subs_a2;
    }

    // 经线跟 line 所在圆相交的纬度
    static float cross_with_lon(const GeoLine &line, float lon) {
        double a, b;
        calc_ab(line, a, b);

        // 由 a * x + b * y = z 平面
        // 得到 tan(B) = a * cos(A) + b * sin(A)
        double A = deg2rad(lon);
        double B = atan(a * cos(A) + b * sin(A));

        return (float)rad2deg(B);   // lat
    }

    // 将圆弧 line 用经线 WE 切割成东西两段
    static uint32_t cut_we(const GeoLine &line, float WE, GeoLine &w, GeoLine &e) {
        GeoLonLat wp = line.src.lon < line.dst.lon ? line.src : line.dst;
        GeoLonLat ep = line.src.lon < line.dst.lon ? line.dst : line.src;

        if (ep.lon <= WE) {
            w = line;
            return D_W;
        }
        if (WE <= wp.lon) {
            e = line;
            return D_E;
        }
        assert(wp.lon != ep.lon);

        GeoLonLat cross(WE, cross_with_lon(line, WE));
        w.src = wp;
        w.dst = cross;
        e.src = cross;
        e.dst = ep;
        return D_W | D_E;
    }

    // wp, ep 分别是圆弧西边和东边的端点
    // 判断圆弧于纬线的交点 (lon, lat) 是否在圆弧内部（不包括圆弧两端）
    // 比较范围时 (lon, lat) 是 double，float 经度不够
    static bool is_cross_on_arc(double lon, double lat, GeoLonLat wp, GeoLonLat ep) {
        GeoLonLat cross((float)lon, (float)lat);
        return wp.lon <= lon && lon <= ep.lon   // 这里用的闭区间
                                                // 因为斜率很大时，lon 与圆弧两端重合，但 lat 不重合
            && cross != wp && cross != ep;      // 再剔除与圆弧两端重合的交点
    }

    // line 是一条不平行于经线的圆弧，可能被纬线 NS 分割。
    // line 所在的圆跟纬线 NS 有 3 种关系：不相交，有两个交点，相切
    // line 这段圆弧被纬线 NS 分割有 5 种情况：
    //      1. 不相交，line 要么在 NS 南边，要么在北边
    //      1. 被割成南北 2 段
    //      2. 被割成 3 段，如果 line 在北半球，两头在南边，中间在北边，如果 line 在南半球，方向相反
    //      3. 相切，当成不相交处理
    //      4. 重合，NS 是赤道的情况，也当成不相交处理
    static void cut_ns_ex(const GeoLine &line, float NS, uint32_t flags[3], GeoLine lines[3])
    {
        flags[0] = flags[1] = flags[2] = D_NONE;

        GeoLonLat np = line.src.lat < line.dst.lat ? line.dst : line.src;
        GeoLonLat sp = line.src.lat < line.dst.lat ? line.src : line.dst;
        GeoLonLat wp = line.src.lon < line.dst.lon ? line.src : line.dst;
        GeoLonLat ep = line.src.lon < line.dst.lon ? line.dst : line.src;

        double a, b;
        calc_ab(line, a, b);

        // 两个交点的经度
        float A_1_deg = NAN, A_2_deg = NAN;
        // 表示 A_1, A_2 是否在 line 上
        bool A_1_ok = false, A_2_ok = false;

        // line 所在平面
        // tan(B) = a * cos(A) + b * sin(A)
        //
        // 通过纬度 B 计算经度 A
        //            b ± sqrt(a*a + b*b - tan(B)*tan(B))     b ± c1
        // tan(A/2) = ----------------------------------- = ----------
        //                        a + tan(B)                a + tan(B)
        double tan_B = tan(deg2rad(NS));
        double c2 = a*a + b*b - tan_B*tan_B;
        // c2 > 0 说明 line 所在的圆跟 NS 有两个交点
        // c2 == 0 说明 line 所在的圆跟 NS 相切，当成不相交处理
        // c2 == 0 也包含了 line 所在的圆跟 NS 重合，NS 是赤道的情况，也当成不相交处理
        if (c2 > 0) {
            // 只有 c2 >= 0 且 line 跟 NS 不重合才能计算 A_1, A_2
            // 否则会得到 NAN
            double c1 = sqrt(c2);
            double A_1 = 2 * atan((b + c1) / (a + tan_B));
            double A_2 = 2 * atan((b - c1) / (a + tan_B));
            assert(!isnan(A_1));
            assert(!isnan(A_2));

            double A_1_deg_d = rad2deg(A_1);
            double A_2_deg_d = rad2deg(A_2);
            // 这里比较范围用 double，float 精度不够
            // line 斜率比较大时，A_1, A_2 的误差比较敏感，有可能 A_1, A_2 跟端点的经度一样了
            // FIXME: handle large slope
            A_1_ok = is_cross_on_arc(A_1_deg_d, NS, wp, ep);
            A_2_ok = is_cross_on_arc(A_2_deg_d, NS, wp, ep);

            A_1_deg = (float)A_1_deg_d;
            A_2_deg = (float)A_2_deg_d;
            // 可能精度不足，导致 A_1_deg == A_2_deg，当成相切处理
            if (A_1_deg == A_2_deg) {
                A_1_ok = A_2_ok = false;
            }
        }

        GeoLonLat cross_1 = GeoLonLat(A_1_deg, NS);
        GeoLonLat cross_2 = GeoLonLat(A_2_deg, NS);

        if (A_1_ok && A_2_ok) {
            // 3 segments
            assert(NS != 0);

            GeoLonLat wap = A_1_deg < A_2_deg ? cross_1 : cross_2;
            GeoLonLat eap = A_1_deg < A_2_deg ? cross_2 : cross_1;

            if (NS > 0) {
                flags[1] = D_N;
                flags[0] = flags[2] = D_S;
            } else {
                flags[1] = D_S;
                flags[0] = flags[2] = D_N;
            }

            lines[0].src = wp;
            lines[0].dst = wap;
            lines[1].src = wap;
            lines[1].dst = eap;
            lines[2].src = eap;
            lines[2].dst = ep;
        } else if (!A_1_ok && !A_2_ok) {
            // 全部在南北，或者全部在北边
            // 由于计算有误差，可能这个圆弧在大部分北边，南边稍微有一段，但是算出交点却不足圆弧内部
            // 导致应该被分割的圆弧没有被分割，所以判断它在哪边时，不能用圆弧两端的纬度
            float mid_ns = (np.lat + sp.lat) / 2;
            flags[0] = NS <= mid_ns ? D_N : D_S;
            lines[0] = line;
        } else {
            // 2 segments
            GeoLonLat cross = A_1_ok ? cross_1 : cross_2;
            flags[0] = D_N;
            flags[1] = D_S;

            lines[0].src = np;
            lines[0].dst = cross;
            lines[1].src = cross;
            lines[1].dst = sp;
        }
    }

    // 没用
    static uint32_t cut_ns(const GeoLine &line, float NS, GeoLine &n, GeoLine &s) {
        GeoLonLat np = line.src.lat < line.dst.lat ? line.dst : line.src;
        GeoLonLat sp = line.src.lat < line.dst.lat ? line.src : line.dst;

        if (NS <= sp.lat) {
            n = line;
            return D_N;
        }
        if (np.lat <= NS) {
            s = line;
            return D_S;
        }
        assert(np.lat != sp.lat);

        double a, b;
        calc_ab(line, a, b);

        // tan(B) = a * cos(A) + b * sin(A)
        //
        //            b ± sqrt(a*a + b*b - tan(B)*tan(B))     b ± c1
        // tan(A/2) = ----------------------------------- = ----------
        //                        a + tan(B)                a + tan(B)
        double tan_B = tan(deg2rad(NS));
        double c2 = a*a + b*b - tan_B*tan_B;
        double c1 = sqrt(c2);
        double A_1 = 2 * atan((b + c1) / (a + tan_B));
        double A_2 = 2 * atan((b - c1) / (a + tan_B));

        GeoLonLat cross((float)rad2deg(A_1), NS);
        GeoLonLat wp = line.src.lon < line.dst.lon ? line.src : line.dst;
        GeoLonLat ep = line.src.lon < line.dst.lon ? line.dst : line.src;
        if (cross.lon < wp.lon || cross.lon > ep.lon) {
            cross.lon = (float)rad2deg(A_2);
        }
        assert(wp.lon <= cross.lon && cross.lon <= ep.lon);

        n.src = np;
        n.dst = cross;
        s.src = cross;
        s.dst = sp;
        return D_N | D_S;
    }

    static bool box_contains_err(const GeoBox &box, GeoLonLat lonlat, float e) {
        return box.W - e <= lonlat.lon && lonlat.lon <= box.E + e
            && box.S - e <= lonlat.lat && lonlat.lat <= box.N + e;
    }

    static void assert_box_contains(const GeoBox &box, const GeoLine &line) {
        // float 里有 23 个 bit 表示小数
        // float 表示经纬度误差度数：180 * 2^(-23) ≈ 2.1e-05
        // 取一个稍微大一点的数 3e-5 作为比较经纬度的允许误差
        // 对应的距离 3e-5 / 360 * 40000 km ≈ 3.3 m
        const float error = 3e-5;
#ifdef TWOBLUECUBES_SINGLE_INCLUDE_CATCH_HPP_INCLUDED
        if (!box_contains_err(box, line.src, error) || !box_contains_err(box, line.dst, error)) {
            CAPTURE(box);
            CAPTURE(line);
            CHECK(false);
        }
#endif
        assert(box_contains_err(box, line.src, error));
        assert(box_contains_err(box, line.dst, error));
    }

    static void insert_rec_sub(
        EdgeNode *func(const EdgeInsertCtx &ctx, EdgeNode *node),
        const EdgeInsertCtx &ctx, EdgeNode *node, const GeoLine &line, uint32_t flag)
    {
        node->child(flag) = func(
            EdgeInsertCtx(line, ctx.box.get(flag), ctx.depth + 1, ctx.max_depth),
            node->child(flag)
        );
    }

    static EdgeNode *insert_rec(const EdgeInsertCtx &ctx, EdgeNode *node);

    static void insert_rec_ns(
        const EdgeInsertCtx &ctx, EdgeNode *node, const GeoLine &wel, uint32_t flag_we)
    {
        assert(flag_we == D_W || flag_we == D_E);
        assert_box_contains(ctx.box, wel);

        float NS = (ctx.box.N + ctx.box.S) / 2.0f;
        uint32_t flags[3];
        GeoLine lines[3];
        cut_ns_ex(wel, NS, flags, lines);
        for (size_t i = 0; i < 3; ++i) {
            if (flags[i] != D_NONE) {
                assert_box_contains(ctx.box, lines[i]);
                insert_rec_sub(insert_rec, ctx, node, lines[i], flags[i] | flag_we);
            }
        }
    }

    static bool is_line_vertical(const GeoLine &line) {
        // TODO: consider slope
        return line.src.lon == line.dst.lon;
    }

    static EdgeNode *insert_rec_vertical(const EdgeInsertCtx &ctx, EdgeNode *node);

    static EdgeNode *insert_rec(const EdgeInsertCtx &ctx, EdgeNode *node) {
        if (is_line_vertical(ctx.line)) {
            // special case
            // 上面很多计算都假设圆弧不平行经线，所以这个要特殊处理
            // cut_ns_ex() 函数由于计算误差，可能把一根斜率比较大的圆弧变的平行于经线
            // 所以这个 case 会随时出现
            return insert_rec_vertical(ctx, node);
        }

        assert_box_contains(ctx.box, ctx.line);

        if (ctx.depth >= ctx.max_depth) {
            // leaf
            if (!node) {
                node = new EdgeNode(EDGENODE_LEAF);
            }
            node->lines.push_back(ctx.line);
        } else {
            float WE = (ctx.box.W + ctx.box.E) / 2.0f;
            GeoLine wl, el;
            uint32_t flag_we = cut_we(ctx.line, WE, wl, el);

            node = make_node(node);
            assert(node->type == EDGENODE_INNER);
            if (flag_we & D_W) {
                insert_rec_ns(ctx, node, wl, D_W);
            }
            if (flag_we & D_E) {
                insert_rec_ns(ctx, node, el, D_E);
            }
        }
        return node;
    }

    static bool is_line_cross_180(const GeoLine &line) {
        return abs(line.src.lon - line.dst.lon) > 180.0;
    }

    static uint32_t split_line_cross_180(const GeoLine &line, GeoLine &w, GeoLine &e) {
        GeoLonLat wp = line.src.lon < line.dst.lon ? line.src : line.dst;
        GeoLonLat ep = line.src.lon < line.dst.lon ? line.dst : line.src;

        float lat = cross_with_lon(line, 180.0);

        w.src = GeoLonLat(-180.0f, lat);
        w.dst = wp;
        e.src = ep;
        e.dst = GeoLonLat(+180.0f, lat);

        uint32_t flag = D_NONE;
        if (w.src != w.dst) {
            flag |= D_W;
        }
        if (e.src != e.dst) {
            flag |= D_E;
        }
        return flag;
    }

    static EdgeNode *insert_rec_vertical(const EdgeInsertCtx &ctx, EdgeNode *node) {
        if (ctx.depth >= ctx.max_depth) {
            // leaf
            if (!node) {
                node = new EdgeNode(EDGENODE_LEAF);
            }
            node->lines.push_back(ctx.line);
        } else {
            float WE = (ctx.box.W + ctx.box.E) / 2.0f;
            float NS = (ctx.box.N + ctx.box.S) / 2.0f;

            const GeoLine &line = ctx.line;
            GeoLonLat np = line.src.lat < line.dst.lat ? line.dst : line.src;
            GeoLonLat sp = line.src.lat < line.dst.lat ? line.src : line.dst;
            assert(ctx.box.contains(np));
            assert(ctx.box.contains(sp));

            // cut ns
            GeoLine nl, sl;
            uint32_t flag_ns;
            if (sp.lat >= NS) {
                flag_ns = D_N;
                nl = line;
            } else if (np.lat <= NS) {
                flag_ns = D_S;
                sl = line;
            } else {
                flag_ns = D_N | D_S;
                GeoLonLat mid(line.src.lon, NS);
                nl.src = np;
                nl.dst = mid;
                sl.src = mid;
                sl.dst = sp;
            }

            // insert ns
            node = make_node(node);
            assert(node->type == EDGENODE_INNER);
            uint32_t flag_we = line.src.lon < WE ? D_W : D_E;
            if (flag_ns & D_N) {
                insert_rec_sub(insert_rec_vertical, ctx, node, nl, flag_we | D_N);
            }
            if (flag_ns & D_S) {
                insert_rec_sub(insert_rec_vertical, ctx, node, sl, flag_we | D_S);
            }
        }
        return node;
    }

    void EdgeTree::insert(const GeoLine &line) {
        assert(line.src.is_valid());
        assert(line.dst.is_valid());

        if (is_line_cross_180(line)) {
            // 如果圆弧 line 两端的经度相差超过 180 度，那么 line 应该穿过 ±180 经线
            // 由于 ±180 经线是边界，line 只能拆成两段
            GeoLine w, e;
            uint32_t flag = split_line_cross_180(line, w, e);
            if (flag & D_W) {
                this->root = insert_rec(
                    EdgeInsertCtx(w, GeoBox(), 0, this->precision),
                    this->root
                );
            }
            if (flag & D_E) {
                this->root = insert_rec(
                    EdgeInsertCtx(e, GeoBox(), 0, this->precision),
                    this->root
                );
            }
        } else {
            this->root = insert_rec(
                EdgeInsertCtx(line, GeoBox(), 0, this->precision),
                this->root
            );
        }
    }

}   // ::geotools
