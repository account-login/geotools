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

    static void calc_ab(const GeoLine &line, double &a, double &b) {
        // a * x1 + b * y1 = z1
        // a * x2 + b * y2 = z2

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

    static float cross_with_lon(const GeoLine &line, float lon) {
        double a, b;
        calc_ab(line, a, b);

        // tan(B) = a * cos(A) + b * sin(A)
        double A = deg2rad(lon);
        double B = atan(a * cos(A) + b * sin(A));

        return (float)rad2deg(B);   // lat
    }

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

    static void cut_ns_ex(const GeoLine &line, float NS, uint32_t flags[3], GeoLine lines[3])
    {
        flags[0] = flags[1] = flags[2] = D_NONE;

        if (line.src.lat == line.dst.lat) {
            // special case: horizontal line
            flags[0] = line.src.lat >= NS ? D_N : D_S;
            lines[0] = line;
            return;
        }

        GeoLonLat np = line.src.lat < line.dst.lat ? line.dst : line.src;
        GeoLonLat sp = line.src.lat < line.dst.lat ? line.src : line.dst;
        assert(np.lat != sp.lat);
        GeoLonLat wp = line.src.lon < line.dst.lon ? line.src : line.dst;
        GeoLonLat ep = line.src.lon < line.dst.lon ? line.dst : line.src;

        double a, b;
        calc_ab(line, a, b);

        float A_1_deg = NAN, A_2_deg = NAN;
        bool A_1_ok = false, A_2_ok = false;

        // tan(B) = a * cos(A) + b * sin(A)
        //
        //            b ± sqrt(a*a + b*b - tan(B)*tan(B))     b ± c1
        // tan(A/2) = ----------------------------------- = ----------
        //                        a + tan(B)                a + tan(B)
        double tan_B = tan(deg2rad(NS));
        double c2 = a*a + b*b - tan_B*tan_B;
        if (c2 > 0) {
            // A have 2 solutions
            // NOTE: the 1 solutions case is discarded
            double c1 = sqrt(c2);
            double A_1 = 2 * atan((b + c1) / (a + tan_B));
            double A_2 = 2 * atan((b - c1) / (a + tan_B));
            double A_1_deg_d = rad2deg(A_1);
            double A_2_deg_d = rad2deg(A_2);
            // is solutions on the arc?
            A_1_ok = wp.lon < A_1_deg_d && A_1_deg_d < ep.lon;
            A_2_ok = wp.lon < A_2_deg_d && A_2_deg_d < ep.lon;
            A_1_deg = (float)A_1_deg_d;
            A_2_deg = (float)A_2_deg_d;
            // NOTE: it is possible that A_1_deg == A_2_deg due to limited precision
            if (A_1_deg == A_2_deg) {
                A_1_ok = A_2_ok = false;
            }
        }

        GeoLonLat cross_1 = GeoLonLat(A_1_deg, NS);
        GeoLonLat cross_2 = GeoLonLat(A_2_deg, NS);

        if (A_1_ok && A_2_ok) {
            // line is splited to 3 segments by NS
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
            // line is not splited or tangent to NS
            // NOTE: using mid_ns instead of sp.lat because sp may be slightly off the box
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
        // 180 * 2^(-23) ≈ 2.1e-05
        // 3e-5 / 360 * 40000 km ≈ 3.3 m
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
        return line.src.lon == line.dst.lon;
    }

    static EdgeNode *insert_rec_vertical(const EdgeInsertCtx &ctx, EdgeNode *node);

    static EdgeNode *insert_rec(const EdgeInsertCtx &ctx, EdgeNode *node) {
        if (is_line_vertical(ctx.line)) {
            // special case: vertical line
            // NOTE: vertical line may be generated from cut_ns_ex due to float pointer error
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
