#include <stdlib.h>

#include <boost/assign/list_of.hpp>

#include "catch.h"
#include "string_fmt.hpp"
#include "../edgetree.cpp"


using namespace std;
using namespace geotools;

using boost::assign::list_of;
using boost::assign::map_list_of;


TEST_CASE("is_line_cross_180") {
    CHECK(!is_line_cross_180(GeoLine(0, 0, 0, 0)));
    CHECK(!is_line_cross_180(GeoLine(0, 0, 180, 0)));
    CHECK(is_line_cross_180(GeoLine(-1, 0, 180, 0)));
}

TEST_CASE("split_line_cross_180.partial") {
    GeoLine w, e;
    uint32_t flag = split_line_cross_180(GeoLine(-1, 0, 180, 0), w, e);
    REQUIRE(flag == D_W);
    CHECK(w == GeoLine(-180, 0, -1, 0));
    CHECK(e == GeoLine(180, 0, 180, 0));
}

static double line_angle(const GeoLine &line) {
    return geo_angle(line.src.lon, line.src.lat, line.dst.lon, line.dst.lat);
}

static void test_same_mirror(const GeoLine &l1, const GeoLine &l2) {
    SECTION("same mirror") {
        double a1, b1, a2, b2;
        calc_ab(l1, a1, b1);
        calc_ab(l2, a2, b2);
        CHECK(Approx(a1) == a2);
        CHECK(Approx(b1) == b2);
    }
}

static void test_same_distance(const GeoLine &sum, const GeoLine &p1, const GeoLine &p2) {
    SECTION("same distance") {
        CHECK(line_angle(p1) + line_angle(p2) == Approx(line_angle(sum)));
    }
}

static void test_split_line_cross_180(const GeoLine &line) {
    REQUIRE(is_line_cross_180(line));

    GeoLine w, e;
    uint32_t flag = split_line_cross_180(line, w, e);
    REQUIRE(flag == (D_W | D_E));

    SECTION("connected with line") {
        // NOTE: assuming line.src is west
        CHECK(w.dst == line.src);
        CHECK(line.dst == e.src);
        CHECK(w.src.lon == -180);
        CHECK(e.dst.lon == +180);
        CHECK(w.src.lat == e.dst.lat);
    }

    test_same_mirror(w, e);
    test_same_distance(line, w, e);
}

TEST_CASE("split_line_cross_180.full") {
    GeoLine line = GeoLine(-2, 20, 179, -70);
    test_split_line_cross_180(line);
}

static float rand_range(float lo, float hi) {
    float r = (float)::random();
    return r / RAND_MAX * (hi - lo) + lo;
}

static unsigned int rand_seed() {
    timespec ts;
    ::clock_gettime(CLOCK_REALTIME, &ts);
    unsigned int seed = (unsigned int)ts.tv_nsec;
    ::srandom(seed);
    return seed;
}

static unsigned int g_seed = rand_seed();

static GeoLine rand_line() {
    GeoLine line(
        rand_range(-179.999f, +180), rand_range(-80, +80),
        rand_range(-179.999f, +180), rand_range(-80, +80)
    );
    if (line.dst.lon < line.src.lon) {
        swap(line.src, line.dst);
    }
    return line;
}

TEST_CASE("split_line_cross_180.random") {
    timespec ts;
    ::clock_gettime(CLOCK_REALTIME, &ts);
    ::srandom((unsigned int)ts.tv_nsec);

    const size_t N = 1024;
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        GeoLine line = rand_line();
        if (line.src.lon == line.dst.lon) {
            continue;
        }

        if (is_line_cross_180(line)) {
            test_split_line_cross_180(line);
            count++;
        }
    }
    CHECK(count > 0);
}

namespace tz {
    template <>
    string str(const GeoLine &line) {
        return strfmt("GeoLine(%f, %f, %f, %f)",
            line.src.lon, line.src.lat, line.dst.lon, line.dst.lat);
    }
}

namespace Catch {
    template <>
    struct StringMaker<GeoBox> {
        static std::string convert(const GeoBox &box) {
            return tz::strfmt("GeoBox(%f, %f, %f, %f)", box.W, box.E, box.N, box.S);
        }
    };

    template <>
    struct StringMaker<GeoLine> {
        static std::string convert(const GeoLine &line) {
            return tz::str(line);
        }
    };

    typedef vector<GeoLine> GeoLineList;
    template <>
    struct StringMaker<GeoLineList> {
        static std::string convert(const GeoLineList &lines) {
            return tz::repr_set(lines, 1024);
        }
    };
}

void test_cut(
    uint32_t cut_func(const GeoLine &line, float NS, GeoLine &n, GeoLine &s),
    const GeoLine &line, float pos, uint32_t flag_expected)
{
    GeoLine l1, l2;
    uint32_t flag = cut_func(line, pos, l1, l2);
    REQUIRE(flag == flag_expected);

    // NOTE: assuming line.src is west
    CHECK(line.src == l1.src);
    CHECK(line.dst == l2.dst);
    CHECK(l1.dst == l2.src);

    test_same_mirror(l1, l2);
    test_same_distance(line, l1, l2);
}

TEST_CASE("cut_we.rand") {
    const size_t N = 1024;
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        GeoLine line = rand_line();
        if (line.src.lon == line.dst.lon) {
            continue;
        }
        if (is_line_cross_180(line)) {
            continue;
        }

        float WE = rand_range(line.src.lon, line.dst.lon);
        if (WE == line.src.lon || WE == line.dst.lon) {
            continue;
        }
        test_cut(cut_we, line, WE, D_W | D_E);
        count++;
    }
    CHECK(count > 0);
}

TEST_CASE("cut_ns.rand") {
    const size_t N = 1024;
    size_t count = 0;
    for (size_t i = 0; i < N; ++i) {
        GeoLine line = rand_line();
        if (line.src.lon == line.dst.lon) {
            continue;
        }
        if (is_line_cross_180(line)) {
            continue;
        }

        if (line.src.lat < line.dst.lat) {
            swap(line.src, line.dst);   // src is north
        }
        float NS = rand_range(line.src.lat, line.dst.lat);
        if (NS == line.src.lat || NS == line.dst.lat) {
            continue;
        }
        test_cut(cut_ns, line, NS, D_N | D_S);
        count++;
    }
    CHECK(count > 0);
}

static size_t check_edgenode(const EdgeNode *node, const GeoBox &box) {
    if (!node) {
        return 0;
    }

    if (node->type == EDGENODE_INNER) {
        return check_edgenode(node->NW, box.get(D_NW)) + check_edgenode(node->NE, box.get(D_NE))
             + check_edgenode(node->SE, box.get(D_SE)) + check_edgenode(node->SW, box.get(D_SW));
    } else {
        // edgenode can not be empty
        CHECK(!node->lines.empty());
        // edgenode do not have child
        bool nochild = !node->NW && !node->NE && !node->SE && !node->SW;
        CHECK(nochild);
        // all endpoints within edgenode
        for (size_t i = 0; i < node->lines.size(); ++i) {
            const GeoLine &line = node->lines[i];
            assert_box_contains(box, line);
//            if (!box.contains(line.src) || !box.contains(line.dst)) {
//                CAPTURE(line);
//                CAPTURE(box);
//                assert(false);
//                REQUIRE(false);
//            }
//            REQUIRE(box.contains(line.src));
//            REQUIRE(box.contains(line.dst));
        }
        return node->lines.size();
    }
}

static void get_lines(const EdgeNode *node, vector<GeoLine> &lines) {
    if (!node) {
        return;
    }

    if (node->type == EDGENODE_INNER) {
        get_lines(node->NW, lines);
        get_lines(node->NE, lines);
        get_lines(node->SE, lines);
        get_lines(node->SW, lines);
    } else {
        lines.insert(lines.end(), node->lines.begin(), node->lines.end());
    }
}

static double get_line_sum(const EdgeNode *node) {
    if (!node) {
        return 0;
    }

    if (node->type == EDGENODE_INNER) {
        return get_line_sum(node->NW) + get_line_sum(node->NE)
             + get_line_sum(node->SE) + get_line_sum(node->SW);
    } else {
        double sum = 0;
        for (size_t i = 0; i < node->lines.size(); ++i) {
            const GeoLine &line = node->lines[i];
            sum += line_angle(line);
        }
        return sum;
    }
}

TEST_CASE("insert.zero-depth") {
    EdgeTree et1(0);
    et1.insert(GeoLine(-100, 80, 80, -80));
    CHECK(et1.root->type == EDGENODE_LEAF);
    CHECK(check_edgenode(et1.root, GeoBox()) == 1);

    // vertical line
    EdgeTree etv(0);
    etv.insert(GeoLine(20, -80, 20, 70));
    CHECK(etv.root->type == EDGENODE_LEAF);
    CHECK(check_edgenode(et1.root, GeoBox()) == 1);
}

TEST_CASE("insert.vertical") {
    EdgeTree etv(2);

    // add a line
    const GeoLine &line1 = GeoLine(20, -80, 20, 30);
    {
        etv.insert(line1);
        CHECK(check_edgenode(etv.root, GeoBox()) == 3);
        vector<GeoLine> expected = list_of
            (GeoLine(20, 30, 20, 0))
            (GeoLine(20, 0, 20, -45))
            (GeoLine(20, -45, 20, -80))
        ;
        vector<GeoLine> lines;
        get_lines(etv.root, lines);
        CHECK(lines == expected);
        CHECK(line_angle(line1) == Approx(get_line_sum(etv.root)));
    }

    // another line
    {
        const GeoLine &line2 = GeoLine(30, 10, 20, 60);
        etv.insert(line2);
        CHECK(check_edgenode(etv.root, GeoBox()) == 5);
        vector<GeoLine> lines;
        get_lines(etv.root, lines);
        vector<GeoLine> expected = list_of
            (GeoLine(20, 60, 24.729235f, 45))
            (GeoLine(20, 30, 20, 0))
            (GeoLine(24.729235f, 45, 30, 10))
            (GeoLine(20, 0, 20, -45))
            (GeoLine(20, -45, 20, -80))
        ;
        CHECK(lines == expected);
        CHECK(line_angle(line1) + line_angle(line2) == Approx(get_line_sum(etv.root)));
    }
}

TEST_CASE("insert.normal") {
    EdgeTree et(2);

    // add a line
    const GeoLine &line1 = GeoLine(20, -80, -100, 30);
    {
        et.insert(line1);
        CHECK(check_edgenode(et.root, GeoBox()) == 5);
        vector<GeoLine> expected = list_of
            (GeoLine(-100, 30, -94.692924f, 0))
            (GeoLine(0, -80.868248f, 20, -80))
            (GeoLine(-94.692924f, 0, -90, -27.053246f))
            (GeoLine(-90, -27.053246f, -85.474213f, -45))
            (GeoLine(-85.474213f, -45, 0, -80.868248f))
        ;
        vector<GeoLine> lines;
        get_lines(et.root, lines);
        CHECK(lines == expected);
        CHECK(line_angle(line1) == Approx(get_line_sum(et.root)));
    }

    // another line
    {
        const GeoLine &line2 = GeoLine(-100, 10, +90, 60);
        et.insert(line2);
        CHECK(check_edgenode(et.root, GeoBox()) == 8);
        vector<GeoLine> lines;
        get_lines(et.root, lines);
        vector<GeoLine> expected = list_of
            (GeoLine(-180, 84.728546f, -104.306900f, 45))
            (GeoLine(-100, 30, -94.692924f, 0))
            (GeoLine(-104.306900f, 45, -100, 10))
            (GeoLine(90, 60, 180, 84.728546f))
            (GeoLine(0, -80.868248f, 20, -80))
            (GeoLine(-94.692924f, 0, -90, -27.053246f))
            (GeoLine(-90, -27.053246f, -85.474213f,-45))
            (GeoLine(-85.474213f, -45, 0, -80.868248f))
        ;
        CHECK(lines == expected);
        CHECK(line_angle(line1) + line_angle(line2) == Approx(get_line_sum(et.root)));
    }
}

TEST_CASE("insert.rand") {
    EdgeTree et(10);

    const size_t N = 1024 * 2;
    size_t count = 0;
    double line_sum = 0;
    for (size_t i = 0; i < N; ++i) {
        GeoLine line = rand_line();
        if (line.src == line.dst) {
            continue;
        }

        line_sum += line_angle(line);
        et.insert(line);

        count++;
    }

    CHECK(count > 0);
    CHECK(check_edgenode(et.root, GeoBox()) > 0);
    CHECK(line_sum == Approx(get_line_sum(et.root)));
}

TEST_CASE("cut_ns_ex.error") {
    GeoLine line(0, 0, 0 + 1e-6f, 80);      // large slope
    uint32_t flags[3];
    GeoLine lines[3];
    cut_ns_ex(line, 1e-6, flags, lines);

    CHECK(flags[0] == D_N);
    CHECK(flags[1] == D_S);
    CHECK(flags[2] == D_NONE);

    CAPTURE(lines[0]);
    CAPTURE(lines[1]);
    CHECK(true);
}

TEST_CASE("cut_ns_ex.error-2") {
    GeoLine line(170, 0, 170 + 1e-5f, 80);      // large slope
    uint32_t flags[3];
    GeoLine lines[3];
    cut_ns_ex(line, 1e-6, flags, lines);

    CHECK(flags[0] == D_N);
    CHECK(flags[1] == D_S);
    CHECK(flags[2] == D_NONE);

    CAPTURE(lines[0]);
    CAPTURE(lines[1]);
    CHECK(true);
}

// TODO: test largest slope
// TODO: test cut_ns on small slope

//TEST_CASE("tmp") {
//    GeoLine line(-172.764801f, -31.8126984f, -172.765045f, -31.9921875f);
//    float NS = -31.8164062f;
//    uint32_t flags[3];
//    GeoLine lines[3];
//    cut_ns_ex(line, NS, flags, lines);
//
//    CAPTURE(flags[0]);
//    CAPTURE(flags[1]);
//    CAPTURE(lines[0]);
//    CAPTURE(lines[1]);
//    CHECK(false);
//}
