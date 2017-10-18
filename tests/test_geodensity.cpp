#include "catch.h"

#include "../geodensity.hpp"


using namespace std;
using namespace geotools;


typedef GeoDensity::KeyType Key;


TEST_CASE("basic") {
    const uint32_t default_radius = 12345;
    GeoDensity den(default_radius);
    CHECK(den.size() == 0);

    // not hit
    CHECK(den.guess_radius(123, 24) == default_radius);

    // pefect hit
    Key key = den.set_radius(123, 24, 100 * 1000);
    CHECK(den.size() == 1);
    CHECK(den.guess_radius(123, 24) == Approx(100 * 1000));

    // not hit
    CHECK(den.guess_radius(100, 24) == default_radius);

    // remove
    CHECK(den.remove(key));
    CHECK(den.size() == 0);
    CHECK(!den.remove(888666));

    // not hit
    CHECK(den.guess_radius(123, 24) == default_radius);

    Key k1 = den.set_radius(123, 24, 100 * 1000);
    den.set_radius(124, 25, 130 * 1000);
    den.set_radius(124, 23, 80 * 1000);
    CHECK(den.size() == 3);

    // similar entry
    Key k4 = den.set_radius(123, 24.1, 102 * 1000);
    CHECK(den.size() == 3);
    CHECK(k1 == k4);

    // avg of 3
    uint32_t r1 = den.guess_radius(123.5, 24);
    CHECK(r1 > 80 * 1000);
    CHECK(r1 < 130 * 1000);

    CHECK(den.remove(k4));
    den.set_radius(110, 30, 500 * 10000);

    // avg of 2
    uint32_t r2 = den.guess_radius(123.5, 24);
    CHECK(r2 > 80 * 1000);
    CHECK(r2 < 130 * 1000);
    CHECK(r2 != Approx(r1));

    CHECK(den.stats.reset().repr() != "");
}
