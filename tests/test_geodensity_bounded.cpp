#include "catch.h"

#include "../geodensity_bounded.hpp"


using namespace std;
using namespace geotools;


typedef GeoDensityBounded::KeyType Key;


TEST_CASE("basic") {
    const uint32_t count = 1001;
    const uint32_t default_radius = 12345;
    GeoDensityBounded den(default_radius, 3);

    Key k1 = den.set_radius(123, 23, 100 * 1000, count);
    Key k2 = den.set_radius(124, 23, 110 * 1000, count);
    Key k3 = den.set_radius(123, 22, 105 * 1000, count);
    // k1 poped
    Key k4 = den.set_radius(123, 24, 130 * 1000, count);

    CHECK(den.size() == 3);
    CHECK(!den.remove(k1));

    Key k5 = den.set_radius(123, 22.1, 106 * 1000, count);
    CHECK(k5 == k3);
    den.set_radius(123.1, 22.5, 118 * 1000, count);

    // k2 popped
    CHECK(den.size() == 3);
    CHECK(!den.remove(k2));

    // query
    uint32_t r1 = den.guess_radius(123.5, 23, count);
    CHECK(r1 > 105 * 1000);
    CHECK(r1 < 130 * 1000);
}
