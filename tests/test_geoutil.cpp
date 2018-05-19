#include "doctest.h"

#include "../geoutil.hpp"


using namespace std;
using namespace geotools;

using doctest::Approx;


TEST_CASE("util") {
    CHECK(geo_distance(-96.276111, 32.726386, -96.809261, 32.770455) == Approx(50114.9810144546));
    CHECK(geo_distance(-111.382765, 39.205074, 133.617180, -26.496858) == Approx(13915095.4801221211));
    CHECK(geo_distance(-111.382765, 39.2, -111.382765, 39.201) == Approx(111.2262999997));
    CHECK(geo_distance(-180, 0, -90, 0) == geo_distance(-180, 0, -90, 10));
}
