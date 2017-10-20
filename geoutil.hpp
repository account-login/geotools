#pragma once

#include <math.h>


namespace geotools {

    static const float LON_MAX = 180.0;
    static const float LON_MIN = -LON_MAX;
    static const float LAT_MAX = 85.0;
    static const float LAT_MIN = -LAT_MAX;

    static const double EARTH_RADIUS_IN_METERS = 6372797.560856;

// shit
#ifndef M_PI
    # define M_PI           3.14159265358979323846  /* pi */
#endif


    inline double deg2rad(double deg) {
        return deg / 180.0 * M_PI;
    }

    inline double geo_distance(double lon1d, double lat1d, double lon2d, double lat2d) {
        double lat1r = deg2rad(lat1d);
        double lon1r = deg2rad(lon1d);
        double lat2r = deg2rad(lat2d);
        double lon2r = deg2rad(lon2d);
        double u = sin((lat2r - lat1r) / 2);
        double v = sin((lon2r - lon1r) / 2);
        return 2.0 * EARTH_RADIUS_IN_METERS *
            asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
    }

    inline int32_t geo_round(double flt) {
        return ceil(flt - 0.5);
    }

}
