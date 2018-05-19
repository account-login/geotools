#pragma once

#include <stdint.h>
#include <math.h>
#include <cassert>


namespace geotools {

    static const float LON_MAX = 180.0;
    static const float LON_MIN = -LON_MAX;
    static const float LAT_MAX = 90.0;
    static const float LAT_MIN = -LAT_MAX;

    static const float GEOTREE_LAT_MAX = 85.0;
    static const float GEOTREE_LAT_MIN = -GEOTREE_LAT_MAX;

    static const double EARTH_RADIUS_IN_METERS = 6372797.560856;

// shit
#ifndef M_PI
#   define M_PI           3.14159265358979323846    /* pi */
#endif


    enum Direction {
        D_NONE = 0,
        D_W = 1 << 0,   // 1
        D_E = 1 << 1,   // 2
        D_N = 1 << 2,   // 4
        D_S = 1 << 3,   // 8
        D_NW = D_N | D_W,
        D_NE = D_N | D_E,
        D_SE = D_S | D_E,
        D_SW = D_S | D_W,
    };

    struct GeoLonLat {
        float lon;  // alpha
        float lat;  // beta

        GeoLonLat(float lon, float lat) : lon(lon), lat(lat) {}
        GeoLonLat() : lon(0), lat(0) {}

        bool operator==(const GeoLonLat &rhs) const {
            return lon == rhs.lon && lat == rhs.lat;
        }

        bool operator!=(const GeoLonLat &rhs) const {
            return !(*this == rhs);
        }

        bool is_valid() const {
            return LON_MIN <= lon && lon <= LON_MAX && LAT_MIN <= lat && lat <= LAT_MAX;
        }
    };

    struct GeoBox {
        float W, E, N, S;

        GeoBox() : W(LON_MIN), E(LON_MAX), N(LAT_MAX), S(LAT_MIN) {}
        GeoBox(float w, float e, float n, float s) : W(w), E(e), N(n), S(s) {}

        bool contains(GeoLonLat lonlat) const {
            return W <= lonlat.lon && lonlat.lon <= E && S <= lonlat.lat && lonlat.lat <= N;
        }

        GeoBox get(int dir) const {
            switch (dir) {
            case D_NW: return GeoBox(W, (W + E) / 2.0f, N, (N + S) / 2.0f);
            case D_NE: return GeoBox((W + E) / 2.0f, E, N, (N + S) / 2.0f);
            case D_SE: return GeoBox((W + E) / 2.0f, E, (N + S) / 2.0f, S);
            case D_SW: return GeoBox(W, (W + E) / 2.0f, (N + S) / 2.0f, S);
            default:
                assert(!"unreachable");
            }
        }

        int locate(GeoLonLat lonlat) const {
            assert(this->contains(lonlat));

            int dir = D_NONE;

            float WE = (W + E) / 2.0f;
            if (lonlat.lon < WE) {
                dir |= D_W;
            } else {
                dir |= D_E;
            }

            float NS =  (N + S) / 2.0f;
            if (lonlat.lat < NS) {
                dir |= D_S;
            } else {
                dir |= D_N;
            }

            return dir;
        }

        int locate_and_move(GeoLonLat lonlat) {
            assert(this->contains(lonlat));

            int dir = D_NONE;

            float WE = (W + E) / 2.0f;
            if (lonlat.lon < WE) {
                dir |= D_W;
                E = WE;
            } else {
                dir |= D_E;
                W = WE;
            }

            float NS =  (N + S) / 2.0f;
            if (lonlat.lat < NS) {
                dir |= D_S;
                N = NS;
            } else {
                dir |= D_N;
                S = NS;
            }

            return dir;
        }
    };


    inline double deg2rad(double deg) {
        return deg / 180.0 * M_PI;
    }

    inline double rad2deg(double rad) {
        return rad * 180.0 / M_PI;
    }

    inline double geo_angle(double lon1d, double lat1d, double lon2d, double lat2d) {
        double lat1r = deg2rad(lat1d);
        double lon1r = deg2rad(lon1d);
        double lat2r = deg2rad(lat2d);
        double lon2r = deg2rad(lon2d);
        double u = sin((lat2r - lat1r) / 2);
        double v = sin((lon2r - lon1r) / 2);
        return 2.0 * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));   // rad
    }

    inline double geo_distance(double lon1d, double lat1d, double lon2d, double lat2d) {
        return geo_angle(lon1d, lat1d, lon2d, lat2d) * EARTH_RADIUS_IN_METERS;
    }

    inline int32_t geo_round(double flt) {
        return (int32_t)ceil(flt - 0.5);
    }

}   // ::geotools
