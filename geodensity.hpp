#pragma once

#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <float.h>

#include <boost/unordered_map.hpp>

#include "geotree.hpp"
#include "geoutil.hpp"


namespace geotools {
    using namespace std;

    class GeoDensity {
    public:
        typedef uint32_t KeyType;
        typedef uint32_t Distance;

        struct Stats {
            size_t perfect_hit;
            size_t guess_hit;
            size_t guess_total;
            size_t set_merged;
            size_t set_total;

            Stats()
                : perfect_hit(0), guess_hit(0), guess_total(0), set_merged(0), set_total(0)
            {}

            Stats reset() {
                Stats rv = *this;
                *this = Stats();
                return rv;
            }

            string repr() const {
                stringstream ss;
                double hit_ratio = double(guess_hit) / guess_total;
                ss << "[hit_ratio:" << hit_ratio << "][perfect_hit:" << perfect_hit << "]"
                   << "[guess_hit:" << guess_hit << "][guess_total:" << guess_total << "]"
                   << "[set_merged:" << set_merged << "][set_total:" << set_total << "]";
                return ss.str();
            }
        } stats;

    public:
        GeoDensity(Distance initial)
            : initial(initial), geotree(3)  // TODO: adjust geotree param
        {}

        // TODO: copy constructor

        Distance guess_radius(float lon, float lat) {
            this->stats.guess_total++;

            vector<GeoType::Item> items = this->geotree.get_nearby(lon, lat, 3);
            for (vector<GeoType::Item>::iterator it = items.begin(); it != items.end(); ) {
                KeyType key = it->value;
                assert(this->key2radius.count(key));
                Distance radius = this->key2radius[key];
                if (!is_valid_entry(radius, it->dist)) {
                    it = items.erase(it);
                } else {
                    ++it;
                }
            }

            if (items.empty()) {
                return this->initial;
            } else {
                this->stats.guess_hit++;

                double weight_sum = 0;
                double value_sum = 0;
                for (size_t i = 0; i < items.size(); ++i) {
                    double dist = items[i].dist;
                    Distance radius = this->key2radius[items[i].value];

                    if (dist <= 100) {
                        // too close
                        this->stats.perfect_hit++;
                        return radius;
                    }

                    double weight = 1.0 / (dist * dist);
                    weight_sum += weight;
                    value_sum += weight * radius;
                }

                return value_sum / weight_sum;
            }
        }

        KeyType set_radius(float lon, float lat, Distance radius) {
            this->stats.set_total++;

            // try to merge similar entries (re-use key)
            vector<GeoType::Item> items = this->geotree.get_nearby(lon, lat, 1);
            if (!items.empty()) {
                KeyType item_key = items[0].value;
                assert(this->key2radius.count(item_key));
                Distance item_radius = this->key2radius[item_key];
                if (is_similar_entry(item_radius, radius, items[0].dist)) {
                    this->stats.set_merged++;
                    return item_key;
                }
            }

            KeyType key = get_key();
            this->key2radius[key] = radius;
            this->geotree.insert(key, lon, lat);
            return key;
        }

        bool remove(KeyType key) {
            bool erased = this->key2radius.erase(key);
            if (erased) {
                bool geo_erased = this->geotree.erase(key);
                assert(geo_erased);
            }
            return erased;
        }

        size_t size() const {
            assert(this->geotree.size() == this->key2radius.size());
            return this->key2radius.size();
        }

    private:
        // XXX: key is not guanranteed to be unique
        KeyType get_key() {
            return this->seq++;
        }

        bool is_valid_entry(Distance radius, Distance dist) {
            // TODO: test this
            return dist < radius * 3;
        }

        bool is_similar_entry(Distance r1, Distance r2, Distance dist) {
            double rr = double(r1) / r2;
            if (rr > 1.2 || rr < 0.8) {
                return false;
            }

            double avg_r = (r1 + r2) / 2.0;
            if (dist > avg_r / 2.0) {
                return false;
            }

            return true;
        }

    private:
        Distance initial;

        typedef boost::unordered_map<KeyType, Distance> RadiusMap;
        RadiusMap key2radius;
        typedef GeoTree<Distance> GeoType;
        GeoType geotree;
        KeyType seq;
    };
}
