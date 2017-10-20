#pragma once

#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <float.h>

#include <boost/unordered_map.hpp>
#include <boost/thread/mutex.hpp>

#include "geotree.hpp"
#include "geoutil.hpp"


namespace geotools {
    using namespace std;

    class GeoDensity {
    public:
        typedef uint32_t KeyType;
        typedef uint32_t Distance;

        struct Entry {
            Distance radius;
            uint32_t count;
        };

        // TODO: update stats
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
        };
        mutable Stats stats;

    public:
        GeoDensity(Distance initial)
            : initial(initial), geotree(3)  // TODO: adjust geotree param
        {}

        // TODO: copy constructor

        Distance guess_radius(float lon, float lat, uint32_t count) const {
            this->stats.guess_total++;

            vector<GeoType::Item> items = this->geotree.get_nearby(lon, lat, 5);
            for (vector<GeoType::Item>::iterator it = items.begin(); it != items.end(); ) {
                KeyType key = it->value;
                assert(this->key2entry.count(key));
                if (!is_valid_entry(this->key2entry.find(key)->second, it->dist, count)) {
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
                    const Entry &e = this->key2entry.find(items[i].value)->second;
                    double ratio = (double)count / e.count;
                    double radius = sqrt(ratio) * e.radius;

                    double dist = items[i].dist;
                    dist = std::max(dist, 10.0);        // avoid divide by zero
                    double weight = 1 / (dist * dist) / (1 + 1.0 * abs(log2(ratio)));
                    weight_sum += weight;
                    value_sum += weight * radius;
                }

                return value_sum / weight_sum;
            }
        }

        KeyType set_radius(float lon, float lat, Distance radius, uint32_t count) {
            this->stats.set_total++;

            Entry e = {radius, count};
            // try to merge similar entries (re-use key)
            vector<GeoType::Item> items = this->geotree.get_nearby(lon, lat, 1);
            if (!items.empty()) {
                KeyType item_key = items[0].value;
                assert(this->key2entry.count(item_key));
                Entry near_entry = this->key2entry[item_key];
                if (is_similar_entry(e, near_entry, items[0].dist)) {
                    this->stats.set_merged++;
                    return item_key;
                }
            }

            KeyType key = get_key();
            this->key2entry.insert(make_pair(key, e));
            this->geotree.insert(key, lon, lat);
            return key;
        }

        bool remove(KeyType key) {
            bool erased = this->key2entry.erase(key);
            if (erased) {
                bool geo_erased = this->geotree.erase(key);
                assert(geo_erased);
            }
            return erased;
        }

        size_t size() const {
            assert(this->geotree.size() == this->key2entry.size());
            return this->key2entry.size();
        }

        Stats get_stats() const {
            return this->stats;
        }

        Stats pop_stats() const {
            return this->stats.reset();
        }

    private:
        // XXX: key is not guanranteed to be unique
        KeyType get_key() {
            return this->seq++;
        }

        static bool is_valid_entry(const Entry &e, Distance dist, uint32_t count) {
            // TODO: test this
            if (e.count > 10 * count || e.count < count / 10) {
                return false;
            }
            if (dist > e.radius * 3) {
                return false;
            }
            return true;
        }

        static bool is_similar_entry(const Entry &e1, const Entry &e2, Distance dist) {
            double r1 = e1.radius;
            double r2 = e2.radius;
            double rr = r1 / r2;
            if (rr > 1.2 || rr < 0.8) {
                return false;
            }

            double cr = double(e1.count) / e2.count;
            if (cr > 1.2 || cr < 0.8) {
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

        typedef boost::unordered_map<KeyType, Entry> EntryMap;
        EntryMap key2entry;
        typedef GeoTree<Distance> GeoType;
        GeoType geotree;
        KeyType seq;
    };


    template <class T>
    class GeoDensityTSAdaptor {
    public:
        typedef typename T::Distance Distance;
        typedef typename T::KeyType KeyType;
        typedef typename T::Stats Stats;

        GeoDensityTSAdaptor(Distance initial)
            : geoden(initial)
        {}

        Distance guess_radius(float lon, float lat, uint32_t count) const {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.guess_radius(lon, lat, count);
        }

        KeyType set_radius(float lon, float lat, Distance radius, uint32_t count) {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.set_radius(lon, lat, radius, count);
        }

        bool remove(KeyType key) {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.remove(key);
        }

        size_t size() const {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.size();
        }

        Stats get_stats() const {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.get_stats();
        }

        Stats pop_stats() const {
            boost::mutex::scoped_lock sl(this->mutex);
            return geoden.pop_stats();
        }

        T &unwrap() {
            return geoden;
        }

    private:
        T geoden;
        mutable boost::mutex mutex;
    };

    typedef GeoDensityTSAdaptor<GeoDensity> GeoDensityTS;
}
