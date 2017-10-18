#pragma once

#include "geodensity.hpp"
#include "lruset.hpp"


namespace geotools {

    class GeoDensityBounded : public GeoDensity {
    public:
        GeoDensityBounded(Distance initial, size_t limit = 0)
            : GeoDensity(initial), limit(limit)
        {}

        KeyType set_radius(float lon, float lat, Distance radius) {
            // remove old entries
            while (this->limit != 0 && this->size() >= this->limit) {
                assert(this->size() == this->lru.size());
                bool erased = this->remove(this->lru.pop());
                assert(erased);
            }

            KeyType key = GeoDensity::set_radius(lon, lat, radius);
            this->lru.insert(key);
            return key;
        }

        void set_limit(size_t limit) {
            this->limit = limit;
        }

        size_t get_limit() const {
            return this->limit;
        }

    private:
        size_t limit = 0;
        LRUSet<KeyType> lru;
    };

}
