#pragma once

#include <list>

#include <boost/unordered_map.hpp>


namespace geotools {

    template <class T>
    class LRUSet {
    private:
        typedef std::list<T> ListType;
        typedef boost::unordered_map<T, typename ListType::iterator> MapType;

        ListType order;     // front: new, back: old
        MapType values;

    public:
        bool insert(const T &value) {
            typename MapType::iterator it = this->values.find(value);
            if (it == this->values.end()) {
                this->order.push_front(value);
                this->values.insert(std::make_pair(value, this->order.begin()));
                return true;
            } else {
                typename ListType::iterator lit = it->second;
                this->order.erase(lit);
                this->order.push_front(value);
                it->second = this->order.begin();
                return false;
            }
        }

        bool remove(const T &value) {
            typename MapType::iterator it = this->values.find(value);
            if (it == this->values.end()) {
                return false;
            }

            this->order.erase(it->second);
            this->values.erase(it);
            return true;
        }

        size_t size() const {
            return this->values.size();
        }

        bool empty() const {
            return this->size() == 0;
        }

        T pop() {
            assert(!this->empty());
            T val = this->order.back();
            bool erased = this->values.erase(val);
            assert(erased);
            this->order.pop_back();
            return val;
        }
    };

}   // namespace geotools
