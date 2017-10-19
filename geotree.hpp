#pragma once

#include <algorithm>
#include <cfloat>
#include <stdint.h>
#include <set>
#include <vector>
#include <map>
#include <math.h>
#include <utility>
#include <cassert>

#include <boost/unordered_map.hpp>

#include "geoutil.hpp"


namespace geotools {
    using namespace std;

    // direction
    enum {
        D_NONE = 0,
        D_W = 1 << 0,
        D_E = 1 << 1,
        D_N = 1 << 2,
        D_S = 1 << 3,
        D_NW = D_N | D_W,
        D_NE = D_N | D_E,
        D_SE = D_S | D_E,
        D_SW = D_S | D_W,
    };


    struct GeoLonLat {
        float lon;
        float lat;

        GeoLonLat(float lon, float lat) : lon(lon), lat(lat) {}

        bool operator==(const GeoLonLat &rhs) const {
            return lon == rhs.lon && lat == rhs.lat;
        }

        bool is_valid() {
            return LON_MIN <= lon && lon <= LON_MAX && LAT_MIN <= lat && lat <= LAT_MAX;
        }
    };

    // node type
    enum { GEONODE_LEAF, GEONODE_INNER };

    template<class T>
    struct GeoNode {
        uint8_t type;
        uint32_t count;
        typedef boost::unordered_map<T, GeoLonLat> MapType;
        MapType values;     // TODO: try skip list here

        GeoNode *NW;
        GeoNode *NE;
        GeoNode *SE;
        GeoNode *SW;

        GeoNode(uint8_t type) : type(type), count(0), NW(NULL), NE(NULL), SE(NULL), SW(NULL) {}

        void destroy() {
            if (NW) NW->destroy();
            if (NE) NE->destroy();
            if (SE) SE->destroy();
            if (SW) SW->destroy();
            delete this;
        }

        bool is_leaf() const {
            if (this->type == GEONODE_LEAF) {
                assert(NW == NULL);
                assert(NE == NULL);
                assert(SE == NULL);
                assert(SW == NULL);
            }
            return this->type == GEONODE_LEAF;
        }

        bool add(const T &value, GeoLonLat lonlat) {
            assert(this->is_leaf());
            pair<typename MapType::iterator, bool> itok = this->values.insert(make_pair(value, lonlat));
            if (itok.second) {
                this->count++;
            }
            assert(this->count == this->values.size());
            return itok.second;
        }

        void must_remove(const T &value) {
            assert(this->is_leaf());
            size_t removed = this->values.erase(value);
            assert(removed == 1);
            this->count--;
            assert(this->count == this->values.size());
        }

        void update_count() {
            assert(!this->is_leaf());
            this->count = count_of(this->NW) + count_of(this->NE) + count_of(this->SE) + count_of(this->SW);
        }

        GeoNode<T> *&get(int dir) {
            switch (dir) {
            case D_NW: return this->NW;
            case D_NE: return this->NE;
            case D_SE: return this->SE;
            case D_SW: return this->SW;
            default:
                assert(!"unreachable");
            }
        }

        const GeoNode<T> *get(int dir) const {
            switch (dir) {
            case D_NW: return this->NW;
            case D_NE: return this->NE;
            case D_SE: return this->SE;
            case D_SW: return this->SW;
            default:
                assert(!"unreachable");
            }
        }

        static uint32_t count_of(const GeoNode<T> *node) {
            return node != NULL ? node->count : 0;
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
            case D_NW: return GeoBox(W, (W + E) / 2.0, N, (N + S) / 2.0);
            case D_NE: return GeoBox((W + E) / 2.0, E, N, (N + S) / 2.0);
            case D_SE: return GeoBox((W + E) / 2.0, E, (N + S) / 2.0, S);
            case D_SW: return GeoBox(W, (W + E) / 2.0, (N + S) / 2.0, S);
            default:
                assert(!"unreachable");
            }
        }

        int locate(GeoLonLat lonlat) const {
            assert(this->contains(lonlat));

            int dir = D_NONE;
            
            float WE = (W + E) / 2.0;
            if (lonlat.lon < WE) {
                dir |= D_W;
            } else {
                dir |= D_E;
            }

            float NS =  (N + S) / 2.0;
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

            float WE = (W + E) / 2.0;
            if (lonlat.lon < WE) {
                dir |= D_W;
                E = WE;
            } else {
                dir |= D_E;
                W = WE;
            }

            float NS =  (N + S) / 2.0;
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


    enum GeoOption {
        GEO_OPT_NONE = 0,
        GEO_NO_SORT = 1 << 0,
    };


    template<class T>
    class GeoTree {
    private:
        typedef GeoNode<T> Node;
        typedef Node *NodePtr;

        GeoNode<T> *root;
        typedef boost::unordered_map<T, GeoLonLat> MapType;
        MapType geos;
        uint32_t split_threshold;
        uint32_t max_depth;     // TODO: setter and getter

        static size_t node_size(const GeoNode<T> *node) {
            return node != NULL ? node->count : 0;
        }

        struct GeoInsertCtx {
            T value;
            GeoLonLat lonlat;
            GeoBox box;
            uint32_t depth;
    
            GeoInsertCtx(const T &value, GeoLonLat lonlat)
                : value(value), lonlat(lonlat), box(), depth(0)
            {}
        };
    
    public:
        GeoTree(uint32_t split_threshold = 128)
            : root(NULL), geos(), split_threshold(split_threshold), max_depth(16)   // less than 1km
        {}

        ~GeoTree() {
            if (this->root != NULL) {
                this->root->destroy();
            }
        }

        // TODO: copy constructor

        struct Item {
            T value;
            float lon;
            float lat;
            uint32_t dist;

            Item(const T &value, float lon, float lat)
                : value(value), lon(lon), lat(lat), dist(0)
            {}

            Item()
                : value(), lon(FLT_MAX), lat(FLT_MAX), dist(~uint32_t(0))
            {}

            bool operator<(const Item &rhs) const {
                return this->dist < rhs.dist;
            }

            bool operator==(const Item &rhs) const {
                return value == rhs.value
                    && lon == rhs.lon && lat == rhs.lat && dist == rhs.dist;
            }
        };

        size_t size() const {
            assert(this->geos.size() == node_size(this->root));
            return this->geos.size();
        }

        static bool is_valid(float lon, float lat) {
            return GeoLonLat(lon, lat).is_valid();
        }

        bool insert(const T &value, float lon, float lat) {
            assert(is_valid(lon, lat));

            GeoInsertCtx ctx(value, GeoLonLat(lon, lat));

            typename MapType::iterator it = this->geos.find(value);
            bool exists = it != this->geos.end();
            if (exists) {
                GeoInsertCtx rem_ctx(value, it->second);
                this->root = this->remove_rec(rem_ctx, this->root);
                it->second = ctx.lonlat;
            } else {
                this->geos.insert(make_pair(value, ctx.lonlat));
            }

            this->root = this->insert_rec(ctx, this->root);
            return !exists;
        }

        bool insert(const Item &item) {
            return insert(item.value, item.lon, item.lat);
        }

        bool erase(const T &value) {
            typename MapType::iterator it = this->geos.find(value);
            if (it == this->geos.end()) {
                return false;
            }

            GeoInsertCtx ctx(value, it->second);
            this->root = this->remove_rec(ctx, this->root);
            this->geos.erase(it);
            return true;
        }

        vector<Item> get_nearby(
            float lon, float lat, size_t count, uint32_t option = GEO_OPT_NONE) const
        {
            return nearby_impl(GeoLonLat(lon, lat), count, option);
        }

        // TODO: optimize
        uint32_t get_nearby_radius_by_count(float lon, float lat, size_t count) const {
            vector<Item> items = get_nearby(lon, lat, count);
            if (items.empty()) {
                return 0;
            } else {
                return items.back().dist;
            }
        }

    private:
        struct NineBox {
            typedef const Node *ConstNodePtr;
            ConstNodePtr NW, N, NE, W, C, E, SW, S, SE;

            NineBox moved(int dir) const {
                switch (dir) {
                case D_NW: return NineBox {
                        get_sub_node(NW, D_SE), // NW
                        get_sub_node(N,  D_SW), // N
                        get_sub_node(N,  D_SE), // NE
                        get_sub_node(W,  D_NE), // W
                        get_sub_node(C,  D_NW), // C
                        get_sub_node(C,  D_NE), // E
                        get_sub_node(W,  D_SE), // SW
                        get_sub_node(C,  D_SW), // S
                        get_sub_node(C,  D_SE), // SE
                    };
                case D_NE: return NineBox {
                        get_sub_node(N,  D_SW), // NW
                        get_sub_node(N,  D_SE), // N
                        get_sub_node(NE, D_SW), // NE
                        get_sub_node(C,  D_NW), // W
                        get_sub_node(C,  D_NE), // C
                        get_sub_node(E,  D_NW), // E
                        get_sub_node(C,  D_SW), // SW
                        get_sub_node(C,  D_SE), // S
                        get_sub_node(E,  D_SW), // SE
                    };
                case D_SE: return NineBox {
                        get_sub_node(C,  D_NW), // NW
                        get_sub_node(C,  D_NE), // N
                        get_sub_node(E,  D_NW), // NE
                        get_sub_node(C,  D_SW), // W
                        get_sub_node(C,  D_SE), // C
                        get_sub_node(E,  D_SW), // E
                        get_sub_node(S,  D_NW), // SW
                        get_sub_node(S,  D_NE), // S
                        get_sub_node(SE, D_NW), // SE
                    };
                case D_SW: return NineBox {
                        get_sub_node(W,  D_NE), // NW
                        get_sub_node(C,  D_NW), // N
                        get_sub_node(C,  D_NE), // NE
                        get_sub_node(W,  D_SE), // W
                        get_sub_node(C,  D_SW), // C
                        get_sub_node(C,  D_SE), // E
                        get_sub_node(SW, D_NE), // SW
                        get_sub_node(S,  D_NW), // S
                        get_sub_node(S,  D_NE), // SE
                    };
                default:
                    assert(!"Unreachable");
                }
            }

            static ConstNodePtr get_sub_node(ConstNodePtr node, int dir) {
                if (node == NULL) {
                    return NULL;
                }
                return node->get(dir);
            }
        };

        vector<Item> nearby_impl(GeoLonLat lonlat, uint32_t count, uint32_t option) const {
            if (count == 0 || this->root == NULL) {
                vector<Item> empty;
                return empty;
            }

            const Node *const r = this->root;
            NineBox ninebox = {r, r, r, r, r, r, r, r, r};
            GeoBox box;

            // find smallest center node that is guaranteed to cover required count
            while (true) {
                assert(ninebox.C != NULL);
                if (ninebox.C->is_leaf()) {
                    break;
                }

                int dir = box.locate_and_move(lonlat);
                NineBox new_ninebox = ninebox.moved(dir);
                if (node_size(new_ninebox.C) < count) {     // new_ninebox.C may be NULL
                    break;
                }

                ninebox = new_ninebox;
            }

            // remove duplicated node
            const Node **arr = &ninebox.NW;
            for (size_t i = 1; i < 9; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    if (arr[i] == arr[j]) {
                        arr[i] = NULL;
                        break;
                    }
                }
            }

            return fetch_items(arr, 9, lonlat, count, option);
        }

        static vector<Item> fetch_items(
            const Node **arr, size_t size, GeoLonLat lonlat, uint32_t count, uint32_t option)
        {
            vector<Item> ans;

            for (size_t i = 0; i < size; ++i) {
                collect_item(arr[i], ans);
            }

            measure_distance(ans, lonlat);
            truncate_by_distance(ans, count);
            if ((option & GEO_NO_SORT) == 0) {
                sort(ans.begin(), ans.end());
            }
            return ans;
        }

        static void collect_item(const Node *node, vector<Item> &data) {
            if (node == NULL) {
                return;
            }

            if (node->is_leaf()) {
                for (typename Node::MapType::const_iterator it = node->values.begin();
                    it != node->values.end(); ++it)
                {
                    const GeoLonLat &lonlat = it->second;
                    data.push_back(Item(it->first, lonlat.lon, lonlat.lat));
                }
            } else {
                collect_item(node->NW, data);
                collect_item(node->NE, data);
                collect_item(node->SE, data);
                collect_item(node->SW, data);
            }
        }

        static void measure_distance(vector<Item> &data, GeoLonLat lonlat) {
            for (size_t i = 0; i < data.size(); ++i) {
                Item &item = data[i];
                item.dist = geo_round(geo_distance(item.lon, item.lat, lonlat.lon, lonlat.lat));
            }
        }

        static void truncate_by_distance(vector<Item> &data, uint32_t count) {
            // partition by count and truncate
            assert(count > 0);
            if (data.size() > count) {
                nth_element(data.begin(), data.begin() + count - 1, data.end());
                data.resize(count);
            }
        }

        GeoNode<T> *insert_rec(GeoInsertCtx &ctx, GeoNode<T> *node) {
            if (node == NULL) {
                node = new GeoNode<T>(GEONODE_LEAF);
            }

            if (node->is_leaf()) {
                node->add(ctx.value, ctx.lonlat);
                if (node->count > this->split_threshold && ctx.depth < this->max_depth) {
                    this->split(ctx.box, node);
                }
            } else {
                ctx.depth++;
                int dir = ctx.box.locate_and_move(ctx.lonlat);
                GeoNode<T> *&c = node->get(dir);
                c = insert_rec(ctx, c);
                node->update_count();
            }

            return node;
        }

        static GeoNode<T> *leaf_add(GeoNode<T> *node, const T &value, GeoLonLat lonlat) {
            if (node == NULL) {
                node = new GeoNode<T>(GEONODE_LEAF);
            }
            node->add(value, lonlat);
            return node;
        }

        static void split(const GeoBox &box, GeoNode<T> *node) {
            assert(node->is_leaf());
            for (typename GeoNode<T>::MapType::iterator it = node->values.begin();
                it != node->values.end(); ++it)
            {
                const T &key = it->first;
                const GeoLonLat &lonlat = it->second;

                int dir = box.locate(lonlat);
                GeoNode<T> *&c = node->get(dir);
                c = leaf_add(c, key, lonlat);
            }
            node->type = GEONODE_INNER;
            node->values.clear();
        }

        GeoNode<T> *remove_rec(GeoInsertCtx &ctx, GeoNode<T> *node) {
            assert(node != NULL);

            if (node->is_leaf()) {
                node->must_remove(ctx.value);
                // remove empty node
                if (node->count == 0) {
                    delete node;
                    node = NULL;
                }
                // TODO: merge
            } else {
                int dir = ctx.box.locate_and_move(ctx.lonlat);
                GeoNode<T> *&c = node->get(dir);
                c = remove_rec(ctx, c);
                node->update_count();
                // remove node when count reaches 0
                if (node->count == 0) {
                    assert(!node->NW && !node->NE && !node->SE && !node->SW);
                    delete node;
                    node = NULL;
                }
            }

            return node;
        }

    // for tests
#ifdef TWOBLUECUBES_SINGLE_INCLUDE_CATCH_HPP_INCLUDED
    public:
        void verify() const {
            this->size();
            this->verify_node(this->root, GeoBox());
        }

        MapType get_all() const {
            return this->geos;
        }

    private:
        void verify_node(const GeoNode<T> *node, const GeoBox &box) const {
            if (node == NULL) {
                return;
            }

            if (node->type == GEONODE_LEAF) {
                CHECK(node->count == node->values.size());
                for (typename GeoNode<T>::MapType::const_iterator it = node->values.begin();
                    it != node->values.end(); ++it)
                {
                    const T &key = it->first;
                    const GeoLonLat &lonlat = it->second;

                    CHECK(box.contains(lonlat));
                    typename MapType::const_iterator geoit = this->geos.find(key);
                    CHECK(geoit != this->geos.end());
                    CHECK(geoit->second == lonlat);
                }
            } else {
                CHECK(node->count
                    == node_size(node->NW) + node_size(node->NE)
                    + node_size(node->SE) + node_size(node->SW));
                this->verify_node(node->NW, box.get(D_NW));
                this->verify_node(node->NE, box.get(D_NE));
                this->verify_node(node->SE, box.get(D_SE));
                this->verify_node(node->SW, box.get(D_SW));
            }
        }
#endif
    };

}   // namespace geotools
