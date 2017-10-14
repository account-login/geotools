#include <set>
#include <map>
#include <utility>
#include <cassert>


namespace geotree {
    using namespace std;

    static const float LON_MAX = 180.0;
    static const float LON_MIN = -LON_MAX;
    static const float LAT_MAX = 85.0;
    static const float LAT_MIN = -LAT_MAX;

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

    enum { GEONODE_LEAF, GEONODE_INNER };

    template<class T>
    struct GeoNode {
        uint8_t type;
        uint32_t count;
        typedef map<T, GeoLonLat> MapType;
        MapType values;   // TODO: hash map

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


    template<class T>
    class GeoTree {
    private:
        GeoNode<T> *root;
        typedef map<T, GeoLonLat> MapType;
        MapType geos;   // TODO: hash map
        uint32_t split_threshold;
        uint32_t max_depth;

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

        size_t size() const {
            assert(this->geos.size() == node_size(this->root));
            return this->geos.size();
        }

        bool insert(const T &value, float lon, float lat) {
            assert(GeoLonLat(lon, lat).is_valid());
            
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

    private:
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

        void split(const GeoBox &box, GeoNode<T> *node) {
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

}   // namespace geotree
