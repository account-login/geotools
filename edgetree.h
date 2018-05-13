#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <vector>

#include "geoutil.hpp"


namespace geotools {

    struct GeoLine {
        GeoLonLat src;
        GeoLonLat dst;

        GeoLine() {}
        GeoLine(float srclon, float srclat, float dstlon, float dstlat)
            : src(srclon, srclat), dst(dstlon, dstlat)
        {}

        bool operator==(const GeoLine &rhs) const {
            return this->src == rhs.src && this->dst == rhs.dst;
        }
    };

    enum { EDGENODE_LEAF, EDGENODE_INNER };

    struct EdgeNode {
        uint8_t type;

        // leaf
        std::vector<GeoLine> lines;

        // inner
        EdgeNode *NW;
        EdgeNode *NE;
        EdgeNode *SE;
        EdgeNode *SW;

        EdgeNode(uint8_t type) : type(type), NW(NULL), NE(NULL), SE(NULL), SW(NULL)
        {}

        void destroy() {
            if (NW) NW->destroy();
            if (NE) NE->destroy();
            if (SE) SE->destroy();
            if (SW) SW->destroy();
            delete this;
        }

        EdgeNode *&child(uint32_t flag) {
            switch (flag) {
            case D_NW: return NW;
            case D_NE: return NE;
            case D_SE: return SE;
            case D_SW: return SW;
            default:
                assert(!"Unreachable");
            }
        }
    };

    struct EdgeTree {
        uint32_t precision;
        EdgeNode *root;

        EdgeTree(uint32_t precision) : precision(precision), root(NULL) {}

        ~EdgeTree() {
            if (this->root != NULL) {
                this->root->destroy();
            }
        }

        void insert(const GeoLine &line);

    private:
        // no copy
        EdgeTree &operator=(const EdgeTree &rhs);
        EdgeTree(const EdgeTree &other);
    };

}   // ::geotools
