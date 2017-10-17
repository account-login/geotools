#include <vector>

#include "catch.h"

#include "../lruset.hpp"


using namespace std;
using namespace geotools;


typedef int T;
typedef LRUSet<T> LRU;


vector<T> flat_lruset(LRU &set) {
    vector<T> ans;
    while (!set.empty()) {
        ans.push_back(set.pop());
    }
    return ans;
}


void inserts(LRU &set, const vector<T> &elems) {
    for (T e : elems) {
        set.insert(e);
    }
}


TEST_CASE("basic") {
    LRU lru;
    vector<T> expect;
    CHECK(flat_lruset(lru) == expect);

    CHECK(lru.insert(1));
    CHECK(lru.size() == 1);
    expect = {1};
    CHECK(flat_lruset(lru) == expect);

    lru.insert(1);
    CHECK(lru.insert(2));
    CHECK(lru.size() == 2);
    expect = {1, 2};
    CHECK(flat_lruset(lru) == expect);

    lru.insert(1);
    lru.insert(2);
    CHECK(!lru.insert(1));
    CHECK(lru.size() == 2);
    expect = {2, 1};
    CHECK(flat_lruset(lru) == expect);

    lru.insert(1);
    lru.insert(2);
    lru.insert(4);
    CHECK(!lru.remove(3));
    CHECK(lru.remove(2));
    expect = {1, 4};
    CHECK(flat_lruset(lru) == expect);
}
