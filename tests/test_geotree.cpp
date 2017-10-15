#include <vector>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>
#include "catch.h"

#define private public
#include "../geotree.hpp"


using namespace std;
using namespace geotree;

using boost::shared_ptr;
using boost::assign::list_of;


typedef int T;
typedef GeoNode<T> Node;
typedef GeoTree<T> Tree;


namespace std {
    ostream &operator<<(ostream &os, const Tree::Item &value) {
        os << "<Item " << value.value << " " << value.lon << " " << value.lat << " " << value.dist << ">";
        return os;
    }
}


Node *I(Node *nw, Node *ne, Node *se, Node *sw) {
    Node *ret = new Node(GEONODE_INNER);
    ret->NW = nw;
    ret->NE = ne;
    ret->SE = se;
    ret->SW = sw;
    return ret;
}

Node *L(const std::vector<T> &values) {
    Node *ret = new Node(GEONODE_LEAF);
    for (size_t i = 0; i < values.size(); ++i) {
        ret->add(values[i], GeoLonLat(i, i));
    }
    return ret;
}


void verify_shape(Node *input, Node *expect) {
    boost::shared_ptr<Node> gc(expect);

    if (input == NULL) {
        REQUIRE(expect == NULL);
        return;
    } else {
        REQUIRE(expect != NULL);
    }

    REQUIRE(input->type == expect->type);
    if (input->type == GEONODE_LEAF) {
        REQUIRE(input->values.size() == expect->values.size());

        Node::MapType::iterator in = input->values.begin();
        Node::MapType::iterator ex = expect->values.begin();
        for (; in != input->values.end(); ++in, ++ex) {
            CHECK(in->first == ex->first);
        }
    } else {
        verify_shape(input->NW, expect->NW);
        verify_shape(input->NE, expect->NE);
        verify_shape(input->SE, expect->SE);
        verify_shape(input->SW, expect->SW);
    }
}


void verify_tree(Tree &tree, Node *expect) {
    tree.verify();
    verify_shape(tree.root, expect);
}


TEST_CASE("util.dist") {
    CHECK(Tree::get_distance(-96.276111, 32.726386, -96.809261, 32.770455) == Approx(50114.9810144546));
    CHECK(Tree::get_distance(-111.382765, 39.205074, 133.617180, -26.496858) == Approx(13915095.4801221211));
    CHECK(Tree::get_distance(-111.382765, 39.2, -111.382765, 39.201) == Approx(111.2262999997));
    CHECK(Tree::get_distance(-180, 0, -90, 0) == Tree::get_distance(-180, 0, -90, 10));
}


TEST_CASE("basic.write") {
    Tree tree(3);
    verify_tree(tree, NULL);

    CHECK(tree.insert(123, -10, 20));
    verify_tree(tree, L(list_of(123)));

    CHECK(tree.insert(124, 10, 20));
    verify_tree(tree, L(list_of(123)(124)));

    CHECK(tree.insert(125, 10, -20));
    verify_tree(tree, L(list_of(123)(124)(125)));

    // split
    CHECK(tree.insert(126, -10, -20));
    verify_tree(tree, I(L(list_of(123)), L(list_of(124)), L(list_of(125)), L(list_of(126))));

    // update existing value
    CHECK(!tree.insert(123, -10, 21));
    verify_tree(tree, I(L(list_of(123)), L(list_of(124)), L(list_of(125)), L(list_of(126))));

    // move to another sub-tree
    CHECK(!tree.insert(123, 10, 21));
    verify_tree(tree, I(NULL, L(list_of(123)(124)), L(list_of(125)), L(list_of(126))));

    // split
    CHECK(tree.insert(131, -10, 21.0));
    CHECK(tree.insert(132, -10, 21.1));
    CHECK(tree.insert(133, -10, 21.2));
    CHECK(tree.insert(134, -10, 21.3));
    verify_tree(tree, I(
        I(NULL, NULL, L(list_of(131)(132)(133)(134)), NULL),
        L(list_of(123)(124)),
        L(list_of(125)),
        L(list_of(126))));

    // earse
    CHECK(!tree.erase(222));

    CHECK(tree.erase(132));
    verify_tree(tree, I(
        I(NULL, NULL, L(list_of(131)(133)(134)), NULL),
        L(list_of(123)(124)),
        L(list_of(125)),
        L(list_of(126))));
    CHECK(tree.erase(125));
    verify_tree(tree, I(
        I(NULL, NULL, L(list_of(131)(133)(134)), NULL),
        L(list_of(123)(124)),
        NULL,
        L(list_of(126))));
}


void sort_by_dist(vector<Tree::Item> &items, float lon, float lat) {
    for (size_t i = 0; i < items.size(); ++i) {
        Tree::Item &item = items[i];
        item.dist = Tree::round(Tree::get_distance(item.lon, item.lat, lon, lat));
    }

    sort(items.begin(), items.end());
}


void test_read(vector<Tree::Item> &items, Tree &tree) {
    for (size_t i = 0; i < items.size(); ++i) {
        float lon = items[i].lon;
        float lat = items[i].lat;

        vector<Tree::Item> sorted = items;
        sort_by_dist(sorted, lon, lat);

        vector<Tree::Item> got;
        vector<Tree::Item> expect;
        for (size_t count = 0; count < sorted.size(); ++count) {
            got = tree.get_nearby(lon, lat, count);
            expect.assign(sorted.begin(), sorted.begin() + count);
            CAPTURE(lon);
            CAPTURE(lat);
            CAPTURE(count);
            CHECK(got == expect);
        }

        got = tree.get_nearby(123, 45, 10000);
        CHECK(got.size() == sorted.size());
    }
}


void test_read_perm(const vector<Tree::Item> &data) {
    Tree tree(1);
    vector<Tree::Item> items;

    for (size_t i = 0; i < data.size(); ++i) {
        T val = data[i].value;
        float lon = data[i].lon;
        float lat = data[i].lat;

        tree.insert(val, lon, lat);
        tree.verify();

        items.push_back(Tree::Item(val, lon, lat));

        test_read(items, tree);
    }
}


TEST_CASE("basic.read") {
    vector<Tree::Item> data = list_of
        (Tree::Item(333, 100, 20))
        (Tree::Item(334, -80, -50))
        (Tree::Item(335, -80, -50.1))
        (Tree::Item(336, -80.05, -50))
        (Tree::Item(337, -80, 50.01))
        (Tree::Item(338, 0, 70))
        (Tree::Item(339, 180, 20))
        (Tree::Item(340, -180, 0));
    test_read_perm(data);
}


// TODO: stress test
