#include <vector>

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
    if (input == NULL) {
        REQUIRE(expect == NULL);
        return;
    } else {
        REQUIRE(expect != NULL);
    }

    boost::shared_ptr<Node> gc(expect);

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


TEST_CASE("basic") {
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
}
