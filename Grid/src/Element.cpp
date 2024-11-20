#include "Element.h"
#include "Integration.h"
#include <cmath>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

Element::Element(const Node& n1, const Node& n2, const Node& n3, const Node& n4) {
    nodes_xy[0] = n1;
    nodes_xy[1] = n2;
    nodes_xy[2] = n3;
    nodes_xy[3] = n4;
    ID.resize(4, 0);
}

const Node* Element::get_nodes() const {
    return nodes_xy;
}

const std::vector<int>& Element::get_ID() const {
    return ID;
}

void Element::set_ID(int node_1, int node_2, int node_3, int node_4) {
    ID[0] = node_1;
    ID[1] = node_2;
    ID[2] = node_3;
    ID[3] = node_4;
    // cout << "Setting Element ID: " << node_1 << ", " << node_2 << ", " << node_3 << ", " << node_4 << endl;
}

void Element::display_ID() const {
    cout << "Element ID: ";
    for (size_t i = 0; i < ID.size(); i++) { 
        cout << ID[i] << " ";
    }
    cout << endl;
}