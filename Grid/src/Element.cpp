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
    Hbc.resize(4, vector<double>(4, 0.0));
    P.resize(4, 0.0);
}

const Node* Element::get_nodes() const {
    return nodes_xy;
}

const vector<int>& Element::get_ID() const {
    return ID;
}

void Element::set_ID(int node_1, int node_2, int node_3, int node_4) {
    ID[0] = node_1;
    ID[1] = node_2;
    ID[2] = node_3;
    ID[3] = node_4;
}

void Element::display_ID() const {
    cout << "Element ID: ";
    for (size_t i = 0; i < ID.size(); i++) { 
        cout << ID[i] << " ";
    }
    cout << endl;
}

void Element::set_Hbc(const vector<vector<double>>& hbc) {
     Hbc = hbc; 
}

const vector<vector<double>>& Element::get_Hbc() const { 
    return Hbc; 
}

void Element::set_P(const vector<double>& p) {
    P = p;
}

const vector<double>& Element::get_P() const {
    return P;
}