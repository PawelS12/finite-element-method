#include "Element.h"
#include <iostream>

using std::cout;
using std::endl;

Element::Element() {
    ID.resize(4, 0);
}

void Element::set_ID(int node_1, int node_2, int node_3, int node_4) {
    ID[0] = node_1;
    ID[1] = node_2;
    ID[2] = node_3;
    ID[3] = node_4;
}

void Element::display_ID() {
    for (size_t i = 0; i < ID.size(); i++) { 
        cout << ID[i] << " ";
    }
    cout << endl;
}