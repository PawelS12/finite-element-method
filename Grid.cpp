#include "Grid.h"
#include <iostream>

using std::endl;
using std::cout;

Grid::Grid() {}

Grid::Grid(double p_nN, double p_nE, double p_nW, double p_nH, double p_height, double p_width) : nN(p_nN), nE(p_nE), nW(p_nW), nH(p_nH), height(p_height), width(p_width) {
    create_nodes();
    create_elements();
}

void Grid::create_nodes() {
    double delta_x = width / (nW - 1);
    double delta_y = height / (nH - 1);

    for (int i = 0; i < nH; ++i) {
        for (int j = 0; j < nW; ++j) {
            double x = i * delta_x;
            double y = j * delta_y;
            nodes.push_back(Node(x, y));
        }
    }
}

void Grid::create_elements() {
    for (int i = 0; i < nH - 1; ++i) {
        for (int j = 0; j < nW - 1; ++j) {
            int node_1 = i * nW + j + 1;
            int node_2 = node_1 + 1;
            int node_3 = node_1 + nW + 1;
            int node_4 = node_1 + nW;

            Element element;
            element.set_ID(node_1, node_2, node_3, node_4);
            elements.push_back(element);
        }
    }
}

void Grid::display_grid_data() {
    cout << "Nodes: " << endl << endl;
    for (size_t i = 0; i < nodes.size(); i++) {
        cout << "Node " << i + 1 << ": ";
        nodes[i].display_node();
    }
    cout << "-----------------------------------" << endl;
    cout << "Elements: " << endl << endl;
    for (size_t i = 0; i < elements.size(); i++) {
        cout << "Element " << i + 1 << ": ";
        elements[i].display_ID();
    }
    cout << "-----------------------------------" << endl;
}