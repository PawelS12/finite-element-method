#include "Grid.h"
#include <iostream>
#include <vector>
#include <math.h>

using std::endl;
using std::cout;
using std::cerr;
using std::vector;

Grid::Grid() {}

Grid::Grid(double p_nN, double p_nE, double p_nW, double p_nH, double p_height, double p_width) : nN(p_nN), nE(p_nE), nW(p_nW), nH(p_nH), height(p_height), width(p_width) {
    create_elements();
}

vector<Node> Grid::nodes_xy;

const vector<Element>& Grid::get_elements() const {
     return elements; 
}

void Grid::create_elements() {
    nodes_xy = GlobalData::get_nodes();  
    size_t nodes_num = nodes_xy.size();

    if (nodes_num != nN) {  
        cerr << "Wrong number of nodes." << endl;
        return;
    }

    int nodes_per_row = sqrt(nN); 
    int nodes_per_col = nN / nodes_per_row;  

    for (int i = 0; i < nodes_per_col - 1; ++i) {
        for (int j = 0; j < nodes_per_row - 1; ++j) {
            int node_1 = i * nodes_per_row + j;   
            int node_2 = node_1 + 1;              
            int node_3 = node_1 + nodes_per_row + 1;   
            int node_4 = node_1 + nodes_per_row;     

            Element new_element(nodes_xy[node_1], nodes_xy[node_2], nodes_xy[node_3], nodes_xy[node_4]);
            new_element.set_ID(node_1, node_2, node_3, node_4);  

            elements.push_back(new_element);  
        }
    }
}


void Grid::display_grid_data() {
    cout << "Nodes:" << endl << endl;
    for (const auto& node : nodes_xy) {
        node.display_node();
    }
    cout << "-----------------------------------" << endl;
    cout << "Elements:" << endl << endl;
    for (const auto& element : elements) {  
        element.display_ID(); 
    }
    cout << "-----------------------------------" << endl;   
}