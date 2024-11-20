#include "Grid.h"
#include <iostream>
#include <vector>

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
    size_t nodes_num = GlobalData::get_nodes().size();
  
    if (nodes_num != 16) {
        cerr << "Wrong number of nodes." << endl;
        return;
    }

    for (int i = 0; i < nH - 1; ++i) { 
        for (int j = 0; j < nW - 1; ++j) {
            int node_1 = i * nW + j;        
            int node_2 = node_1 + 1;        
            int node_3 = node_1 + nW + 1;    
            int node_4 = node_1 + nW;       

            //Element new_element(nodes_xy[node_1], nodes_xy[node_2], nodes_xy[node_3], nodes_xy[node_4]);
            //new_element.set_ID(node_1, node_2, node_3, node_4); 

            Element new_element(nodes_xy[node_3], nodes_xy[node_4], nodes_xy[node_1], nodes_xy[node_2]);
            new_element.set_ID(node_3, node_4, node_1, node_2); 
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