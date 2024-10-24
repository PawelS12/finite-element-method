#ifndef GRID_H
#define GRID_H

#include <vector>
#include "Node.h"
#include "Elem4.h"
#include "Element.h"

using std::vector;

class Grid {
private:
    double nN, nE, nW, nH, height, width;
    vector<Node> nodes;
    vector<Node> nodes_xy;
    vector<Element> elements;

public:
    Grid();
    Grid(double p_nN, double p_nE, double p_nW, double p_nH, double p_height, double p_width);
    void create_nodes();
    void create_elements();
    void display_grid_data();
};

#endif // GRID_H