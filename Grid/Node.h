#ifndef NODE_H
#define NODE_H

#include <iostream>

class Node {
private:
    double x, y;

public:
    Node();
    Node(double p_x, double p_y);
    double get_x();
    double get_y();
    void display_node();
};

#endif // NODE_H