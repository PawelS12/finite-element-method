#ifndef NODE_H
#define NODE_H

#include <iostream>

class Node {
private:
    double x, y;
    int BC;

public:
    Node();
    Node(double p_x, double p_y, int BC);
    double get_x() const;
    double get_y() const;
    int get_BC() const;
    void display_node() const;
};

#endif // NODE_H