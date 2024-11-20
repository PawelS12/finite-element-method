#ifndef NODE_H
#define NODE_H

#include <iostream>

class Node {
private:
    double x, y;

public:
    Node();
    Node(double p_x, double p_y);
    double get_x() const;
    double get_y() const;
    void display_node() const;
};

#endif // NODE_H