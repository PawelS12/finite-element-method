#include "Node.h"

Node::Node() {}

Node::Node(double p_x, double p_y) : x(p_x), y(p_y) {}

double Node::get_x() { return x; }
double Node::get_y() { return y; }

void Node::display_node() {
    std::cout << "(" << x << ", " << y << ")" << std::endl;
}