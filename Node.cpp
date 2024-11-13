#include "Node.h"

using std::cout;
using std::endl;

Node::Node() : x(0), y(0) {}

Node::Node(double p_x, double p_y) : x(p_x), y(p_y) {}

double Node::get_x() const { return x; }
double Node::get_y() const { return y; }

void Node::display_node() const {
    cout << "(" << x << ", " << y << ")" << endl;
}