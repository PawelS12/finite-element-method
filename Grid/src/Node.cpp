#include "Node.h"

using std::cout;
using std::endl;

Node::Node() : x(0), y(0), BC(0) {}

Node::Node(double p_x, double p_y, int p_BC) : x(p_x), y(p_y), BC(p_BC) {}

double Node::get_x() const { return x; }
double Node::get_y() const { return y; }
int Node::get_BC() const { return BC; }

void Node::display_node() const {
    cout << "Node: " << "(" << x << " ; " << y << "), status = " << (BC == 1 ? "contact" : "no contact") << endl;
}