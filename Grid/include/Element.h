#ifndef ELEMENT_H 
#define ELEMENT_H 

#include <iostream>
#include <vector>
#include "Node.h"

using std::vector;

class Element {
private:
    Node nodes_xy[4]; 
    vector<int> ID;
    vector<vector<double>> Hbc;

public:
    Element(const Node& n1, const Node& n2, const Node& n3, const Node& n4);
    void set_ID(int node_1, int node_2, int node_3, int node_4);
    void display_ID() const;
    const Node* get_nodes() const;
    const std::vector<int>& get_ID() const;
    void set_Hbc(const vector<vector<double>>& hbc);
    const vector<vector<double>>& get_Hbc() const;
};

#endif // ELEMENT_H