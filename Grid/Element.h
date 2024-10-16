#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>

class Element {
private:
    int ID[4];

public:
    Element();
    void set_ID(int node_1, int node_2, int node_3, int node_4);
    void display_ID();
};

#endif // ELEMENT_H