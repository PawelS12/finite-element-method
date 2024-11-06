#ifndef ELEM4_H 
#define ELEM4_H 

#include <iostream>
#include <vector>
#include "Node.h"

class Elem4 {
private:
    Node nodes_xy[4]; 

public:
    Elem4(const Node& n1, const Node& n2, const Node& n3, const Node& n4);
    void compute_jacobian(double xi, double eta, double J[2][2]) const;
    double compute_jacobian_determinant(double J[2][2]) const; 
    void compute_inverse_jacobian(double J[2][2], double invJ[2][2]) const; 
    void calculate_H_matrix(double conductivity) const;
    double calculate_H_integrand(double conductivity, int i, int j, double xi, double eta) const;
};

#endif 
