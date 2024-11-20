#ifndef FEMSOLVER_H
#define FEMSOLVER_H 

#include <iostream>
#include <vector>
#include "Grid.h"   
#include "Element.h" 

class FEMSolver {
private:
    Grid& grid;
    vector<vector<double>> local_H_matrices;
public:
    explicit FEMSolver(Grid& grid);
    void compute_jacobian(const Element& element, double xi, double eta, double J[2][2]) const;
    double compute_jacobian_determinant(double J[2][2]) const; 
    void compute_inverse_jacobian(double J[2][2], double invJ[2][2]) const; 
    void calculate_H_matrix(double conductivity);
    double calculate_H_integrand(const Element& element, double conductivity, int i, int j, double xi, double eta) const;
    void aggregate_H_matrix(vector<vector<double>>& H_global, int nodes_num) const;
};

#endif // FEMSOLVER_H
