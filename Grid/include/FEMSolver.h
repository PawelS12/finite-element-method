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
    explicit FEMSolver(Grid& grid, double alpha);
    void display_matrix(vector<vector<double>>& matrix);
    void compute_jacobian(const Element& element, double xi, double eta, double J[2][2]) const;
    double compute_jacobian_determinant(double J[2][2]) const; 
    void compute_inverse_jacobian(double J[2][2], double invJ[2][2]) const; 
    void calculate_H_matrix(double conductivity);
    double calculate_H_integrand(const Element& element, double conductivity, int i, int j, double xi, double eta) const;
    void aggregate_H_matrix(vector<vector<double>>& H_global, int nodes_num) const;
    void calculate_Hbc_matrix(double alpha);
    void integrate_Hbc_on_edge(const Node& node1, const Node& node2, double alpha, vector<vector<double>>& Hbc) const;
    void compute_edge_jacobian(const Node& node1, const Node& node2, double xi, double& detJ) const;
};

#endif // FEMSOLVER_H
