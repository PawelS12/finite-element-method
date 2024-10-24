#ifndef GLOBALDATA_H
#define GLOBALDATA_H

#include <iostream>
#include "Node.h"
#include "Elem4.h"
#include <fstream>
#include <vector>

using std::vector;
using std::string;

class GlobalData {
private:
    double simulation_time, simulation_step_time, conductivity, alfa, tot, initial_temp, density, specific_heat, nN, nE, nH, nW, H, W;
    vector<Node> nodes_xy;
    vector<Elem4> elements;

public:
    void read_file();
    void read_node_coordinates();
    double get_simulation_time();
    double get_simulation_step_time();
    double get_conductivity();
    double get_alfa();
    double get_tot();
    double get_initial_temp();
    double get_density();
    double get_specific_heat();
    double get_nN();
    double get_nE();
    double get_nH();
    double get_nW();
    double get_height();
    double get_width();
    const std::vector<Elem4>& get_elements() const;
    void display_simulation_data();
    void create_elements_from_nodes();

};

#endif // GLOBALDATA_H