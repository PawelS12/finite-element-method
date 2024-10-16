#ifndef GLOBALDATA_H
#define GLOBALDATA_H

#include <iostream>
#include <fstream>

class GlobalData {
private:
    double simulation_time, simulation_step_time, conductivity, alfa, tot, initial_temp, density, specific_heat, nN, nE, nH, nW, H, W;

public:
    void read_file();
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
    void display_simulation_data();
};

#endif // GLOBALDATA_H