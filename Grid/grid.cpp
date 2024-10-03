#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::ifstream;


class Node {
    private:
        int x, y;
};

class Element {
    private:
        int ID[4];
};

class Grid {
    private:
        int nN, nE;
        // element[nE]
        // node [nN]
};

class GlobalData {
    private:
        double simulation_time, simulation_step_time, conductivity, alfa, tot, initial_temp, density, specific_heat, nN, nE;
    public:
        void read_file() {
            ifstream file("data.txt");
            if (!file.is_open()) {
                cout << "File not found." << endl;
            }
            
            file >> simulation_time;
            file >> simulation_step_time;
            file >> conductivity;
            file >> alfa;
            file >> tot;
            file >> initial_temp;
            file >> density;
            file >> specific_heat;
            file >> nN;
            file >> nE;
            
            file.close();

            cout << "Simulation time: " << simulation_time << endl;
            cout << "Simulation step time: " << simulation_step_time << endl;
            cout << "Conductivity: " << conductivity << endl;
            cout << "Alfa: " << alfa << endl;
            cout << "Tot: " << tot << endl;
            cout << "Initial temperature: " << initial_temp << endl;
            cout << "Density: " << density << endl;
            cout << "Specific heat: " << specific_heat << endl;
            cout << "Number of nodes: " << nN << endl;
            cout << "Number of elements: " << nE << endl;
        }
};

int main() {
    GlobalData data;
    data.read_file();


    return 0;
}