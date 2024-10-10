#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

class Node { // węzły siatki oraz ich współrzędne
    private:
        double x, y;
    public:
        Node() {}
        Node(double p_x, double p_y) : x(p_x), y(p_y) {}

        double get_x() { return x; }
        double get_y() { return y; }

        void display_node() {
            cout << "(" << x << ", " << y  << ")" << endl;
        }
};

class Element { // elementy tworzone przez 4 węzły (obiekty Node)
    private:
        int ID[4];
    public:
        Element() {
            for (int i = 0; i < 4; i++) {
                ID[i] = 0;
            }
        }

        void set_ID(int node_1, int node_2, int node_3, int node_4) {
            ID[0] = node_1;
            ID[1] = node_2;
            ID[2] = node_3;
            ID[3] = node_4;
        }

        void display_ID() {
            for (int i = 0; i < sizeof(ID) / sizeof(ID[0]); i++) {
                cout << ID[i] << " ";
            }
            cout << endl;
        }
};

class Grid { // siatka tworzona z elementów i węzłów
    private:
        double nN, nE, nW, nH, height, width;
        vector<Node> nodes;
        vector<Element> elements;

    public:
        Grid() {}
        Grid(double p_nN, double p_nE, double p_nW, double p_nH, double p_height, double p_width) : nN(p_nN), nE(p_nE), nW(p_nW), nH(p_nH), height(p_height), width(p_width) {
            create_nodes();
            create_elements();
        }

        void create_nodes() {
            double delta_x = width / (nW - 1);   // odstępy pomiędzy kolejnymi węzłami
            double delta_y = height / (nH - 1); 

            for (int i = 0; i < nH; ++i) {
                for (int j = 0; j < nW; ++j) { 
                    double x = i * delta_x;  
                    double y = j * delta_y; 
                    nodes.push_back(Node(x, y)); 
                }
            }
        }

        void create_elements() {
            for (int i = 0; i < nH - 1; ++i) {  
                for (int j = 0; j < nW - 1; ++j) {  
                    int node_1 = i * nW + j + 1;    
                    int node_2 = node_1 + 1;             
                    int node_3 = node_1 + nW + 1;       
                    int node_4 = node_1 + nW;          

                    Element element;                  
                    element.set_ID(node_1, node_2, node_3, node_4); 
                    elements.push_back(element);       
                }
            }
        }

        void display_grid_data() {
            //cout << "Grid data: nN = " << nN << ", nE = " << nE << ", nH = " << nH << ", nW = " << nH <<  ", H = " << height << ", W = " << width << endl << endl; 

            cout << "Nodes: " << endl << endl;
            for (size_t i = 0; i < nodes.size(); i++) {
                cout << "Node " << i + 1 << ": ";
                nodes[i].display_node();
            }
            cout << "-----------------------------------" << endl;
            cout << "Elements: " << endl << endl;
            for (size_t i = 0; i < elements.size(); i++) {
                cout << "Element " << i + 1 << ": ";
                elements[i].display_ID();
            }
            cout << "-----------------------------------" << endl;
        }
};

class GlobalData { // dane do symulacji
    private:
        double simulation_time, simulation_step_time, conductivity, alfa, tot, initial_temp, density, specific_heat, nN, nE, nH, nW, H, W;
    public:
        void read_file() {
            ifstream file("data.txt");
            if (!file.is_open()) {
                cout << "File not found." << endl;
            }
            
            file >> simulation_time >> simulation_step_time >> conductivity >> alfa >> tot >> initial_temp >> density >> specific_heat >> nN >> nE >> nH >> nW >> H >> W;
            
            file.close();
        }

        double get_simulation_time() { return simulation_time; }
        double get_simulation_step_time() { return simulation_step_time; }
        double get_conductivity() { return conductivity; }
        double get_alfa() { return alfa; }
        double get_tot() { return tot; }
        double get_initial_temp() { return initial_temp; }
        double get_density() { return density; }
        double get_specific_heat() { return specific_heat; }
        double get_nN() { return nN; }
        double get_nE() { return nE; }
        double get_nH() { return nH; }
        double get_nW() { return nW; }
        double get_height() { return H; }
        double get_width() { return W; }

        void display_simulation_data() {
            cout << "-----------------------------------" << endl;
            cout << "Simulation data: " << endl << endl;
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
            cout << "Nodes height: " << nH << endl;
            cout << "Nodes width: " << nW << endl;
            cout << "Height: " << H << endl;
            cout << "Width: " << W << endl;
            cout << "-----------------------------------" << endl;
        }
};

int main() {
    GlobalData data;
    data.read_file();
    Grid grid_1(data.get_nN(), data.get_nE(), data.get_nW(), data. get_nH(), data.get_height(), data.get_width());
    data.display_simulation_data();
    grid_1.display_grid_data();

    return 0;
}