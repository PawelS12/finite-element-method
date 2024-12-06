#include "Integration.h"
#include <iostream>
#include <math.h>
#include <functional>

using std::cout;
using std::endl;
using std::cerr;
using std::function;

const double Integration::x2[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
const double Integration::w2[2] = { 1.0, 1.0 };

const double Integration::x3[3] = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
const double Integration::w3[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

const double Integration::x4[4] = { -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0) };
const double Integration::w4[4] = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };

Integration::Integration() {}

void Integration::display_results(double result) {
    cout << "Result: " << result << endl << endl;
}

double Integration::gauss_integration_2D(function<double(double, double)> f, int n, double a, double b, double c, double d) {
    double suma = 0.0;
    const double* x; 
    const double* w;

    switch (n) {
        case 2:
            x = x2;
            w = w2;
            break;
        case 3:
            x = x3;
            w = w3;
            break;
        case 4:
            x = x4;
            w = w4;
            break;
        default:
            cerr << "Wrong number of points." << endl;
            return 0; 
    }

    double c1 = (b + a) / 2.0;
    double d1 = (b - a) / 2.0; 
    double c2 = (d + c) / 2.0; 
    double d2 = (d - c) / 2.0; 

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            suma += w[i] * w[j] * f(d1 * x[i] + c1, d2 * x[j] + c2);
        }
    }
    
    return d1 * d2 * suma; 
}

vector<double> Integration::get_weights_1D(int order)  {
    if (order == 2) {
        return  { 1.0, 1.0 };
    } else if (order == 4) {
        return { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
    } else if(order == 3) {
        return { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    } else {
        cerr << "Bad number of points!" << endl;
    }
    
    return {};
}

vector<double> Integration::get_points_1D(int order)  {
    if (order == 2) {
        return { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    } else if (order == 4) {
        return { -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0) };
    } else if(order == 3) {
        return { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
    } else {
        cerr << "Bad number of points!" << endl;
    }

    return {};
}