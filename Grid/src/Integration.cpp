#include "Integration.h"
#include <iostream>
#include <math.h>
#include <functional>

using std::cout;
using std::endl;
using std::cerr;
using std::function;

const double Integration::x2[] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
const double Integration::w2[] = { 1.0, 1.0 };

const double Integration::x3[] = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
const double Integration::w3[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

const double Integration::x4[] = { -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0) };
const double Integration::w4[] = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };

const double Integration::x9[9] = { -0.968160, -0.836031, -0.613371, -0.324253, 0.0, 0.324253, 0.613371, 0.836031, 0.968160 };
const double Integration::w9[9] = { 0.081274, 0.180648, 0.260611, 0.312347, 0.330239, 0.312347, 0.260611, 0.180648, 0.081274 };

const double Integration::x16[16] = {-0.989400934991649932596154173450,-0.944575023073232576077988415535, -0.865631202387831743880467897712, -0.755404408355003033895101194847, -0.617876244402643748446671764049, -0.458016777657227386342419442984, -0.281603550779258913230460501460, -0.095012509837637440185319335425, 0.095012509837637440185319335425,  0.281603550779258913230460501460, 0.458016777657227386342419442984,  0.617876244402643748446671764049, 0.755404408355003033895101194847,  0.865631202387831743880467897712, 0.944575023073232576077988415535,  0.989400934991649932596154173450 };

const double Integration::w16[16] = { 0.027152459411754094851780572456, 0.062253523938647892862843836994, 0.095158511682492784809925107602, 0.124628971255533872052476282192, 0.149595988816576732081501730547, 0.169156519395002538189312079030, 0.182603415044923588866763667969, 0.189450610455068496285396723208, 0.189450610455068496285396723208,0.182603415044923588866763667969,0.169156519395002538189312079030, 0.149595988816576732081501730547,  0.124628971255533872052476282192,  0.095158511682492784809925107602,  0.062253523938647892862843836994, 0.027152459411754094851780572456 };

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
        case 9:
            x = x9;
            w = w9;
            break;
        case 16:
            x = x16;
            w = w16;
            break;
        default:
            cerr << "Error in gauss_integration_2D: Unsupported number of points (n = " << n << "). Supported values: 2, 3, 4, 9, 16." << endl;
            throw std::invalid_argument("Unsupported number of points for Gauss integration");
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

vector<double> Integration::get_weights(int order)  {
    if (order == 2) {
        return  { 1.0, 1.0 };
    } else if (order == 4) {
        return { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
    } else if(order == 3) {
        return { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    } else {
        cerr << "Wrong number of points!" << endl;
        return {};
    }
}

vector<double> Integration::get_points(int order)  {
    if (order == 2) {
        return { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    } else if (order == 4) {
        return { -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0), -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0), sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0) };
    } else if(order == 3) {
        return { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
    } else {
        cerr << "Wrong number of points!" << endl;
        return {};
    }
}