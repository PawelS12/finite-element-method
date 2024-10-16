#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <iostream>

class Integration {
private:
    static const double x2[2];
    static const double w2[2];
    static const double x3[3];
    static const double w3[3];

public:
    Integration();
    double gauss_integration_2D(double (*f)(double, double), int n, double a, double b, double c, double d); // funkcja | liczba punktów | granice całkowania
    void display_results(double result);
};

#endif // INTEGRATION_H