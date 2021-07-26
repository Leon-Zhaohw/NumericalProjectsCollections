#include <vector>
#include <iostream>
#include "ceres/jet.h"

// Completely arbitrary function that takes two inputs
// and returns two outputs
template<typename T>
std::vector<T> simpleFunc(const T &x,
                          const T &y)
{
    std::vector<T> output(2);

    output[0] = x * y;            // f_1
    output[1] = sin(2.0 * x + y); // f_2

    return output;
}

// For this simple function, we can compute the
// analytical Jacobian
std::vector<double> jacobian(double x,
                             double y)
{
    std::vector<double> output(4);

    output[0] = y;                       // df_1 / dx
    output[1] = x;                       // df_1 / dy
    output[2] = std::cos(2.0*x+y) * 2.0; // df_2 / dx
    output[3] = std::cos(2.0*x+y);       // df_2 / dy

    return output;
}

// The Jet type will keep track of derivatives for us
// The 2 tells it to keep track of two derivatives at once
// The Jet class has two member variables: "a" holds the variable value,
// and "v" is the vector of all the derivatives it is keeping track of.
typedef ceres::Jet<double, 2> MyJet;

int main(int argc, char** argv)
{
    double x = 5;
    double y = 3;

    // Can call this function and the analytical jacobian function
    std::vector<double> output = simpleFunc(x, y);
    std::vector<double> jac = jacobian(x, y);

    // Can also call it with jet variables
    // The x derivative is in the 0th position,
    // and the y derivative is in the 1st position
    MyJet xJet(x, 0);
    MyJet yJet(y, 1);

    std::vector<MyJet> jetOutput = simpleFunc(xJet, yJet);

    printf("x: %f, y: %f\n", x, y);

    // The a member variable holds the value
    printf("f_1 = x*y = %f\nf_2 = sin(2x+y) = %f\n\n", jetOutput[0].a, jetOutput[1].a);

    // The v member variable holds the derivatives
    printf("          Analytical Automatic\n");
    printf("df_1/dx: %9f  %9f\n", jac[0], jetOutput[0].v(0));
    printf("df_1/dy: %9f  %9f\n", jac[1], jetOutput[0].v(1));
    printf("df_2/dx: %9f  %9f\n", jac[2], jetOutput[1].v(0));
    printf("df_2/dy: %9f  %9f\n", jac[3], jetOutput[1].v(1));

    return 0;
}

