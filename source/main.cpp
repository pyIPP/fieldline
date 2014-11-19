#ifndef main_cpp
#define main_cpp

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fieldline.hpp>

void testMagneticField() {
    fieldline::core::magneticField magneticField;
    magneticField.BR = 0.1;
    magneticField.Bz = -0.1;
    magneticField.Btor = -2.5;
    std::cout << magneticField;
    std::cout << magneticField.Btot() << '\n';
    std::cout << magneticField.Bpol() << '\n';
    fieldline::axiSymmetric::magneticField equilibrium(-2.5, 1.65, "psi.txt");
    std::cout << equilibrium.get_magnetic_field(1.0, 0.0);
}

int main() {
    std::cout << "Fieldline\n";
    testMagneticField();
    return EXIT_SUCCESS;
}

#endif
