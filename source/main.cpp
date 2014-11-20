#ifndef main_cpp
#define main_cpp

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fieldline.hpp>

void testMagneticField() {
    Fieldline::core::magneticField magneticField;
    magneticField.BR = 0.1;
    magneticField.Bz = -0.1;
    magneticField.Btor = -2.5;
    Fieldline::axiSymmetric::magneticField equilibrium(-2.5, 1.65, "psi.txt");
    equilibrium.write_ASCII_matrix("dumpMatrix.txt");
    std::fstream file("output.txt", std::ios::out);
    for(uint32_t i = 0; i < 10000; ++i) {
        magneticField = equilibrium.get_magnetic_field(1.1+1.2/10000.0*(double)i, 0.0, 0.0);
        file << 1.1+1.2/10000.0*(double)i << '\t' << magneticField.BR << '\t' << magneticField.Bz << '\t' << magneticField.Btor << '\n';
    }
    file.close();
    Fieldline::trace::stopCriterionSteps stopCriterion(10000);
    Fieldline::trace::fieldline fieldline(&equilibrium, 2.1183, 0.0, 0.0, 0.01, &stopCriterion);
    fieldline.write_ASCII("fieldlineTest.txt");
    Fieldline::trace::fieldline fieldline2(&equilibrium, 2.1183, 0.0, 0.0, -0.01, &stopCriterion);
    fieldline2.write_ASCII("fieldlineTest2.txt");

    Fieldline::trace::rungeKutta tracer(&equilibrium);
    tracer.init(2.0, 0.0, 0.0, -0.01);
    file.open("trace.txt", std::ios::out);
    file << tracer.get_R() << '\t' << tracer.get_z() << '\t' << tracer.get_phi() << '\t' << tracer.get_Btot() << '\n';
    for(uint32_t i = 0; i < 10000; ++i) {
        tracer.next(1, -0.01);
        file << tracer.get_R() << '\t' << tracer.get_z() << '\t' << tracer.get_phi() << '\t' << tracer.get_Btot() << '\n';
    }
    file.close();
}

int main() {
    testMagneticField();
    return EXIT_SUCCESS;
}

#endif
