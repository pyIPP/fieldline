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

void testTarget() {
    Fieldline::core::target target("outer_target_curve.txt");
    target.save_to_file("newTarget.txt");
    Fieldline::axiSymmetric::magneticField equilibrium(-2.5, 1.65, "psi.txt");
    Fieldline::trace::stopCriterionTarget stopCriterion(target, 10000);
    Fieldline::trace::fieldline fieldline(&equilibrium, 2.1183, 0.0, 0.0, -0.01, &stopCriterion);
    fieldline.write_ASCII("fieldlineTestTarget.txt");
}

void testLine() {
    Fieldline::core::line line1(0,0,1,1);
    Fieldline::core::line line2(0,1,1,0);
    Fieldline::core::point point = line1.get_intersection(line2);
    std::cout << point << '\n';
}

void testSingularity() {
    Fieldline::axiSymmetric::magneticField magField(-2.5, 1.65, "psi.txt");
    Fieldline::core::point pointS = magField.get_singularity(1.6,-0.8);
    std::cout << "Singularity:" << '\t' <<"R:" << '\t' << pointS.R << '\t' <<"z:" << '\t'<< pointS.z << '\t' <<"exist:" << '\t' <<pointS.hit << '\n';
}


int main() {
    testMagneticField();
    testTarget();
    testLine();
    testSingularity();
    return EXIT_SUCCESS;
}

#endif
