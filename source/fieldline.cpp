#ifndef source_fieldline_cpp
#define source_fieldline_cpp

#define WITH_PYTHON 1

#include <boost/python.hpp>
#include <fieldline.hpp>

class DummyCore{};

using namespace boost::python;
BOOST_PYTHON_MODULE(Fieldline) {
    object coreModule(handle<>(borrowed(PyImport_AddModule("Fieldline.core"))));
    scope().attr("core") = coreModule;
    scope core_scope = coreModule;

    class_<Fieldline::core::point>("point")
        .def(init<double, double, double>())
        .def(init<Fieldline::core::point>())
        .def("get_distance", &Fieldline::core::point::get_distance)
        .add_property("R", &Fieldline::core::point::R, &Fieldline::core::point::R)
        .add_property("z", &Fieldline::core::point::z, &Fieldline::core::point::z)
        .add_property("hit", &Fieldline::core::point::hit, &Fieldline::core::point::hit)
        ;

    class_<Fieldline::core::line>("line", init<double, double, double, double>())
        .def(init<Fieldline::core::point, Fieldline::core::point>())
        .def(init<Fieldline::core::line>())
        .def("get_intersection", &Fieldline::core::line::get_intersection)
        .add_property("R0", &Fieldline::core::line::get_R0)
        .add_property("z0", &Fieldline::core::line::get_z0)
        .add_property("R1", &Fieldline::core::line::get_R1)
        .add_property("z1", &Fieldline::core::line::get_z1)
        .add_property("length", &Fieldline::core::line::get_length)
        ;

    class_<Fieldline::core::magneticFlux>("magneticFlux")
        .def(init<double, double, double>())
        .def(init<Fieldline::core::magneticFlux>())
        .add_property("psi", &Fieldline::core::magneticFlux::psi, &Fieldline::core::magneticFlux::psi)
        .add_property("dpsi_dR", &Fieldline::core::magneticFlux::dpsi_dR, &Fieldline::core::magneticFlux::dpsi_dR)
        .add_property("dpsi_dz", &Fieldline::core::magneticFlux::dpsi_dz, &Fieldline::core::magneticFlux::dpsi_dz)
        ;

    class_<Fieldline::core::magneticField>("magneticField")
        .def(init<double, double, double>())
        .def(init<Fieldline::core::magneticField>())
        .add_property("BR", &Fieldline::core::magneticField::BR, &Fieldline::core::magneticField::BR)
        .add_property("Bz", &Fieldline::core::magneticField::Bz, &Fieldline::core::magneticField::Bz)
        .add_property("Btor", &Fieldline::core::magneticField::Btor, &Fieldline::core::magneticField::Btor)
        .add_property("Bpol", &Fieldline::core::magneticField::Bpol)
        .add_property("Btot", &Fieldline::core::magneticField::Btot)
        ;

    class_<Fieldline::core::target>("target", init<std::string>())
        .def(init<Fieldline::core::target>())
        .def(init<boost::python::list, boost::python::list>())
        .def("get_intersection", &Fieldline::core::target::get_intersection)
        .def("save_to_file", &Fieldline::core::target::save_to_file)
        .def("load_from_file", &Fieldline::core::target::load_from_file)
        .add_property("R", &Fieldline::core::target::get_R_python)
        .add_property("z", &Fieldline::core::target::get_z_python)
        ;


}

#endif 
