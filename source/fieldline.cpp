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

    class_<Fieldline::core::fieldline>("fieldline")
        .def(init<Fieldline::core::fieldline>())
        .def(init<boost::python::list, boost::python::list, boost::python::list, boost::python::list>())
        .add_property("size", &Fieldline::core::fieldline::size)
        .def("push_back", &Fieldline::core::fieldline::push_back)
        .add_property("length", &Fieldline::core::fieldline::get_length)
        .add_property("lastR", &Fieldline::core::fieldline::get_last_R)
        .add_property("lastz", &Fieldline::core::fieldline::get_last_z)
        .add_property("end", &Fieldline::core::fieldline::get_end)
        .add_property("x", &Fieldline::core::fieldline::get_x_python)
        .add_property("y", &Fieldline::core::fieldline::get_y_python)
        .add_property("z", &Fieldline::core::fieldline::get_z_python)
        .add_property("R", &Fieldline::core::fieldline::get_R_python)
        .add_property("phi", &Fieldline::core::fieldline::get_phi_python)
        .add_property("Btot", &Fieldline::core::fieldline::get_Btot_python)
        .def("write_to_file", &Fieldline::core::fieldline::write_to_file)
        .def("read_from_file", &Fieldline::core::fieldline::read_from_file)
        .def("replace_end", &Fieldline::core::fieldline::replace_end)
        ;

    object axiSymmetricModule(handle<>(borrowed(PyImport_AddModule("Fieldline.axiSymmetric"))));
    scope().attr("axiSymmetric") = axiSymmetricModule;
    scope axiSymmetric_scope = axiSymmetricModule;

    class_<Fieldline::axiSymmetric::magneticField>("magneticField")
        ;

}

#endif 
