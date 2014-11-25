#ifndef source_fieldline_cpp
#define source_fieldline_cpp

#define WITH_PYTHON 1

#include <boost/python.hpp>
#include <fieldline.hpp>

class DummyCore{};

using namespace boost::python;
BOOST_PYTHON_MODULE(Fieldline) {
    // Fieldline.core
    {
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
        .add_property("turns", &Fieldline::core::fieldline::get_number_of_turns)
        .def("write_to_file", &Fieldline::core::fieldline::write_to_file)
        .def("read_from_file", &Fieldline::core::fieldline::read_from_file)
        .def("replace_end", &Fieldline::core::fieldline::replace_end)
        ;

        class_<Fieldline::core::stopCriterion, boost::noncopyable>("stopCriterion", no_init)
            .def("reached", pure_virtual(&Fieldline::core::stopCriterion::reached))
            ;
        register_ptr_to_python<Fieldline::core::stopCriterion *>();

        class_<Fieldline::core::rhoPoloidal>("rhoPoloidal", init<double, double>())
            .def(init<Fieldline::core::rhoPoloidal>())
            .def("getRhoPoloidal", &Fieldline::core::rhoPoloidal::get_rho_poloidal)
            .def("getPoloidalFlux", &Fieldline::core::rhoPoloidal::get_poloidal_flux)
            ;
    }

    // Fieldline.exceptions
    {
        object exceptionsModule(handle<>(borrowed(PyImport_AddModule("Fieldline.exceptions"))));
        scope().attr("exceptions") = exceptionsModule;
        scope exceptions_scope = exceptionsModule;

        class_<Fieldline::exceptions::undefinedMagneticField>("undefinedMagneticField")
            .def("__str__", &Fieldline::exceptions::undefinedMagneticField::what)
            .def("what", &Fieldline::exceptions::undefinedMagneticField::what)
            ;

        class_<Fieldline::exceptions::traceNotInitialized>("traceNotInitialized")
            .def("__str__", &Fieldline::exceptions::traceNotInitialized::what)
            .def("what", &Fieldline::exceptions::traceNotInitialized::what)
            ;

    
    }

    // Fieldline.axiSymmetric
    {
        object axiSymmetricModule(handle<>(borrowed(PyImport_AddModule("Fieldline.axiSymmetric"))));
        scope().attr("axiSymmetric") = axiSymmetricModule;
        scope axiSymmetric_scope = axiSymmetricModule;

        class_<Fieldline::axiSymmetric::magneticField>("magneticField")
        .def(init<double, double, std::string>())
        .def("get_magnetic_flux", &Fieldline::axiSymmetric::magneticField::get_magnetic_flux)
        .def("get_magnetic_field", &Fieldline::axiSymmetric::magneticField::get_magnetic_field)
        .def("write_ASCII_matrix", &Fieldline::axiSymmetric::magneticField::write_ASCII_matrix)
        .def("get_singularity", &Fieldline::axiSymmetric::magneticField::get_singularity)
        ;
        register_ptr_to_python<Fieldline::axiSymmetric::magneticField *>();
    }

    // Fieldline.trace
    {
        object traceModule(handle<>(borrowed(PyImport_AddModule("Fieldline.trace"))));
        scope().attr("trace") = traceModule;
        scope trace_scope = traceModule;

        class_<Fieldline::trace::fieldline, bases<Fieldline::core::fieldline> >("fieldline", init<Fieldline::axiSymmetric::magneticField *, const double, const double, const double, const double, Fieldline::core::stopCriterion* >())
            ;

        class_<Fieldline::trace::stopCriterionSteps, bases<Fieldline::core::stopCriterion> >("stopCriterionSteps", init<uint32_t>())
            .def("reached", &Fieldline::trace::stopCriterionSteps::reached)
            ;

        class_<Fieldline::trace::stopCriterionTarget, bases<Fieldline::core::stopCriterion> >("stopCriterionTarget", init<Fieldline::core::target, uint32_t>())
            .def("reached", &Fieldline::trace::stopCriterionTarget::reached)
            ;

        class_<Fieldline::trace::rungeKutta>("rungeKutta", init<Fieldline::axiSymmetric::magneticField *>())
            .def("init", &Fieldline::trace::rungeKutta::init)
            .def("next", &Fieldline::trace::rungeKutta::next)
            .add_property("R", &Fieldline::trace::rungeKutta::get_R)
            .add_property("z", &Fieldline::trace::rungeKutta::get_z)
            .add_property("phi", &Fieldline::trace::rungeKutta::get_phi)
            .add_property("Btot", &Fieldline::trace::rungeKutta::get_Btot)
            ;

        
    }

}

#endif 
