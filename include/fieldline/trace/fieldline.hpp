#ifndef include_fieldline_trace_fieldline_hpp
#define include_fieldline_trace_fieldline_hpp

#include <fieldline/core.hpp>
#include <fieldline/exceptions.hpp>
#include <fieldline/axiSymmetric/magneticField.hpp>
#include <fieldline/trace/rungeKutta.hpp>

namespace Fieldline {
    namespace trace {
        class fieldline : public Fieldline::core::fieldline {
            public:
                fieldline(Fieldline::axiSymmetric::magneticField * equilibrium, 
                          const double R, const double z, const double phi, const double dphi, 
                          Fieldline::core::stopCriterion * stopCriterion) : Fieldline::core::fieldline() {
                    Fieldline::trace::rungeKutta tracer(equilibrium);
                    tracer.init(R,z,phi,dphi);
                    Fieldline::core::fieldline::push_back(R,z,phi,tracer.get_Btot());
                    try {
                        while(!stopCriterion->reached(this)) {
                            tracer.next(1, dphi);
                            Fieldline::core::fieldline::push_back(tracer.get_R(), tracer.get_z(), tracer.get_phi(), tracer.get_Btot());
                        }
                    }
                    catch(Fieldline::exceptions::undefinedMagneticField & e) {
                    }
                }
        };
    }
}


#endif
