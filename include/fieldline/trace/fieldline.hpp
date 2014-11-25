#ifndef include_fieldline_trace_fieldline_hpp
#define include_fieldline_trace_fieldline_hpp

#include <fieldline/core.hpp>
#include <fieldline/exceptions.hpp>
#include <fieldline/axiSymmetric/magneticField.hpp>
#include <fieldline/trace/rungeKutta.hpp>

#define FIELDLINE_FORWARDS          1
#define FIELDLINE_BACKWARDS         2
#define FIELDLINE_BIDIRECTIONAL     3


namespace Fieldline {
    namespace trace {
        /*! \brief Field line trace on a magnetic field provided a starting point.
         *
         *  This class is derived from Fieldline::core::fieldline.
         *  Given a magnetic field, a starting point \f$(R,z,\phi)\f$, a step width \f$\mathrm{d}\phi\f$ and a stop criterion,
         *  the field line trace will be performed and the result will be stored in the resulting fieldline instance.
         */
        class fieldline : public Fieldline::core::fieldline {
            public:
                /*! \brief Constructor
                 *
                 *  The constructor requires the pointer to a Fieldline::axiSymmetric::magneticField instance or an instance derived from it.
                 *  In addition a starting point \f$(R,z,\phi)\f$ is required where to start the field line trace.
                 *  \f$\mathrm{d}\phi\f$ is specified to give the step width
                 *  A pointer to an instance derived from Fieldline::core::stopCriterion is required. 
                 *  The stopCriterion can be any arbitrary criterion, e.g. field line length, target intersection ...
                 *  \param equilibrium Pointer to the magnetic field instance.
                 *  \param R Start position \f$R\f$ of the field line trace.
                 *  \param z Start position \f$z\f$ of the field line trace.
                 *  \param phi Toroidal angle \f$\phi\f$ where the field line trace starts.
                 *  \param dphi Toroidal step with of the field line trace \f$\mathrm{d}\phi\f$.
                 *  \param stopCriterion Pointer to the stop criterion of the field line trace.
                 *  \param direction Specified the direction the field line will be traced in. 1 forwards, 2 backwards, 3 both directions.
                 */
                fieldline(Fieldline::axiSymmetric::magneticField * equilibrium,
                          const double R,
                          const double z,
                          const double phi,
                          const double dphi,
                          Fieldline::core::stopCriterion * stopCriterion,
                          const uint32_t direction
                          ) : Fieldline::core::fieldline() {
                    switch(direction)
                    {
                        case FIELDLINE_FORWARDS:
                        {
                            Fieldline::trace::rungeKutta tracer(equilibrium);
                            tracer.init(R,z,phi,std::abs(dphi));
                            Fieldline::core::fieldline::push_back(R,z,phi,tracer.get_Btot());
                            try {
                                while(!stopCriterion->reached(this)) {
                                    tracer.next(1, std::abs(dphi));
                                    Fieldline::core::fieldline::push_back(tracer.get_R(), tracer.get_z(), tracer.get_phi(), tracer.get_Btot());
                                }
                            }
                            catch(Fieldline::exceptions::undefinedMagneticField & e) {
                            }
                            break;
                        }
                        case FIELDLINE_BACKWARDS:
                        {
                            Fieldline::trace::rungeKutta tracer(equilibrium);
                            tracer.init(R,z,phi,-std::abs(dphi));
                            Fieldline::core::fieldline::push_back(R,z,phi,tracer.get_Btot());
                            try {
                                while(!stopCriterion->reached(this)) {
                                    tracer.next(1, -std::abs(dphi));
                                    Fieldline::core::fieldline::push_back(tracer.get_R(), tracer.get_z(), tracer.get_phi(), tracer.get_Btot());
                                }
                            }
                            catch(Fieldline::exceptions::undefinedMagneticField & e) {
                            }
                            break;
                        }
                        case FIELDLINE_BIDIRECTIONAL:
                        {
                            Fieldline::trace::rungeKutta tracer(equilibrium);
                            tracer.init(R,z,phi,std::abs(dphi));
                            Fieldline::core::fieldline::push_back(R,z,phi,tracer.get_Btot());
                            try {
                                while(!stopCriterion->reached(this)) {
                                    tracer.next(1, std::abs(dphi));
                                    Fieldline::core::fieldline::push_back(tracer.get_R(), tracer.get_z(), tracer.get_phi(), tracer.get_Btot());
                                }
                            }
                            catch(Fieldline::exceptions::undefinedMagneticField & e) {
                            }
                            Fieldline::core::fieldline::reverse();
                            tracer.init(R,z,phi,-std::abs(dphi));
                            try {
                                while(!stopCriterion->reached(this)) {
                                    tracer.next(1, -std::abs(dphi));
                                    Fieldline::core::fieldline::push_back(tracer.get_R(), tracer.get_z(), tracer.get_phi(), tracer.get_Btot());
                                }
                            }
                            catch(Fieldline::exceptions::undefinedMagneticField & e) {
                            }
                            break;
                        }
                    }
                }
        };
    }
}


#endif
