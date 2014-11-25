#ifndef include_fieldline_trace_stopCriterionSteps_hpp
#define include_fieldline_trace_stopCriterionSteps_hpp

#include <fieldline/core/stopCriterion.hpp>

/*! Stop criterion after N steps.
    This is a stop criterion derived from Fieldline::core::stopCriterion.
    The stop criterion is meet if the field line has N or more points.
*/

namespace Fieldline {
    namespace trace {
        class stopCriterionSteps : public Fieldline::core::stopCriterion {
            public:
                /*! Constructor
                    The constructor requires the number of points after which the stop criterion is met.
                */
                stopCriterionSteps(const uint32_t N) : m_N(N){}
                /*! Stop criterion is reached.
                    The stop criterion is reached if \f$N_{field line} \ge N\f$.
                */
                virtual bool reached(Fieldline::core::fieldline * rhs) {
                    return rhs->size() >= m_N;
                }
            protected:
                /*! Number of points until stop criterion is reached. */
                uint32_t m_N;

        };
    }
}

#endif
