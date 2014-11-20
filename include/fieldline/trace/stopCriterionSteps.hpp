#ifndef include_fieldline_trace_stopCriterionSteps_hpp
#define include_fieldline_trace_stopCriterionSteps_hpp

#include <fieldline/core/stopCriterion.hpp>

namespace Fieldline {
    namespace trace {
        class stopCriterionSteps : public Fieldline::core::stopCriterion {
            public:
                stopCriterionSteps(const uint32_t N) : m_N(N){}
                virtual bool reached(const Fieldline::core::fieldline * rhs) {
                    return rhs->size() >= m_N;
                }
            protected:
                uint32_t m_N;

        };
    }
}

#endif
