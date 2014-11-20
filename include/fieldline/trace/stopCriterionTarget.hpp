#ifndef include_fieldline_trace_stopCriterionTarget_hpp
#define include_fieldline_trace_stopCriterionTarget_hpp

#include <stdint.h>
#include <fieldline/core/stopCriterion.hpp>
#include <fieldline/core/target.hpp>

namespace Fieldline {
    namespace trace {
        class stopCriterionTarget : public Fieldline::core::stopCriterion {
            public:
                stopCriterionTarget(const Fieldline::core::target & target, const uint32_t N) : m_target(target), m_N(N) {}
                virtual bool reached(Fieldline::core::fieldline * fieldline)  {
                    if(fieldline->size() < 2) return false;
                    Fieldline::core::point point = m_target.intersection(fieldline->get_end());
                    if(point.hit) {
                        fieldline->replace_end(point.R, point.z);
                        return true;
                    }
                    return fieldline->size() >= m_N;
                }

            protected:
                Fieldline::core::target m_target;
                uint32_t m_N;
                
        };
    }
}


#endif
