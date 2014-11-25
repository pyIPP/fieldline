#ifndef include_fieldline_trace_stopCriterionTarget_hpp
#define include_fieldline_trace_stopCriterionTarget_hpp

#include <stdint.h>
#include <fieldline/core/stopCriterion.hpp>
#include <fieldline/core/target.hpp>


namespace Fieldline {
    namespace trace {
        /*! \brief This class implements the stop criterion for a target.
         *
         *  The stop criterion is reached if the specified target is hit.
         *  If the target is hit the last point on the field line is replaced by the hit point on target.
         *  If no target is reached after N points on the field line the stop criterion is also met.
         */
        class stopCriterionTarget : public Fieldline::core::stopCriterion {
            public:
                /*! \brief Constructor
                 *
                 *  The constructor requires a Fieldline::core::target instance representing the target.
                 *  The maximum number of points allowed on the field line is specified with N.
                 */
                stopCriterionTarget(const Fieldline::core::target & target, const uint32_t N) : m_target(target), m_N(N) {}
                /*! \brief Stop criterion
                 *
                 *  The stop criterion is reached if the last line segment of the field line intersects with the target.
                 *  The last point on the field line is replaced by the hit point on target.
                 *  If \f$N_{field line} \ge N\f$ the stop criterion is also met.
                 *  This is done to prevent an infinity loop in case the field line does not hit the target, e.g, confined region.
                */
                virtual bool reached(Fieldline::core::fieldline * fieldline)  {
                    if(fieldline->size() < 2) return false;
                    Fieldline::core::point point = m_target.get_intersection(fieldline->get_end());
                    if(point.hit) {
                        fieldline->replace_end(point.R, point.z);
                        return true;
                    }
                    return fieldline->size() >= m_N;
                }

            protected:
                Fieldline::core::target m_target; /*!< \brief Target */
                uint32_t m_N; /*!< \brief Maximum number of points on the field line */
                
        };
    }
}


#endif
