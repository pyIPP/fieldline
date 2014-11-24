#ifndef include_fieldline_core_stopCriterion_hpp
#define include_fieldline_core_stopCriterion_hpp

#include <fieldline/core/fieldline.hpp>

/*! Stop Criterion
    Abstract class for the definition of stop criteria for the field line trace.
    This class should only be used as a parent for the designated stop criteria.
    Never use an instance of this class, only use its children.

*/

namespace Fieldline {
    namespace core {
        class stopCriterion {
            public:
                /*! Stop criterion reached.
                    This function is defined as an abstract function.
                    In the derived classes this function has to be overridden with the stop criteria.
                    If the stop criteria is reached the function has to return true.
                    If the criteria is not met the field line trace has to continue and the function has to return false.
                */
                virtual bool reached(fieldline * rhs) = 0;
        };
    }
}

#endif
