#ifndef include_fieldline_core_stopCriterion_hpp
#define include_fieldline_core_stopCriterion_hpp

#include <fieldline/core/fieldline.hpp>

namespace Fieldline {
    namespace core {
        class stopCriterion {
            public:
                virtual bool reached(const fieldline * rhs) = 0;
        };
    }
}

#endif
