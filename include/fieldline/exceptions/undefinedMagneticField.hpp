#ifndef include_fieldline_exceptions_undefinedMagneticField_hpp
#define include_fieldline_exceptions_undefinedMagneticField_hpp

#include <exception>

namespace Fieldline {
    namespace exceptions {
        class undefinedMagneticField : public std::exception {
            public:
                virtual const char * what() const throw() {
                    return "Undefined magnetic field encountered.";
                }
        };
    }
}

#endif
