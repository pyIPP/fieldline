#ifndef include_fieldline_exceptions_undefinedMagneticField_hpp
#define include_fieldline_exceptions_undefinedMagneticField_hpp

#include <exception>

/*! Undefined magnetic field
    This exceptions is thrown if the magnetic field at the desired position is undefined.
    This can happen if a position outside the poloidal flux matrix is specified.
*/

namespace Fieldline {
    namespace exceptions {
        class undefinedMagneticField : public std::exception {
            public:
                /*! Error message */
                virtual const char * what() const throw() {
                    return "Undefined magnetic field encountered.";
                }
        };
    }
}

#endif
