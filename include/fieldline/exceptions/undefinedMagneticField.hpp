#ifndef include_fieldline_exceptions_undefinedMagneticField_hpp
#define include_fieldline_exceptions_undefinedMagneticField_hpp

#include <exception>


namespace Fieldline {
    namespace exceptions {
        /*! \brief Undefined magnetic field
         *
         *  This exceptions is thrown if the magnetic field at the desired position is undefined.
         *  This can happen if a position outside the poloidal flux matrix is specified.
         */
        class undefinedMagneticField : public std::exception {
            public:
                /*! \brief Error message */
                virtual const char * what() const throw() {
                    return "Undefined magnetic field encountered.";
                }
        };
    }
}

#endif
