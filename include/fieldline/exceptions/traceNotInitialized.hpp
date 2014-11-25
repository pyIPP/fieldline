#ifndef include_fieldline_exceptions_traceNotInitialized_hpp
#define include_fieldline_exceptions_traceNotInitialized_hpp

#include <exception>

/*! Exception if the field line trace is not initialized.
    This exception is thrown in the Runge Kutta field line trace if no start point is defined.
*/

namespace Fieldline {
    namespace exceptions {
        class traceNotInitialized : public std::exception {
            public:
                /*! Error message */
                virtual const char * what() const throw() {
                    return "Trace not initialized.";
                }
        };
    }
}

#endif
