#ifndef include_fieldline_exceptions_traceNotInitialized_hpp
#define include_fieldline_exceptions_traceNotInitialized_hpp

#include <exception>


namespace Fieldline {
    namespace exceptions {
        /*! \brief Exception if the field line trace is not initialized.
         *
         *  This exception is thrown in the Runge Kutta field line trace if no start point is defined.
         */
        class traceNotInitialized : public std::exception {
            public:
                /*! \brief Error message */
                virtual const char * what() const throw() {
                    return "Trace not initialized.";
                }
        };
    }
}

#endif
