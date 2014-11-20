#ifndef include_fieldline_exceptions_traceNotInitialized_hpp
#define include_fieldline_exceptions_traceNotInitialized_hpp

#include <exception>

namespace fieldline {
    namespace exceptions {
        class traceNotInitialized : public std::exception {
            public:
                virtual const char * what() const throw() {
                    return "Trace not initialized.";
                }
        }traceNotInitialized;
    }
}

#endif
