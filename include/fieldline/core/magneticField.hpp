#ifndef include_fieldline_core_magneticField_hpp
#define include_fieldline_core_magneticField_hpp

#include <stdlib.h>
#include <iostream>
#include <math.h>

namespace fieldline {
    namespace core {
        class magneticField {
            public:
                magneticField() : BR(0.0), Bz(0.0), Btor(0.0) {}
                magneticField(const double br, const double bz, const double btor) : BR(br), Bz(bz), Btor(btor) {}
                magneticField(const magneticField & rhs) : BR(rhs.BR), Bz(rhs.Bz), Btor(rhs.Btor) {}
                virtual ~magneticField() {}

                magneticField & operator= (const magneticField & rhs) {
                    if(this != &rhs) {
                        BR = rhs.BR;
                        Bz = rhs.Bz;
                        Btor = rhs.Btor;
                    }
                    return *this;
                }
                double Btot() const {
                    return sqrt(BR*BR + Bz*Bz + Btor*Btor);
                }
                double Bpol() const {
                    return sqrt(BR*BR + Bz*Bz);
                }
                double BR;
                double Bz;
                double Btor;

        };

        std::ostream & operator<<(std::ostream & ostr, const magneticField & field) {
            ostr << "BR: " << field.BR << "\tBz: " << field.Bz << "\tBtor: " << field.Btor << '\n';
            return ostr;
        }


    }
}

#endif
