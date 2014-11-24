#ifndef include_fieldline_core_magneticField_hpp
#define include_fieldline_core_magneticField_hpp

#include <stdlib.h>
#include <iostream>
#include <math.h>

/*! Class containing one magnetic field vector
    This class contains the magnetic field vector.
*/ 

namespace Fieldline {
    namespace core {
        class magneticField {
            public:
                /*! Default constructor
                    The default constructor initialized the instance without magnetic field. All components are zero.
                */
                magneticField() : BR(0.0), Bz(0.0), Btor(0.0) {}
                /*! Initializing constructor
                    This constructor takes the magnetic field components in R,z and toroidal direction.
                */
                magneticField(const double br, const double bz, const double btor) : BR(br), Bz(bz), Btor(btor) {}
                /*! Copy constructor */
                magneticField(const magneticField & rhs) : BR(rhs.BR), Bz(rhs.Bz), Btor(rhs.Btor) {}
                /*! Destructor */
                virtual ~magneticField() {}

                /*! Assignment operator
                    This operator copies the magnetic field vector from the given instance.
                */
                magneticField & operator= (const magneticField & rhs) {
                    if(this != &rhs) {
                        BR = rhs.BR;
                        Bz = rhs.Bz;
                        Btor = rhs.Btor;
                    }
                    return *this;
                }
                /*! Get total magnetic field
                    This function returns the total magnetic field \f$B_{tot} = \sqrt{B_R^2 + B_z^2 + B_{tor}^2}\f$
                */
                double Btot() const {
                    return sqrt(BR*BR + Bz*Bz + Btor*Btor);
                }
                /*! Get poloidal magnetic field
                    This function returns the poloidal magnetic field \f$B_{pol} = \sqrt{B_R^2 + B_z^2} \f$ 
                */
                double Bpol() const {
                    return sqrt(BR*BR + Bz*Bz);
                }
                /*! Magnetic field component in R direction. */
                double BR;
                /*! Magnetic field component in z direction. */
                double Bz;
                /*! Magnetic field component in toroidal direction. */
                double Btor;

        };

        std::ostream & operator<<(std::ostream & ostr, const magneticField & field) {
            ostr << "BR: " << field.BR << "\tBz: " << field.Bz << "\tBtor: " << field.Btor << '\n';
            return ostr;
        }


    }
}

#endif
