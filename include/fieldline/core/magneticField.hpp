#ifndef include_fieldline_core_magneticField_hpp
#define include_fieldline_core_magneticField_hpp

#include <stdlib.h>
#include <iostream>
#include <math.h>


namespace Fieldline {
    namespace core {
        /*! \brief Class containing one magnetic field vector
         *
         * This class contains the magnetic field vector.
         */ 
        class magneticField {
            public:
                /*! \brief Default constructor
                 *
                 *  The default constructor initialized the instance without magnetic field. All components are zero.
                 */
                magneticField() : BR(0.0), Bz(0.0), Btor(0.0) {}
                /*! \brief Initializing constructor
                 *
                 * This constructor takes the magnetic field components in R,z and toroidal direction.
                 */
                magneticField(const double br, const double bz, const double btor) : BR(br), Bz(bz), Btor(btor) {}
                /*! \brief Copy constructor */
                magneticField(const magneticField & rhs) : BR(rhs.BR), Bz(rhs.Bz), Btor(rhs.Btor) {}
                /*! \brief Destructor */
                virtual ~magneticField() {}

                /*! \brief Assignment operator
                 *
                 * This operator copies the magnetic field vector from the given instance.
                 */
                magneticField & operator= (const magneticField & rhs) {
                    if(this != &rhs) {
                        BR = rhs.BR;
                        Bz = rhs.Bz;
                        Btor = rhs.Btor;
                    }
                    return *this;
                }
                /*! \brief Get total magnetic field \f$B_{tot}\f$.
                 *
                 *  This function returns the total magnetic field \f$B_{tot} = \sqrt{B_R^2 + B_z^2 + B_{tor}^2}\f$
                 */
                double Btot() const {
                    return sqrt(BR*BR + Bz*Bz + Btor*Btor);
                }
                /*! \brief Get poloidal magnetic field \f$B_{pol}\f$.
                 *
                 *  This function returns the poloidal magnetic field \f$B_{pol} = \sqrt{B_R^2 + B_z^2} \f$ 
                 */
                double Bpol() const {
                    return sqrt(BR*BR + Bz*Bz);
                }
                double BR; /*!< Magnetic field component in R direction \f$B_R\f$.*/
                double Bz; /*!< Magnetic field component in z direction \f$B_z\f$.*/
                double Btor; /*!< Magnetic field component in toroidal direction \f$B_{tor}\f$.*/

                /*! \brief Get toroidal inclination of the magnetic field.
                 *
                 *  This function returns the toroidal inclination of the magnetic field
                 *  \f$ \alpha = \pi - \arctan\left(\frac{B_{pol}}{B_{tor}}\right) \f$
                 *  The angle is zero if only the toroidal field is present.
                 */
                double toroidal_inclination() const {
                    return M_PI - atan2(Bpol(), Btor);
                }

        };

        std::ostream & operator<<(std::ostream & ostr, const magneticField & field) {
            ostr << "BR: " << field.BR << "\tBz: " << field.Bz << "\tBtor: " << field.Btor << '\n';
            return ostr;
        }


    }
}

#endif
