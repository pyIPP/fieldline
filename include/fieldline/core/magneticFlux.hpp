#ifndef include_fieldline_core_magneticFlux_hpp
#define include_fieldline_core_magneticFlux_hpp


namespace Fieldline {
    namespace core {
        /*! \brief Class containing the poloidal magnetic flux \f$\psi\f$ and its derivatives. */
        class magneticFlux {
            public:
                /*! \brief Default constructor 
                 *
                 * Sets the magnetic flux and its derivatives to zero.
                */
                magneticFlux() : psi(0.0), dpsi_dR(0.0), dpsi_dz(0.0) {}
                /*! \brief Constructor
                 *
                 *  Initializes the magnetic flux with the given value and derivatives.
                 */
                magneticFlux(const double Psi, const double dPsi_dR, const double dPsi_dz) : psi(Psi), dpsi_dR(dPsi_dR), dpsi_dz(dPsi_dz) {}
                /*! \brief Copy constructor */
                magneticFlux(const magneticFlux & rhs) : psi(rhs.psi), dpsi_dR(rhs.dpsi_dR), dpsi_dz(rhs.dpsi_dz) {}
                /*! \brief Assignment operator
                 *
                 *  Copies the given magnetic flux to the current instance.
                 */
                magneticFlux & operator= (const magneticFlux & rhs) {
                    if(this != &rhs) {
                        psi = rhs.psi;
                        dpsi_dR = rhs.dpsi_dR;
                        dpsi_dz = rhs.dpsi_dz;
                    }
                    return *this;
                }

                double psi; /*!< Poloidal magnetic flux \f$\psi\f$.*/
                double dpsi_dR; /*!< Gradient of the magnetic flux in R direction \f$ \frac{\partial \psi}{\partial R} \f$.*/
                double dpsi_dz; /*!< Gradient of the magnetic flux in z direction \f$ \frac{\partial \psi}{\partial z} \f$.*/
        };
    }
}


#endif
