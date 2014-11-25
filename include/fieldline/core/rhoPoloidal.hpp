#ifndef include_fieldline_core_rhoPoloidal_hpp
#define include_fieldline_core_rhoPoloidal_hpp

#include <math.h>

namespace Fieldline {
    namespace core {
        /*! \brief Class for the calculation of rho poloidal.
         *
         * This class is a generic class for the calculation of rho poloidal.
         */
        class rhoPoloidal {
            public:
                /*! \brief Constructor
                 *  
                 *  The constructor takes the poloidal magnetic flux \f$\psi\f$ at the core and the edge of the confined plasma.
                 *  \param psi_core poloidal magnetic flux \f$\psi_{core}\f$ at the core of the confined plasma
                 *  \param psi_edge poloidal magnetic flux \f$\psi_{edge}\f$ at the edge of the confined plasma
                 */ 
                rhoPoloidal(const double psi_core, const double psi_edge) : m_psiCore(psi_core), m_psiEdge(psi_edge), m_dpsi(psi_edge-psi_core) {}
                /*! \brief Copy constructor */
                rhoPoloidal(const rhoPoloidal & rhs) : m_psiCore(rhs.m_psiCore), m_psiEdge(rhs.m_psiEdge), m_dpsi(rhs.m_dpsi) {}
                /*! \brief Destructor */
                virtual ~rhoPoloidal() {}

                /*! \brief Get rho poloidal
                 *
                 * This function calculates rho poloidal for the given poloidal magnetic flux 
                 * \f$ \rho_{pol} = \sqrt{\frac{\psi - \psi_{core}}{\psi_{edge}-\psi_{core}}}\f$
                 */
                double get_rho_poloidal(const double psi) const {
                    return sqrt((psi-m_psiCore)/m_dpsi);
                }
                /*! \brief Get poloidal magnetic flux
                 *
                 * This function calculates the poloidal magnetic flux for a given rho poloidal
                 * \f$ \psi = \rho_{pol}^2 \cdot (\psi_{edge}-\psi_{core}) + \psi_{core}\f$
                 */
                double get_poloidal_flux(const double rho) const {
                    return rho*rho*m_dpsi + m_psiCore;
                }

            protected:
                double m_psiCore; /*!< \brief Poloidal magnetic flux \f$\psi_{core}\f$ at the core of the plasma. */
                double m_psiEdge; /*!< \brief Poloidal magnetic flux \f$\psi_{edge}\f$ at the edge of the plasma. */
                double m_dpsi; /*!< \brief Difference in poloidal magnetic flux between the edge and the core \f$\psi_{edge} - \psi_{core}\f$. */
        };
    }
}

#endif
