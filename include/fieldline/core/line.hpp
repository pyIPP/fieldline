#ifndef include_fieldline_core_line_hpp
#define include_fieldline_core_line_hpp

#include <math.h>
#include <fieldline/core/point.hpp>


namespace Fieldline {
    namespace core {
        /*! \brief Class representing a line in the poloidal plain
         *
         * This class represents the line with the two points \f$(R_0, z_0\f$ and \f$(R_1, z_1)\f$.
         */
        class line {
            public:
                /*! \brief Constructor
                 *
                 *  This constructor creates the line with explicit R and z values.
                 */
                line(const double R0, const double z0, const double R1, const double z1) : m_R0(R0), m_z0(z0), m_R1(R1), m_z1(z1), m_dR(R1-R0), m_dz(z1-z0) {}
                /*! \brief Copy constructor */
                line(const line & rhs) : m_R0(rhs.m_R0), m_z0(rhs.m_z0), m_R1(rhs.m_R1), m_z1(rhs.m_z1), m_dR(rhs.m_dR), m_dz(rhs.m_dz) {}
                /*! \brief Constructor
                 *
                 * This constructor creates the line represented by two points.
                 */
                line(const Fieldline::core::point & p0, const Fieldline::core::point & p1) : m_R0(p0.R), m_z0(p0.z), m_R1(p1.R), m_z1(p1.z), m_dR(p1.R-p0.R), m_dz(p1.z-p0.z) {}
                /*! \brief Destructor */
                virtual ~line() {}
                /*! \brief Copy operator
                 *
                 *  This operator copies the line to the given instance.
                 */
                line & operator= (const line & rhs) {
                    if(this != &rhs) {
                        m_R0 = rhs.m_R0;
                        m_z0 = rhs.m_z0;
                        m_R1 = rhs.m_R1;
                        m_z1 = rhs.m_z1;
                        m_dR = rhs.m_dR;
                        m_dz = rhs.m_dz;
                    }
                    return *this;
                }
                /*! \brief Get the line intersection
                 *
                 *  This function returns the intersection point between the two lines.
                 *  If the lines are parallel the origin is returned with the information that the lines did not hit.
                 *  If the intersections lies between the points defining each line, point.hit is set to true.
                 */
                Fieldline::core::point get_intersection(const line & rhs) const {
                    double dR = rhs.m_R0 - m_R0;
                    double dz = rhs.m_z0 - m_z0;
                    double D = m_dR*(-rhs.m_dz) + rhs.m_dR*m_dz;
                    if(D==0) {
                        return Fieldline::core::point(0.0, 0.0, false);
                    }
                    double Dx = dR*(-rhs.m_dz) + dz*rhs.m_dR;
                    double Dy = m_dR*dz - m_dz*dR;
                    double x = Dx/D;
                    double y = Dy/D;
                    return Fieldline::core::point(m_R0 + x*m_dR, m_z0 + x*m_dz, (x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0));
                }
                /*! \brief Get the R value of the first point */
                inline double get_R0() const { return m_R0; }
                /*! \brief Get the z value of the first point */
                inline double get_z0() const { return m_z0; }
                /*! \brief Get the R value of the second point */
                inline double get_R1() const { return m_R1; }
                /*! \brief Get the z value of the second point */
                inline double get_z1() const { return m_z1; }
                /*! \brief Get the length of the line
                 *
                 * This function returns the length of the line: \f$ d = \sqrt{(R_1 - R_0)^2 + (z_1 - z_0)^2} \f$
                 */
                inline double get_length() const {
                    return sqrt(m_dR*m_dR + m_dz*m_dz);
                }
            protected:
                double m_R0; /*!< \brief R value of the first point. */
                double m_z0; /*!< \brief z value of the first point. */
                double m_R1; /*!< \brief R value of the second point. */
                double m_z1; /*!< \brief z value of the second point. */
                double m_dR; /*!< \brief Distance between the two points in R: \f$ R_1 - R_0 \f$ */
                double m_dz; /*!< \brief Distance between the two points in z: \f$ z_1 - z_0 \f$ */
        };
    }
}



#endif
