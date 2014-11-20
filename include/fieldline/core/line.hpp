#ifndef include_fieldline_core_line_hpp
#define include_fieldline_core_line_hpp

#include <fieldline/core/point.hpp>

namespace Fieldline {
    namespace core {
        class line {
            public:
                line(const double R0, const double z0, const double R1, const double z1) : m_R0(R0), m_z0(z0), m_R1(R1), m_z1(z1), m_dR(R1-R0), m_dz(z1-z0) {}
                line(const line & rhs) : m_R0(rhs.m_R0), m_z0(rhs.m_z0), m_R1(rhs.m_R1), m_z1(rhs.m_z1), m_dR(rhs.m_dR), m_dz(rhs.m_dz) {}
                virtual ~line() {}
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
                Fieldline::core::point intersection(const line & rhs) const {
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
                inline double get_R0() const { return m_R0; }
                inline double get_z0() const { return m_z0; }
            protected:
                double m_R0;
                double m_z0;
                double m_R1;
                double m_z1;
                double m_dR;
                double m_dz;
        };
    }
}



#endif
