#ifndef include_fieldline_trace_rungeKutta_hpp
#define include_fieldline_trace_rungeKutta_hpp

#include <stdint.h>
#include <fieldline/core.hpp>
#include <fieldline/exceptions.hpp>
#include <fieldline/axiSymmetric.hpp>

namespace fieldline {
    namespace trace {
        class rungeKutta {
            private:
                const double m_c[4] = {1.0/3.0, 2.0/3.0, 1.0/3.0, 1.0/6.0};
                const double m_a[4] = {0.5, 0.5, 1.0, 1.0/6.0};
                const double m_cm[4] = {-3.0/8.0, 3.7e1/2.4e1, -5.9e1/2.4e1,5.5e1/2.4e1};
                const double m_fm[5] = {-1.9/7.2e1, 5.3/3.6e1,-1.1/3.0, 1.615e1/1.8e1, 2.51e1/7.2e1};
            public:
                rungeKutta(fieldline::axiSymmetric::magneticField * magneticField) : m_magneticField(magneticField), m_initialized(false) {}
                virtual ~rungeKutta() {}

                void init(const double R, const double z, const double phi, const double dphi) {
                    fieldline::core::magneticField magneticField;
                    double qR[4];
                    double qz[4];
                    double fa[4];
                    m_phi = phi;
                    m_R = R;
                    m_z = z;
                    for(uint32_t i = 0; i < 3; ++i) {
                        for(uint32_t j = 0; j < 3; ++j) {
                            magneticField = m_magneticField->get_magnetic_field(m_R, m_z, phi);
                            fa[j] = R*dphi*m_a[j]/magneticField.Btor;
                            qR[j] = magneticField.BR*fa[j];
                            qz[j] = magneticField.Bz*fa[j];
                            m_R = R + qR[j];
                            m_z = z + qz[j];
                            m_phi = phi + m_a[j]*dphi;
                        }
                        magneticField = m_magneticField->get_magnetic_field(m_R, m_z, phi);
                        m_dR[i] = 2.0*qR[0];
                        m_dz[i] = 2.0*qz[0];
                        fa[3] = m_R*m_a[3]*dphi/magneticField.Btor;
                        m_R += magneticField.BR*fa[3];
                        m_z += magneticField.Bz*fa[3];
                        for(uint32_t j = 0; j < 3; ++j) {
                            m_R += qR[j]*m_c[j];
                            m_z += qz[j]*m_c[j];
                        }
                        m_phi = phi;
                    }
                    magneticField = m_magneticField->get_magnetic_field(m_R, m_z, m_phi);
                    m_dR[3] = magneticField.BR*m_R*dphi/magneticField.Btor;
                    m_dz[3] = magneticField.Bz*m_z*dphi/magneticField.Btor;
                    m_Btot = magneticField.Btot();
                    m_initialized = true;
                }
                void next(const uint32_t nSteps, const double dphi) {
                    if(!m_initialized) {
                        throw fieldline::exceptions::traceNotInitialized;
                    }
                    fieldline::core::magneticField magneticField;
                    double tempR;
                    double tempz;
                    for(uint32_t i = 0; i < nSteps; ++i) {
                        tempR = m_R;
                        tempz = m_z;
                        for(uint32_t j = 0; j < 4; ++j) {
                            tempR += m_cm[j]*m_dR[j];
                            tempz += m_cm[j]*m_dz[j];
                        }
                        m_phi += dphi;
                        magneticField = m_magneticField->get_magnetic_field(tempR, tempz, m_phi);
                        m_dR[4] = magneticField.BR*dphi*m_R/magneticField.Btor;
                        m_dz[4] = magneticField.Bz*dphi*m_R/magneticField.Btor;
                        for(uint32_t j = 0; j < 5; ++j) {
                            m_R += m_fm[j]*m_dR[j];
                            m_z += m_fm[j]*m_dz[j];
                        }
                        for(uint32_t j = 0; j < 4; ++j) {
                            m_dR[j] = m_dR[j+1];
                            m_dz[j] = m_dz[j+1];
                        }
                    }
                    magneticField = m_magneticField->get_magnetic_field(m_R, m_z, m_phi);
                    m_Btot = magneticField.Btot();
                }
                double get_R() const { return m_R; };
                double get_z() const { return m_z; };
                double get_phi() const { return m_phi; };
                double get_Btot() const { return m_Btot; };



            protected:
                fieldline::axiSymmetric::magneticField * m_magneticField;
                bool m_initialized;
                double m_phi;
                double m_R;
                double m_z;
                double m_dR[5];
                double m_dz[5];
                double m_Btot;


        };
    }
}


#endif