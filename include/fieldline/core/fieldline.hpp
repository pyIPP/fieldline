#ifndef include_fieldline_core_fieldline_hpp
#define include_fieldline_core_fieldline_hpp

#include <vector>
#include <math.h>
#include <fstream>
#include <string.h>

namespace Fieldline {
    namespace core {
        class fieldline {
            public:
                fieldline() {}
                fieldline(const fieldline & rhs) : m_R(rhs.m_R), m_z(rhs.m_z), m_phi(rhs.m_phi), m_Btot(rhs.m_Btot) {}
                fieldline(const std::vector<double> R, const std::vector<double> z, const std::vector<double> phi, const std::vector<double> Btot) :
                    m_R(R), m_z(z), m_phi(phi), m_Btot(Btot) {}
                virtual ~fieldline() {}
                fieldline & operator= (const fieldline & rhs) {
                    if(this != &rhs){
                        m_R = rhs.m_R;
                        m_z = rhs.m_z;
                        m_phi = rhs.m_phi;
                        m_Btot = rhs.m_Btot;
                    }
                    return *this;
                }

                uint32_t size() const { return m_R.size(); }
                void push_back(const double R, const double z, const double phi, const double Btot) {
                    m_R.push_back(R);
                    m_z.push_back(z);
                    m_phi.push_back(phi);
                    m_Btot.push_back(Btot);
                }
                std::vector<double> get_x() const {
                    std::vector<double> x(m_R.size());
                    auto R = m_R.cbegin();
                    auto phi = m_phi.cbegin();
                    auto iter = x.begin();
                    for(;R != m_R.cend(); ++R, ++phi, ++iter) {
                        *iter = *R*cos(*phi);
                    }
                    return x;
                }
                std::vector<double> get_y() const {
                    std::vector<double> y(m_R.size());
                    auto R = m_R.cbegin();
                    auto phi = m_phi.cbegin();
                    auto iter = y.begin();
                    for(;R != m_R.cend(); ++R, ++phi, ++iter) {
                        *iter = *R*sin(*phi);
                    }
                    return y;
                }
                std::vector<double> get_z() const {
                    return m_z;
                }
                std::vector<double> get_R() const {
                    return m_R;
                }
                std::vector<double> get_phi() const {
                    return m_phi;
                }
                std::vector<double> get_Btot() const {
                    return m_Btot;
                }
                double get_length() const {
                    double output = 0.0;
                    std::vector<double> x = get_x();
                    std::vector<double> y = get_y();
                    for(uint32_t i = 0; i < x.size()-1; ++i) {
                        output += (x[i]*x[i] + y[i]*y[i] + m_z[i]*m_z[i]);
                    }
                    return fabs(output);
                }
                double get_last_R() const {
                    return m_R[m_R.size()-1];
                }
                double get_last_z() const {
                    return m_z[m_z.size()-1];
                }
                void write_ASCII(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    for(uint32_t i = 0; i < m_R.size(); ++i) {
                        file << m_R[i] << '\t' << m_z[i] << '\t' << m_phi[i] << '\t' << m_Btot[i] << '\n';
                    }
                    file.close();
                }

                

            protected:
                std::vector<double> m_R;
                std::vector<double> m_z;
                std::vector<double> m_phi;
                std::vector<double> m_Btot;
        };
    }
}

#endif
