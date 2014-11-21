#ifndef include_fieldline_core_fieldline_hpp
#define include_fieldline_core_fieldline_hpp

#include <vector>
#include <math.h>
#include <fstream>
#include <string.h>
#include <fieldline/core/line.hpp>

#ifdef WITH_PYTHON
    #include <boost/python/list.hpp>
#endif

namespace Fieldline {
    namespace core {
        class fieldline {
            public:
                fieldline() {}
                fieldline(const fieldline & rhs) : m_R(rhs.m_R), m_z(rhs.m_z), m_phi(rhs.m_phi), m_Btot(rhs.m_Btot) {}
                fieldline(const std::vector<double> R, const std::vector<double> z, const std::vector<double> phi, const std::vector<double> Btot) :
                    m_R(R), m_z(z), m_phi(phi), m_Btot(Btot) {}
                fieldline(const std::string & filename) : m_R(), m_z(), m_phi(), m_Btot() {
                    read_from_file(filename);
                }
#ifdef WITH_PYTHON
                fieldline(const boost::python::list & R, const boost::python::list & z, const boost::python::list & phi, const boost::python::list & Btot) :
                    m_R(), m_z(), m_phi(), m_Btot() {
                    uint32_t N = boost::python::len(R);
                    m_R.resize(N);
                    m_z.resize(N);
                    m_phi.resize(N);
                    m_Btot.resize(N);
                    for(uint32_t i = 0; i < N; ++i) {
                        m_R[i] = boost::python::extract<double>(R[i]);
                        m_z[i] = boost::python::extract<double>(z[i]);
                        m_phi[i] = boost::python::extract<double>(phi[i]);
                        m_Btot[i] = boost::python::extract<double>(Btot[i]);
                    }
                }
#endif
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
#ifdef WITH_PYTHON
                boost::python::list get_x_python() const {
                    std::vector<double> temp = get_x();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                boost::python::list get_y_python() const {
                    std::vector<double> temp = get_y();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                boost::python::list get_z_python() const {
                    std::vector<double> temp = get_z();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                boost::python::list get_R_python() const{
                    boost::python::list output;
                    for(auto iter = m_R.cbegin(); iter != m_R.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                boost::python::list get_phi_python() const{
                    boost::python::list output;
                    for(auto iter = m_phi.cbegin(); iter != m_phi.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                boost::python::list get_Btot_python() const{
                    boost::python::list output;
                    for(auto iter = m_Btot.cbegin(); iter != m_Btot.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

#endif

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
                void write_to_file(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    file << m_R.size() << '\n';
                    for(uint32_t i = 0; i < m_R.size(); ++i) {
                        file << m_R[i] << '\t' << m_z[i] << '\t' << m_phi[i] << '\t' << m_Btot[i] << '\n';
                    }
                    file.close();
                }
                void read_from_file(const std::string & filename) {
                    std::fstream file(filename.c_str(), std::ios::in);
                    uint32_t N;
                    file >> N;
                    m_R.resize(N);
                    m_z.resize(N);
                    m_phi.resize(N);
                    m_Btot.resize(N);
                    auto R = m_R.begin();
                    auto z = m_z.begin();
                    auto phi = m_phi.begin();
                    auto Btot = m_Btot.begin();
                    for( ; R != m_R.end(); ++R, ++z, ++phi, ++Btot) {
                        file >> *R >> *z >> *phi >> *Btot;
                    }
                    file.close();
                }

                Fieldline::core::line get_end() const {
                    uint32_t N = m_R.size();
                    return Fieldline::core::line(m_R[N-2], m_z[N-2], m_R[N-1], m_z[N-1]);
                }
                void replace_end(const double R, const double z) {
                    uint32_t N = m_R.size();
                    double dR0 = m_R[N-2]-R;
                    double dz0 = m_z[N-2]-z;
                    double dR1 = m_R[N-2]-m_R[N-1];
                    double dz1 = m_z[N-2]-m_z[N-1];
                    double d0 = sqrt(dR0*dR0 + dz0*dz0);
                    double d1 = sqrt(dR1*dR1 + dz1*dz1);
                    m_R[N-1] = R;
                    m_z[N-1] = z;
                    m_Btot[N-1] = m_Btot[N-2] + d0/d1*(m_Btot[N-1]-m_Btot[N-2]);
                    m_phi[N-1] = m_phi[N-2] + d0/d1*(m_phi[N-1]-m_phi[N-2]);
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
