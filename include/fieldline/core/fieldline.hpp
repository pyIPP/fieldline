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

/*! Class containing the field line information.
    This class stores the information of the field line.
    The trace of the field line is stored in $\f R,z,\phi$\f.
    In addition the total magnetic field at each position is stored.
*/

namespace Fieldline {
    namespace core {
        class fieldline {
            public:
                /*! Default constructor
                    This constructor initializes an empty field line.
                */
                fieldline() {}
                /*! Copy constructor */
                fieldline(const fieldline & rhs) : m_R(rhs.m_R), m_z(rhs.m_z), m_phi(rhs.m_phi), m_Btot(rhs.m_Btot) {}
                /*! Constructor
                    This constructor initializes the field line with the given coordinates \f$ R,z,\hp\f$ and the total magnetic field.
                */
                fieldline(const std::vector<double> R, const std::vector<double> z, const std::vector<double> phi, const std::vector<double> Btot) :
                    m_R(R), m_z(z), m_phi(phi), m_Btot(Btot) {}
                /*! File constructor
                    This constructor read the field line from the specified file.
                */
                fieldline(const std::string & filename) : m_R(), m_z(), m_phi(), m_Btot() {
                    read_from_file(filename);
                }
#ifdef WITH_PYTHON
                /*! Python constructor
                    This constructor initializes the field line with the given coordinates \f$ R,z,\hp\f$ and the total magnetic field.
                    This constructor is intended as an interface to python and should never be used from within C++.
                */
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
                /*! Destructor */
                virtual ~fieldline() {}
                /*! Assignment operator.
                    This operator copies the field line from the given instance.
                */
                fieldline & operator= (const fieldline & rhs) {
                    if(this != &rhs){
                        m_R = rhs.m_R;
                        m_z = rhs.m_z;
                        m_phi = rhs.m_phi;
                        m_Btot = rhs.m_Btot;
                    }
                    return *this;
                }
   
                /*! Size
                    This function returns the number of points on the field line.
                */
                uint32_t size() const { return m_R.size(); }
                /*! Append point at the end of the field line
                    This function adds the point \f$(R,z,\phi)\f$ with the magnetic field \f$B_{tot}\f$ at the end of the field line.
                */
                void push_back(const double R, const double z, const double phi, const double Btot) {
                    m_R.push_back(R);
                    m_z.push_back(z);
                    m_phi.push_back(phi);
                    m_Btot.push_back(Btot);
                }
                /*! Return the x coordinates of the field line.
                    This function returns the x coordinates of the field line: \f$ x = R \cos(\phi)\f$
                */
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
                /*! Return the y coordinates of the field line.
                    This function returns the y coordinates of the field line: \f$ y = R \sin(\phi)\f$
                */
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
                /*! Return the z coordinates of the field line 
                    This function returns the z coordinates of the field line
                */
                std::vector<double> get_z() const {
                    return m_z;
                }
                /*! Return the R coordinates of the field line 
                    This function returns the R coordinates of the field line
                */
                std::vector<double> get_R() const {
                    return m_R;
                }
                /*! Return the toroidal angle of the field line
                    This function returns the toroidal angle \f$\phi\f$ of the field line.
                */
                std::vector<double> get_phi() const {
                    return m_phi;
                }
                /*! Return the total magnetic field of the field line
                    This function returns the total magnetic field \f$B_{tot}\f$ of the field line.
                */
                std::vector<double> get_Btot() const {
                    return m_Btot;
                }
#ifdef WITH_PYTHON
                /*! Return the x coordinates of the field line as python list.
                    This function returns the x coordinates of the field line: \f$ x = R \cos(\phi)\f$
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_x_python() const {
                    std::vector<double> temp = get_x();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! Return the y coordinates of the field line as python list.
                    This function returns the y coordinates of the field line: \f$ y = R \sin(\phi)\f$
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_y_python() const {
                    std::vector<double> temp = get_y();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! Return the z coordinates of the field line as python list.
                    This function returns the z coordinates of the field line.
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_z_python() const {
                    std::vector<double> temp = get_z();
                    boost::python::list output;
                    for(auto iter = temp.cbegin(); iter != temp.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! Return the R coordinates of the field line as python list.
                    This function returns the R coordinates of the field line.
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_R_python() const{
                    boost::python::list output;
                    for(auto iter = m_R.cbegin(); iter != m_R.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! Return the toroidal angle of the field line as python list.
                    This function returns the toroidal angle \f$\phi\f$ of the field line.
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_phi_python() const{
                    boost::python::list output;
                    for(auto iter = m_phi.cbegin(); iter != m_phi.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! Return the total magnetic field of the field line as python list
                    This function returns the total magnetic field \f$B_{tot}\f$ of the field line.
                    This function is intended as a python interface, never use this function from within python.
                */
                boost::python::list get_Btot_python() const{
                    boost::python::list output;
                    for(auto iter = m_Btot.cbegin(); iter != m_Btot.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

#endif
                /*! Field line length
                    This function returns the length of the field line.
                    \f$L = \sum\limits_0^{N-1} \sqrt{(x[i+1]-x[i])^2 + (y[i+1]-y[i])^2 + (z[i+1]-z[y])^2} \f$ 
                */
                double get_length() const {
                    double output = 0.0;
                    std::vector<double> x = get_x();
                    std::vector<double> y = get_y();
                    for(uint32_t i = 0; i < x.size()-1; ++i) {
                        output += sqrt(pow(x[i+1]-x[i],2) + pow(y[i+1]-y[i],2) + pow(m_z[i+1]-m_z[i],2));
                    }
                    return fabs(output);
                }
                /*! Last R value in field line.
                    This function returns the last R value in the field line.
                */
                double get_last_R() const {
                    return m_R[m_R.size()-1];
                }
                /*! Last z value in field line.
                    This function returns the last z value in the field line.
                */
                double get_last_z() const {
                    return m_z[m_z.size()-1];
                }
                /*! Write field line to file.
                    This function writes the field line to an ASCII file.
                */
                void write_to_file(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    file << m_R.size() << '\n';
                    for(uint32_t i = 0; i < m_R.size(); ++i) {
                        file << m_R[i] << '\t' << m_z[i] << '\t' << m_phi[i] << '\t' << m_Btot[i] << '\n';
                    }
                    file.close();
                }
                /*! Read field line from file.
                    This function reads the field line from an ASCII file.
                */
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

                /*! Return the last line segment of the field line.
                    This function returns the last line segment of the field line.
                */
                Fieldline::core::line get_end() const {
                    uint32_t N = m_R.size();
                    return Fieldline::core::line(m_R[N-2], m_z[N-2], m_R[N-1], m_z[N-1]);
                }
                /*! Replace the end of the field line
                    This functions replaces the last point on the field line with the specified R,z coordinates.
                    The total magnetic field and the toroidal angle are interpolated accordingly.
                */
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
                /*! R coordinates of the field line */
                std::vector<double> m_R;
                /*! z coordinates of the field line */
                std::vector<double> m_z;
                /*! Toroidal angle of the field line */
                std::vector<double> m_phi;
                /*! Total magnetic field of the field line */
                std::vector<double> m_Btot;
        };
    }
}

#endif
