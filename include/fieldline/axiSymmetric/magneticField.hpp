#ifndef include_fieldline_axiSymmetric_magneticField_hpp
#define include_fieldline_axiSymmetric_magneticField_hpp

#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <fieldline/exceptions.hpp>
#include <fieldline/core/magneticField.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include <string.h>

namespace fieldline {
    namespace axiSymmetric {
        class magneticField {
            private: 
                const double m_a1[3] = {3.0, -1.5, 1.0/3.0};
                const double m_a2[3] = {-2.5, 2.5, 0.5};
                const double m_a3[3] = {0.5, -0.5, 1.0/6.0};
            public:
                magneticField() {}
                magneticField(const double Btor, const double R0, const double Rmin, const double Rmax, 
                    const double zmin, const double zmax, const boost::numeric::ublas::matrix<double> psi) :
                        m_Btor(Btor), m_R0(R0), m_Rmin(Rmin), m_Rmax(Rmax), m_zmin(zmin), m_zmax(zmax),
                        m_dR((Rmax-Rmin)/(psi.size1()-1)), m_dz((zmax-zmin)/(psi.size2()-1)), m_psi(psi) {
                }
                magneticField(const double Btor, const double R0, const std::string & filename) : m_Btor(Btor), m_R0(R0) {
                    std::fstream file(filename.c_str(), std::ios::in);
                    load_from_file(file);
                }
                magneticField(const double Btor, const double R0, std::fstream & file) : m_Btor(Btor), m_R0(R0) {
                    load_from_file(file);
                }
                magneticField(const magneticField & rhs) : m_Btor(rhs.m_Btor), m_R0(rhs.m_R0), m_NR(rhs.m_NR), m_Nz(rhs.m_Nz),
                    m_Rmin(rhs.m_Rmin), m_Rmax(rhs.m_Rmax), m_zmin(rhs.m_zmin), m_zmax(rhs.m_zmax), m_dR(rhs.m_dR), m_dz(rhs.m_dz), m_psi(rhs.m_psi) {
                }

                virtual ~magneticField() {}

                fieldline::core::magneticField get_magnetic_field(const double R, const double z) const {
                    if(R < m_Rmin || R > m_Rmax || z < m_zmin || z > m_zmax) {
                        throw fieldline::exceptions::undefinedMagneticField;
                    }
                    double x = (R - m_Rmin)/m_dR;
                    double y = (z - m_zmin)/m_dz;
                    uint32_t jr = std::min(std::max((uint32_t)0, (uint32_t)x-1), m_NR-4);
                    uint32_t jz = std::min(std::max((uint32_t)0, (uint32_t)y-1), m_Nz-4);
                    x -= jr;
                    y -= jz;
                    double PF[4];
                    double x2 = x*x;
                    double y2 = y*y;
                    double x3 = x2*x;
                    double y3 = y3*y;
                    double F[4];
                    double dzdy[4];
                    for(uint32_t i = 0; i < 4; ++i) {

                    }
                    return fieldline::core::magneticField(0.0, 0.0, m_R0/R*m_Btor);
                }


            protected:

                double m_Btor;
                double m_R0;
                uint32_t m_NR;
                uint32_t m_Nz;
                double m_Rmin;
                double m_Rmax;
                double m_zmin;
                double m_zmax;
                double m_dR;
                double m_dz;
                boost::numeric::ublas::matrix<double> m_psi;
                void load_from_file(std::fstream & file) {
                    std::string temp;
                    file >> m_NR >> m_Nz;
                    m_psi = boost::numeric::ublas::matrix<double>(m_NR, m_Nz);
                    std::getline(file, temp);
                    file >> m_Rmin >> m_Rmax >> m_zmin >> m_zmax;
                    m_dR = (m_Rmax - m_Rmin)/(m_NR-1);
                    m_dz = (m_zmax - m_zmin)/(m_Nz-1);
                    for(uint32_t i = 0; i < m_NR; ++i) {
                        for(uint32_t j = 0; j < m_Nz; ++j) {
                            file >> m_psi(i,j);
                        }
                    }
                }

        };
    }
}


#endif 
