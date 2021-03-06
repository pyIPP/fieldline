#ifndef include_fieldline_axiSymmetric_magneticField_hpp
#define include_fieldline_axiSymmetric_magneticField_hpp

#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <exception>
#include <fieldline/exceptions.hpp>
#include <fieldline/core/magneticField.hpp>
#include <fieldline/core/point.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include <string.h>


namespace Fieldline {
    namespace axiSymmetric {
        /*! \brief Class representing an axial symmetric magnetic equilibrium.
         *
         *  This class represents an axial symmetric magnetic equilibrium.
         *  The equilibrium information is stored as a rectangular grid of the poloidal magnetic flux \f$\psi\f$.
         *  In addition the toroidal magnetic field \f$B_{tor}\f$ and the major radius of the magnetic axis \f$R_0\f$ is stored.
         */
        class magneticField {
            private: 
                const double m_a1[3] = {3.0,  -1.5,    1.0/3.0};
                const double m_a2[3] = {-2.5,   2.0,   -0.5};
                const double m_a3[3] = {0.5,  -0.5,    1.0/6.0};
            public:
                /*! \brief Default constructor
                 *
                 *  The default constructor leaves the class empty.
                 */
                magneticField() {}
                /*! \brief Constructor
                 *
                 *  The constructor initializes the poloidal magnetic flux matrix and the corresponding information, 
                 *  such as spatial extent and the toroidal magnetic field \f$B_{tor}\f$ and the major radius of the magnetic axis
                 *  \f$R_0\f$.
                 */
                magneticField(const double Btor, const double R0, const double Rmin, const double Rmax, 
                    const double zmin, const double zmax, const boost::numeric::ublas::matrix<double> psi) :
                        m_Btor(Btor), m_R0(R0), m_Rmin(Rmin), m_Rmax(Rmax), m_zmin(zmin), m_zmax(zmax),
                        m_dR((Rmax-Rmin)/(psi.size1()-1)), m_dz((zmax-zmin)/(psi.size2()-1)), m_psi(psi) {
                }
                /*! \brief Python constructor
                 *
                 *  The constructor initializes the poloidal magnetic flux matrix and the corresponding information, 
                 *  such as spatial extent and the toroidal magnetic field \f$B_{tor}\f$ and the major radius of the magnetic axis
                 *  \f$R_0\f$.
                 *  This constructor is intended as python interface.
                 *  Do not use this constructor from within C++.
                 */
                magneticField(const double Btor, const double R0, const double Rmin, const double Rmax,
                    const double zmin, const double zmax, const uint32_t NR, const uint32_t Nz, const boost::python::list & psi):
                        m_Btor(Btor), m_R0(R0), m_Rmin(Rmin), m_Rmax(Rmax), m_zmin(zmin), m_zmax(zmax),
                        m_dR((Rmax-Rmin)/(NR-1)), m_dz((zmax-zmin)/(Nz-1)), m_psi(NR,Nz) {
                    uint32_t i = 0;
                    for(uint32_t x = 0; x < NR; ++x) {
                        for(uint32_t y = 0; y < Nz; ++y, ++i) {
                            m_psi(x,y) = boost::python::extract<double>(psi[i]);
                        }
                    }           
                }

                /*! \brief Constructor
                 *
                 *  The constructor reads the equilibrium from file.
                 *  In addition the constructor requires the information about the toroidal magnetic field \f$B_{tor}\f$ and the 
                 *  major radius of the magnetic axis \f$R_0\f$
                 */
                magneticField(const double Btor, const double R0, const std::string & filename) : m_Btor(Btor), m_R0(R0) {
                    std::fstream file(filename.c_str(), std::ios::in);
                    load_from_file(file);
                }
                /*! \brief Constructor
                 *
                 *  The constructor reads the equilibrium from a file instance.
                 *  In addition the constructor requires the information about the toroidal magnetic field \f$B_{tor}\f$ and the 
                 *  major radius of the magnetic axis \f$R_0\f$
                 */
                magneticField(const double Btor, const double R0, std::fstream & file) : m_Btor(Btor), m_R0(R0) {
                    load_from_file(file);
                }
                /*! \brief Copy constructor
                 *
                 *  Copies the equilibrium from the given instance to the new instance.
                 */
                magneticField(const magneticField & rhs) : m_Btor(rhs.m_Btor), m_R0(rhs.m_R0), m_NR(rhs.m_NR), m_Nz(rhs.m_Nz),
                    m_Rmin(rhs.m_Rmin), m_Rmax(rhs.m_Rmax), m_zmin(rhs.m_zmin), m_zmax(rhs.m_zmax), m_dR(rhs.m_dR), m_dz(rhs.m_dz), m_psi(rhs.m_psi) {
                }

                /*! \brief Destructor */
                virtual ~magneticField() {}

                /*! \brief Get the poloidal magnetic flux.
                 *
                 *  This function returns the interpolates magnetic flux at the specified position R,z at the toroidal angle \f$\psi\f$.
                 *  The toroidal angle is a dummy in this class, but is specified for further non axial symmetric cases derived from this class.
                 *  The function returns the poloidal magnetic flux \f$\psi\f$ together with its derivatives \f$\frac{\partial \psi}{\partial R}\f$ and
                 *  \f$\frac{\partial \psi}{\partial z}\f$
                 */
                Fieldline::core::magneticFlux get_magnetic_flux(const double R, const double z, const double phi) const {
                    if(R < m_Rmin || R > m_Rmax || z < m_zmin || z > m_zmax) {
                        throw Fieldline::exceptions::undefinedMagneticField();
                    }
                    double x = (R - m_Rmin)/m_dR;
                    double y = (z - m_zmin)/m_dz;
                    uint32_t jr = std::min((uint32_t)std::max((int32_t)0, (int32_t)x-1), m_NR-4);
                    uint32_t jz = std::min((uint32_t)std::max((int32_t)0, (int32_t)y-1), m_Nz-4);
                    x -= jr;
                    y -= jz;
                    double PF[3];
                    double x2 = x*x;
                    double y2 = y*y;
                    double x3 = x2*x;
                    double y3 = y2*y;
                    double F[4];
                    double dzdy[4];
                    double b1(0.0),b2(0.0),b3(0.0);
                    for(uint32_t i = 0; i < 4; ++i) {
                        b1 = b2 = b3 = 0.0;
                        for(uint32_t j = 0; j < 3; ++j) {
                            b1 += m_a1[j]*m_psi(i+jr,j+jz+1)-m_a1[j]*m_psi(i+jr,jz);
                            b2 += m_a2[j]*(m_psi(i+jr,j+jz+1)-m_psi(i+jr,jz));
                            b3 += m_a3[j]*(m_psi(i+jr,j+jz+1)-m_psi(i+jr,jz));
                        }
                        F[i] = m_psi(jr+i,jz) + b1*y + b2*y2 + b3*y3;
                        dzdy[i] = b1 + 2.0*b2*y + 3.0*b3*y2;
                    }
                    b1 = b2 = b3 = 0.0;
                    for(uint32_t i = 0; i < 3; ++i) {
                        b1 += m_a1[i]*(F[i+1]-F[0]);
                        b2 += m_a2[i]*(F[i+1]-F[0]);
                        b3 += m_a3[i]*(F[i+1]-F[0]);
                    }
                    PF[0] = F[0] + b1*x + b2*x2 + b3*x3;
                    PF[1] = b1 + 2.0*b2*x + 3.0*b3*x2;
                    b1 = b2 = b3 = 0.0;
                    for(uint32_t i = 0; i < 3; ++i) {
                        b1 += m_a1[i]*(dzdy[i+1]-dzdy[0]);
                        b2 += m_a2[i]*(dzdy[i+1]-dzdy[0]);
                        b3 += m_a3[i]*(dzdy[i+1]-dzdy[0]);
                    }
                    PF[2] = dzdy[0] + b1*x + b2*x2 + b3*x3;
                    return Fieldline::core::magneticFlux(PF[0], PF[2], PF[1]);
                }

                /*! \brief Get singularity
                 *
                 *  This function returns the singular point (x- or o-Point) closest to the specified R,z position.
                 */
                Fieldline::core::point get_singularity(const double R, const double z, const double phi) const {
                    if(R < m_Rmin || R > m_Rmax || z < m_zmin || z > m_zmax) {
                        throw Fieldline::exceptions::undefinedMagneticField();
                    }
                    double x = (R - m_Rmin)/m_dR;
                    double y = (z - m_zmin)/m_dz;
                    uint32_t jr = std::min((uint32_t)std::max((int32_t)0, (int32_t)x-1), m_NR-4);
                    uint32_t jz = std::min((uint32_t)std::max((int32_t)0, (int32_t)y-1), m_Nz-4);
                    x -= jr;
                    y -= jz;
                    double PF[3];
                    double x2 = x*x;
                    double y2 = y*y;
                    double x3 = x2*x;
                    double y3 = y2*y;
                    double F[4];
                    double dFdy1[4];
                    double dFdy2[4];
                    double dFdR1;
                    double dFdR2;
                    double dFdz1;
                    double dFdz2;
                    double dFdRdz;
                    double D;
                    double DDDR;
                    double DDDz;
                    double DL;
                    double deltaX;
                    double deltaY;
                    double b1(0.0),b2(0.0),b3(0.0);
                    for(uint32_t i = 0; i < 4; ++i) {
                        b1 = b2 = b3 = 0.0;
                        for(uint32_t j = 0; j < 3; ++j) {
                            b1 += m_a1[j]*m_psi(i+jr,j+jz+1)-m_a1[j]*m_psi(i+jr,jz);
                            b2 += m_a2[j]*(m_psi(i+jr,j+jz+1)-m_psi(i+jr,jz));
                            b3 += m_a3[j]*(m_psi(i+jr,j+jz+1)-m_psi(i+jr,jz));
                        }
                        F[i]     = m_psi(jr+i,jz) + b1*y +     b2*y2 +     b3*y3;
                        dFdy1[i] =                  b1   + 2.0*b2*y  + 3.0*b3*y2;
                        dFdy2[i] =                         2.0*b2    + 6.0*b3*y;
                    }
                    b1 = b2 = b3 = 0.0;
                    for(uint32_t i = 0; i < 3; ++i) {
                        b1 += m_a1[i]*(F[i+1]-F[0]);
                        b2 += m_a2[i]*(F[i+1]-F[0]);
                        b3 += m_a3[i]*(F[i+1]-F[0]);
                    }
                    dFdR1 = b1 + 2.0*b2*x + 3.0*b3*x2;
                    dFdR2 =      2.0*b2   + 6.0*b3*x;
                    b1 = b2 = b3 = 0.0;
                    for(uint32_t i = 0; i < 3; ++i) {
                        b1 += m_a1[i]*(dFdy1[i+1]-dFdy1[0]);
                        b2 += m_a2[i]*(dFdy1[i+1]-dFdy1[0]);
                        b3 += m_a3[i]*(dFdy1[i+1]-dFdy1[0]);
                    }
                    dFdz1  = dFdy1[0] + b1*x +     b2*x2 +     b3*x3;
                    dFdRdz =            b1   + 2.0*b2*x  + 3.0*b3*x2;
                    b1 = b2 = b3 = 0.0;
                    for(uint32_t i = 0; i < 3; ++i) {
                        b1 += m_a1[i]*(dFdy2[i+1]-dFdy2[0]);
                        b2 += m_a2[i]*(dFdy2[i+1]-dFdy2[0]);
                        b3 += m_a3[i]*(dFdy2[i+1]-dFdy2[0]);
                    }
                    dFdz2 = dFdy2[0] + b1 *x + b2*x2 + b3*x2;
                    D = 0.5*(pow(dFdR1,2) + pow(dFdz1,2));
                    DDDR = dFdz1*dFdRdz + dFdR1*dFdR2;
                    DDDz = dFdR1*dFdRdz + dFdz1*dFdz2;
                    DL = -D/(pow(DDDR,2)+pow(DDDz,2));
                    deltaX = DL*DDDR;
                    deltaY = DL*DDDz;
                    if(pow((deltaX*m_dR),2)+pow((deltaY*m_dz),2) < 1.0e-17) {
                        return Fieldline::core::point(R+deltaX*m_dR, z+deltaY*m_dz, true);
                    }
                    return this->get_singularity(R+deltaX*m_dR, z+deltaY*m_dz, 0.0);
                }

                /*! \brief Get magnetic field
                 *
                 *  This function returns the magnetic field at the position R,z at the toroidal angle \f$\phi\f$.
                 *  \f$ \mathbf{B} = (B_R, B_z, B_{tor})\f$ 
                 */
                Fieldline::core::magneticField get_magnetic_field(const double R, const double z, const double phi) const {
                    Fieldline::core::magneticFlux temp = get_magnetic_flux(R,z,phi);
                    return Fieldline::core::magneticField(-(temp.dpsi_dR/m_dz/(2.0*M_PI*R)), temp.dpsi_dz/m_dR/(2.0*M_PI*R), m_R0/R*m_Btor);
                }

                /*! \brief Write matrix to ASCII file
                 *
                 *  This function writes the poloidal flux matrix to an ASCII file.
                 *  This is helpful if one wants to check if the matrix is read in correcly.
                 */
                void write_ASCII_matrix(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    for(uint32_t i = 0; i < m_NR; ++i){ 
                        for(uint32_t j = 0; j < m_Nz; ++j) {
                            file << m_psi(i,j) << '\t';
                        }
                        file << '\n';
                    }
                    file.close();
                }

                /*! \brief Get the position of the magnetic axis
                 *
                 * This function returns the position of the magnetic axis using get_singularity.
                 */
                Fieldline::core::point get_magnetic_axis(const double phi = 0.0) const {
                    return this->get_singularity(m_Rmax/2.0 + m_Rmin/2.0, 0.0, phi);
                }

                /*! \brief Get the position of the lower x-point.
                 *
                 * This function returns the position of the lower x-point using get_singularity.
                 */
                Fieldline::core::point get_lower_x_point(const double phi = 0.0) const {
                    return this->get_singularity(m_Rmax*0.6 + m_Rmin*0.4, m_zmin*0.8, phi);
                }

                /*! \brief Get the position of the upper x-point.
                 *
                 * This function returns the position of the upper x-point using get_singularity.
                 */
                Fieldline::core::point get_upper_x_point(const double phi = 0.0) const {
                    return this->get_singularity(m_Rmax*0.6 + m_Rmin*0.4, m_zmax*0.8, phi);
                }


            protected:
                double m_Btor;      /*!< \brief Toroidal magnetic field \f$B_{tor}\f$. */
                double m_R0;        /*!< \brief Major radius of magnetic axis \f$R_0\f$. */
                uint32_t m_NR;      /*!< \brief Number of poloidal flux points in R direction. */
                uint32_t m_Nz;      /*!< \brief Number of poloidal flux points in z direction. */
                double m_Rmin;      /*!< \brief Minimum major radius \f$R_{min}\f$ of the poloidal flux matrix. */
                double m_Rmax;      /*!< \brief Maximum major radius \f$R_{max}\f$ of the poloidal flux matrix. */
                double m_zmin;      /*!< \brief Minimum vertical position \f$z_{min}\f$ of the poloidal flux matrix. */
                double m_zmax;      /*!< \brief Maximum vertical position \f$z_{max}\f$ of the poloidal flux matrix. */
                double m_dR;        /*!< \brief Spatial distance \f$\Delta R\f$ between points of the poloidal flux matrix in R direction. */
                double m_dz;        /*!< \brief Spatial distance \f$\Delta z\f$ between points of the poloidal flux matrix in z direction. */
                boost::numeric::ublas::matrix<double> m_psi;    /*!< \brief Poloidal flux matrix \f$\psi\f$. */

                /*! \brief Load poloidal flux matrix from file.
                 *
                 *  This file reads the poloidal flux matrix from file.
                 *  The file format is that of the FORTRAN code.
                 */
                void load_from_file(std::fstream & file) {
                    std::string temp;
                    file >> m_NR >> m_Nz;
                    m_psi = boost::numeric::ublas::matrix<double>(m_NR, m_Nz);
                    std::getline(file, temp);
                    file >> m_Rmin >> m_Rmax >> m_zmin >> m_zmax;
                    m_dR = (m_Rmax - m_Rmin)/(m_NR-1);
                    m_dz = (m_zmax - m_zmin)/(m_Nz-1);
                    for(uint32_t j = 0; j < m_Nz; ++j) {
                        for(uint32_t i = 0; i < m_NR; ++i) {
                            file >> m_psi(i,j);
                        }
                    }
                }

        };
    }
}


#endif 
