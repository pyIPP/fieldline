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
        class magneticField {
            private: 
                const double m_a1[3] = {3.0,  -1.5,    1.0/3.0};
                const double m_a2[3] = {-2.5,   2.0,   -0.5};
                const double m_a3[3] = {0.5,  -0.5,    1.0/6.0};
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

               Fieldline::core::point get_secondDerivation(const double x, const double y, const double jr, const double jz) const {
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
                    bool exist;
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

                    exist = pow((deltaX*m_dR),2)+pow((deltaY*m_dz),2) < 1.0e-17;
                    return Fieldline::core::point(deltaX, deltaY, exist);
                }


                Fieldline::core::magneticField get_magnetic_field(const double R, const double z, const double phi) const {
                    Fieldline::core::magneticFlux temp = get_magnetic_flux(R,z,phi);
                    return Fieldline::core::magneticField(-(temp.dpsi_dR/m_dz/(2.0*M_PI*R)), temp.dpsi_dz/m_dR/(2.0*M_PI*R), m_R0/R*m_Btor);
                }

                Fieldline::core::point get_singularity(double R, double z) {
                    if(R < m_Rmin || R > m_Rmax || z < m_zmin || z > m_zmax) {
                        throw Fieldline::exceptions::undefinedMagneticField();
                    }
                    double x = (R - m_Rmin)/m_dR;
                    double y = (z - m_zmin)/m_dz;
                    uint32_t jr = std::min((uint32_t)std::max((int32_t)0, (int32_t)x-1), m_NR-4);
                    uint32_t jz = std::min((uint32_t)std::max((int32_t)0, (int32_t)y-1), m_Nz-4);
                    x -= jr;
                    y -= jz;
                    bool stop=false;
                    double f;
                    double fY;
                    uint32_t tempCase = 0;
                    double deltaX;
                    double deltaY;
                    bool exist;
                    Fieldline::core::point tempPoint;
                    while(!stop){
                        tempPoint  = get_secondDerivation(x,y,jr,jz);
                        deltaX = tempPoint.R;
                        deltaY = tempPoint.z;
                        exist = tempPoint.hit;
                        if(exist){
                            R = (x+deltaX+jr)*m_dR+m_Rmin;
                            z = (y+deltaY+jz)*m_dz+m_zmin;
                            return Fieldline::core::point(R,z,true);
                        }
			            tempCase = 0;
                        if(x+deltaX<1.0){
                            f = (1.0-x)/deltaX;
                            tempCase = 1;
                        }
                        else if(x+deltaX>2.0){
                            f = (2.0-x)/deltaX;
                            tempCase = 2;
                        }
                        if(y+deltaY<1.0){
                            fY = (1.0 - y)/deltaY;
                            if(tempCase==0 || fY<f){
                                f = fY;
                                tempCase = 3;
                            }
                        }
                        else if(y+deltaY>2.0){
                            fY = (2.0 - y)/deltaY;
                            if(tempCase==0 || fY<f){
                                f = fY;
                                tempCase = 4;
                            }
                        }
                        switch(tempCase) {
                            case 0: x = x+deltaX; y = y+deltaY; break;
                            case 1: jr--; x = 2.0; y += f*deltaY;break;
                            case 2: jr++; x = 1.0; y += f*deltaY;break;
                            case 3: jz--; y = 2.0; x += f*deltaX;break;
                            case 4: jz++; y = 1.0; x += f*deltaX;break;
                        }
                        if(jz>m_Nz-3 || jz<1 || jr>m_NR-3 || jr<1){
                            stop = true;
                        }
                    }
                    return Fieldline::core::point(R,z,exist);
                }


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
