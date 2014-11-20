#ifndef include_fieldline_core_magneticFlux_hpp
#define include_fieldline_core_magneticFlux_hpp

namespace Fieldline {
    namespace core {
        class magneticFlux {
            public:
                magneticFlux() : psi(0.0), dpsi_dR(0.0), dpsi_dz(0.0) {}
                magneticFlux(const double Psi, const double dPsi_dR, const double dPsi_dz) : psi(Psi), dpsi_dR(dPsi_dR), dpsi_dz(dPsi_dz) {}
                magneticFlux(const magneticFlux & rhs) : psi(rhs.psi), dpsi_dR(rhs.dpsi_dR), dpsi_dz(rhs.dpsi_dz) {}
                magneticFlux & operator= (const magneticFlux & rhs) {
                    if(this != &rhs) {
                        psi = rhs.psi;
                        dpsi_dR = rhs.dpsi_dR;
                        dpsi_dz = rhs.dpsi_dz;
                    }
                    return *this;
                }

                double psi;
                double dpsi_dR;
                double dpsi_dz;
        };
    }
}


#endif
