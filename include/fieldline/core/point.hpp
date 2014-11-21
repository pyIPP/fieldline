#ifndef include_fieldline_core_point_hpp
#define include_fieldline_core_point_hpp

#include <iostream>
#include <math.h>

namespace Fieldline {
    namespace core {
        class point {
            public:
                point() : R(0.0), z(0.0), hit(false) {}
                point(const double RIn, const double zIn, const bool hitIn) : R(RIn), z(zIn), hit(hitIn) {}
                point(const point & rhs) : R(rhs.R), z(rhs.z), hit(rhs.hit) {}
                virtual ~point() {}
                point & operator= (const point & rhs) {
                    if(this != &rhs) {
                        R = rhs.R;
                        z = rhs.z;
                        hit = rhs.hit;
                    }
                    return *this;
                }


                double R;
                double z;
                bool hit;
                double get_distance(const point & rhs) const {
                    return sqrt(pow(R-rhs.R,2) + pow(z-rhs.z,2));
                }
        };

        std::ostream & operator<< (std::ostream & out, const point & rhs) {
            out << rhs.R << '\t' << rhs.z << '\t' << rhs.hit << '\n';
            return out;
        }
    }
}


#endif 
