#ifndef include_fieldline_core_point_hpp
#define include_fieldline_core_point_hpp

#include <iostream>
#include <math.h>


namespace Fieldline {
    namespace core {
        /*! \brief Point class
         *
         * This class contains the information about a point in the poloidal plain \f$(R,z)\f$.
         */
        class point {
            public:
                /*! \brief Default constructor. */
                point() : R(0.0), z(0.0), hit(false) {}
                /*! /brief Constructor.
                 *
                 *  Initialized the point with the given \f$(R,z)\f$ coordinate and the information if it got hit.
                 */
                point(const double RIn, const double zIn, const bool hitIn) : R(RIn), z(zIn), hit(hitIn) {}
                /*! \brief Copy constructor. */
                point(const point & rhs) : R(rhs.R), z(rhs.z), hit(rhs.hit) {}
                /*! \brief Destructor. */
                virtual ~point() {}
                /*! \brief Assignment operator
                 *
                 *  This operator copies the information of the point into the given instance.
                 */
                point & operator= (const point & rhs) {
                    if(this != &rhs) {
                        R = rhs.R;
                        z = rhs.z;
                        hit = rhs.hit;
                    }
                    return *this;
                }
                double R; /*!< \brief R value of the point. */
                double z; /*!< \brief z value of the point. */
                bool hit; /*!< \brief Information if point was hit. */
                /*! \brief Distance between two points.
                 *
                 *  This function calculates the distance in the poloidal plain between two points.
                 *  \f$ d = \sqrt{(R_1 - R_0)^2 + (z_1 - z_0)^2} \f$
                */
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
