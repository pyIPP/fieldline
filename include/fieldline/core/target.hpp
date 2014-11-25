#ifndef include_fieldline_core_target_hpp
#define include_fieldline_core_target_hpp

#include <vector>
#include <fstream>
#include <algorithm>
#include <fieldline/core/point.hpp>
#include <fieldline/core/line.hpp>

#ifdef WITH_PYTHON
    #include <boost/python/list.hpp>
#endif 


namespace Fieldline {
    namespace core {
        /*! \brief Class containing the target contour.
         *
         *  The contour of the target is represented in the poloidal plane as a list of R,z values.
         */
        class target {
            public:
                /*! \brief Constructor
                 *
                 *  This constructor sets the target according to the given R,z values.
                 */
                target(const std::vector<double> & R, const std::vector<double> & z) : m_R(R), m_z(z) {}
                /*! \brief Constructor
                 *
                 *  This constructor loads the target contour from the specified location.
                 */
                target(const std::string & filename) : m_R(), m_z() {
                    load_from_file(filename);
                }
                /*! \brief Copy constructor */
                target(const target & rhs) : m_R(rhs.m_R), m_z(rhs.m_z) {}
#ifdef WITH_PYTHON
                /*! \brief Python constructor
                 *
                 *  This constructor sets the target according to the given R,z values.
                 *  It is only intended to interface with python, do not use this constructor from within C++.
                 */
                target(const boost::python::list & R, const boost::python::list & z) : m_R(), m_z() {
                    uint32_t N = boost::python::len(R);
                    m_R.resize(N);
                    m_z.resize(N);
                    for(uint32_t i = 0; i < N; ++i) {
                        m_R[i] = boost::python::extract<double>(R[i]);
                        m_z[i] = boost::python::extract<double>(z[i]);
                    }
                }
#endif                
                /*! \brief Destructor */
                virtual ~target() {}
                /*! \brief Assignment operator
                 *
                 *  Copies the target contour from the given instance.
                 */
                target & operator=(const target & rhs) {
                    if(this != &rhs) {
                        m_R = rhs.m_R;
                        m_z = rhs.m_z;
                    }
                    return *this;
                }
                /*! \brief Target size.
                 *
                 *  Returns the number of points forming the target contour.
                 */
                uint32_t size() const { return m_R.size(); }
                /*! \brief Save target contour to file.
                 *
                 *  This function saves the target contour to an ASCII file.
                 */
                void save_to_file(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    auto R = m_R.cbegin();
                    auto z = m_z.cbegin();
                    file << m_R.size() << '\n';
                    for(; R != m_R.cend(); ++R, ++z) {
                        file << *R << '\t' << *z << '\n';
                    }
                    file.close();
                }

                /*! \brief Load target contour from file.
                 *
                 *  This functions load the target contour from the specified ASCII file.
                */
                void load_from_file(const std::string & filename) {
                    std::fstream file(filename.c_str(), std::ios::in);
                    uint32_t N;
                    file >> N;
                    m_R.resize(N);
                    m_z.resize(N);
                    auto R = m_R.begin();
                    auto z = m_z.begin();
                    for( ; R != m_R.end(); ++R, ++z) {
                        file >> *R >> *z;
                    }
                    file.close();
                }

                /*! \brief Get field line intersection with the target.
                 *
                 *  This function calculates the intersection of the field line with the target contour.
                 *  If no intersection is found the origin is returned with hit set to false.
                 *  If two or more intersections are found the one closest to the beginning of the given line is returned.
                 */
                Fieldline::core::point get_intersection(const Fieldline::core::line & line) const {
                    Fieldline::core::point temp;
                    std::vector<Fieldline::core::point> points;
                    for(uint32_t i = 0; i < m_R.size()-1; ++i) {
                        temp = line.get_intersection(Fieldline::core::line(m_R[i], m_z[i], m_R[i+1], m_z[i+1]));
                        if(temp.hit) {
                            points.push_back(temp);
                        }
                    }
                    switch(points.size()) {
                        case 0:
                            return Fieldline::core::point(0,0,false);
                        case 1:
                            return points[0];
                        default:
                            temp.R = line.get_R0();
                            temp.z = line.get_z0();
                            std::vector<double> distance;
                            for(auto iter = points.cbegin(); iter != points.cend(); ++iter) {
                                distance.push_back(iter->get_distance(temp));
                            }
                            uint32_t i = std::distance(distance.cbegin(), std::min_element(distance.cbegin(), distance.cend()));
                            return points[i];
                    }
                }

#ifdef WITH_PYTHON
                /*! \brief Return R values as python list.
                 *
                 *  This function returns the R values of the target contour as a python list.
                 *  This function is only intended as an interface to python, do not use this function from within C++.
                 */
                boost::python::list get_R_python() const {
                    boost::python::list output;
                    for(auto iter = m_R.cbegin(); iter != m_R.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

                /*! \brief Return z values as python list.
                 *
                 *  This function returns the z values of the target contour as a python list.
                 *  This function is only intended as an interface to python, do not use this function from within C++.
                 */
                boost::python::list get_z_python() const {
                    boost::python::list output;
                    for(auto iter = m_z.cbegin(); iter != m_z.cend(); ++iter) {
                        output.append(*iter);
                    }
                    return output;
                }

#endif                

            protected:
                std::vector<double> m_R; /*!< \brief R values of the target contour. */
                std::vector<double> m_z; /*!< \brief z values of the target contour. */
        };
    }
}


#endif
