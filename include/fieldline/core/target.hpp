#ifndef include_fieldline_core_target_hpp
#define include_fieldline_core_target_hpp

#include <vector>
#include <fstream>
#include <algorithm>
#include <fieldline/core/point.hpp>
#include <fieldline/core/line.hpp>

namespace Fieldline {
    namespace core {
        class target {
            public:
                target(const std::vector<double> & R, const std::vector<double> & z) : m_R(R), m_z(z) {}
                target(const std::string & filename) : m_R(), m_z() {
                    std::fstream file(filename.c_str(), std::ios::in);
                    load_from_file(file);
                }
                target(std::fstream & file) : m_R(), m_z() {
                    load_from_file(file);
                }
                virtual ~target() {}
                target & operator=(const target & rhs) {
                    if(this != &rhs) {
                        m_R = rhs.m_R;
                        m_z = rhs.m_z;
                    }
                    return *this;
                }
                uint32_t size() const { return m_R.size(); }
                void save_to_file(const std::string & filename) const {
                    std::fstream file(filename.c_str(), std::ios::out);
                    save_to_file(file);
                }
                void save_to_file(std::fstream & file) const {
                    auto R = m_R.cbegin();
                    auto z = m_z.cbegin();
                    file << m_R.size() << '\n';
                    for(; R != m_R.cend(); ++R, ++z) {
                        file << *R << '\t' << *z << '\n';
                    }
                    file.close();
                }

                void load_from_file(std::fstream & file) {
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

                Fieldline::core::point intersection(const double R0, const double z0, const double R1, const double z1) const {
                    return intersection(Fieldline::core::line(R0, z0, R1, z1));
                }

                Fieldline::core::point intersection(const Fieldline::core::line & line) const {
                    Fieldline::core::point temp;
                    std::vector<Fieldline::core::point> points;
                    for(uint32_t i = 0; i < m_R.size()-1; ++i) {
                        temp = line.intersection(Fieldline::core::line(m_R[i], m_z[i], m_R[i+1], m_z[i+1]));
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
                                distance.push_back(iter->distance(temp));
                            }
                            uint32_t i = std::distance(distance.cbegin(), std::min_element(distance.cbegin(), distance.cend()));
                            return points[i];
                    }
                }

            protected:
                std::vector<double> m_R;
                std::vector<double> m_z;
        };
    }
}


#endif
