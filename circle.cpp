#include "circle.hpp"
#include <cmath>

float circle::distance(float x, float y) const {
    return std::sqrt(std::pow(x - m_x, 2) + std::pow(y - m_y, 2)) - m_radius;
}

float circle::distance(cpoint<2> p) const {
    return distance(p[0], p[1]);
}

bool circle::in(float x, float y) const {
    return distance(x, y) <= 0;
}

bool circle::in(cpoint<2> p) const {
    return in(p[0], p[1]);
}

float circle2::distance(float x, float y) const {
    return std::min(c1.distance(x, y), c2.distance(x, y));
}
float circle2::distance(cpoint<2> p) const {
    return distance(p[0], p[1]);
}
bool circle2::in(float x, float y) const {
    return distance(x, y) <= 0;
}
bool circle2::in(cpoint<2> p) const {
    return in(p[0], p[1]);
}

