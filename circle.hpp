#pragma once

#include "object.hpp"

class circle : public object<2> {
    private:
        float m_radius;
        float m_x, m_y;
    public:
        circle(float x, float y, float radius) : m_radius(radius), m_x(x), m_y(y) {}

        virtual float distance(float x, float y) const override;
        virtual float distance(cpoint<2> p) const override;
        virtual bool in(float x, float y) const override;
        virtual bool in(cpoint<2> p) const override;

        inline float get_radius() const { return m_radius; }
        inline float get_x() const {  return m_x; }
        inline float get_y() const { return m_y; }
};

class circle2 : public object<2> {
    private:
        circle c1, c2;
    public:
        circle2(circle c1, circle c2) : c1(c1), c2(c2) {}
        circle2(float x1, float y1, float r1, float x2, float y2, float r2) : 
            c1(x1, y1, r1), c2(x2, y2, r2) {}

        virtual float distance(float x, float y) const override;
        virtual float distance(cpoint<2> p) const override;
        virtual bool in(float x, float y) const override;
        virtual bool in(cpoint<2> p) const override;
};