#pragma once

#include "types.hpp"

template<std::uint32_t N>
class object;

template<>
class object<2> {
    public:
        virtual bool in(float x, float y) const = 0;
        virtual bool in(cpoint<2> p) const = 0;
        virtual float distance(float x, float y) const = 0;
        virtual float distance(cpoint<2> p) const = 0;
};

template<>
class object<3> {
    public:
        virtual bool in(float x, float y, float z) const = 0;
        virtual bool in(cpoint<3> p) const = 0;
        virtual float distance(float x, float y, float z) const = 0;
        virtual float distance(cpoint<3> p) const = 0;
};