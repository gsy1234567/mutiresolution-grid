#pragma once

#include <span>
#include <cstdint>

template<std::uint32_t N>
using cpoint = std::span<const float, N>;

template<std::uint32_t N>
using point = std::span<float, N>;
