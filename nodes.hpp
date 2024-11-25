#pragma once

#include <cstdint>

class node {
    public:
        std::uint32_t up, down;
        std::uint32_t parent, child;
        constexpr bool is_obstacle() const { return (args & obstacleMask) != 0; }
        inline void set_obstacle(bool exist) { exist ? args |= obstacleMask : args &= ~obstacleMask; }
        constexpr bool has_left() const { return (args & leftMask) != 0; }
        inline void set_left(bool exist) { exist ? args |= leftMask : args &= ~leftMask; } 
        constexpr bool has_right() const { return (args & rightMask) != 0; }
        inline void set_right(bool exist) { exist ? args |= rightMask : args &= ~rightMask; }
        constexpr bool has_up() const { return up != static_cast<std::uint32_t>(-1); }
        constexpr bool has_down() const { return down != static_cast<std::uint32_t>(-1); }
        constexpr bool has_child() const { return child != static_cast<std::uint32_t>(-1); }
        inline void invalid_child() { child = static_cast<std::uint32_t>(-1); }
        inline std::uint32_t get_up() const { return up; }
        inline std::uint32_t get_down() const { return down; }
    protected:
        std::uint32_t args;
    private:
        static constexpr std::uint32_t obstacleMask = 0x0000'0001;
        static constexpr std::uint32_t leftMask = 0x0000'0002;
        static constexpr std::uint32_t rightMask = 0x0000'0004;
};

//Cell-Vertex Node
class nodeCV : public node {
    public:
        /**
         * \brief Check whether the node is a `regular fine node`.
        */
        constexpr bool is_rfn() const { return (args & rfnMask) != 0; }
        /**
         * \brief If `exist` is true, set the node to be a `regular fine node`.
         *        Otherwise, set the node to be `non-regular coarse node`.
        */
        inline void set_rfn(bool exist) { exist ? args |= rfnMask : args &= ~rfnMask; }
        /**
         * \brief Check whether the node is a `fine interface node with co-located coarse partner`.
        */
        constexpr bool is_fin() const { return (args & finMask) != 0; }
        /**
         * \brief If `exist` is true, set the node to be a `fine interface node`.
         *        Otherwise, set the node to be `non-fine interface node`.
        */
        inline void set_fin(bool exist) {exist ? args |= finMask : args &= ~finMask; }
        /**
         * \brief Check whether the node is a `hanging middle node`.
        */
        constexpr bool is_hmn() const { return (args & hmnMask) != 0; }
        /**
         * \brief If `exist` is true, set the node to be a `hanging middle node`.
         *        Otherwise, set the node to be `non-hanging middle node`.
        */
        inline void set_hmn(bool exist) { exist ? args |= hmnMask : args &= ~hmnMask; }
        /**
         * \brief Check whether the node is a `regular coarse node`.
        */
        constexpr bool is_rcn() const { return (args & rcnMask) != 0; }
        /**
         * \brief If `exist` is true, set the node to be a `regular coarse node`.
         *        Otherwise, set the node to be `non-regular coarse node`.
        */
        inline void set_rcn(bool exist) { exist ? args |= rcnMask : args &= ~rcnMask; }
        /**
         * \brief Check whether the node is a `coarse interface node`.
        */
        constexpr bool is_cin() const { return (args & cinMask) != 0; }
        /**
         * \brief If `exist` is true, set the node to be a `coarse interface node`.
         *        Otherwise, set the node to be `non-coarse interface node`.
        */
        inline void set_cin(bool exist) { exist ? args |= cinMask : args &= ~cinMask; }
    private:
        //`regular fine node` mask
        static constexpr std::uint32_t rfnMask = 0x0000'0008;
        //`fine interface node with co-located coarse partner` mask
        static constexpr std::uint32_t finMask = 0x0000'0010;
        //`hanging middle node` mask
        static constexpr std::uint32_t hmnMask = 0x0000'0020;
        //`regular coarse node` mask
        static constexpr std::uint32_t rcnMask = 0x0000'0040;
        //`coarse interface node with co-located fine partner` mask
        static constexpr std::uint32_t cinMask = 0x0000'0080;
};