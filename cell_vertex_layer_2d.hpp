#pragma once

#include "nodes.hpp"
#include "object.hpp"
#include <vector>
#include <utility>
#include <cstdint>
#include <array>
#include <iostream>
#include <functional>
#include <memory>
#include <span>

class cell_vertex_layer_2d {
    private:
        using u32 = std::uint32_t;
        using f32 = float;

        struct _cluster {
            u32 begin, size;
            f32 x, y;
        };

        std::shared_ptr<std::vector<nodeCV>> m_nodes;
        f32 m_delta;
        u32 m_layer;
        std::vector<_cluster> m_clusters;
    public:

        class cluster {
            private:
                u32 m_begin;
                u32 m_layer;
                f32 m_x, m_y;
                f32 m_delta;
                std::span<const nodeCV> m_nodes;

                inline cluster(const cell_vertex_layer_2d& layer, u32 cluIdx) : 
                    m_begin(layer.m_clusters[cluIdx].begin),
                    m_layer(layer.m_layer),
                    m_x(layer.m_clusters[cluIdx].x),
                    m_y(layer.m_clusters[cluIdx].y),
                    m_delta(layer.m_delta),
                    m_nodes(std::span<const nodeCV>(*layer.m_nodes).subspan(m_begin, layer.m_clusters[cluIdx].size)) {}

                friend class cell_vertex_layer_2d;

                void traverse(std::function<void(long, long, long)> f) const;
                    
            public:
                cluster() = default;
                inline u32 get_layer() const { return m_layer; }
                inline std::pair<f32, f32> get_pos_bottom_left() const { return std::make_pair(m_x, m_y); }
                inline void get_pos_bottom_left(std::span<f32, 2> point) const { point[0] = m_x; point[1] = m_y; }
                inline f32 get_delta() const { return m_delta; }
                inline f32 get_pos_x(long xoff) const { return m_x + m_delta * static_cast<f32>(xoff); }
                inline f32 get_pos_y(long yoff) const { return m_y + m_delta * static_cast<f32>(yoff); }

                void to_point2d(std::ostream& os) const;
                void to_dense(std::vector<bool>& denseGrid, u32 nx, u32 ny, u32 x0, u32 y0) const;
                void obstacle_to_dense(std::vector<bool>& denseGrid, u32 nx, u32 ny, u32 x0, u32 y0) const;
        };
        
        cell_vertex_layer_2d() = default;
        cell_vertex_layer_2d(float xmin, float ymin, float xmax, float ymax, float delta, const object<2>& obj);

        cell_vertex_layer_2d refine(const object<2>& obj);
        inline u32 get_cluster_size() const { return m_clusters.size(); }

        inline bool has_left(u32 idx) const { return m_nodes->at(idx).has_left(); }
        inline bool has_right(u32 idx) const { return m_nodes->at(idx).has_right(); } 
        inline bool has_up(u32 idx) const { return m_nodes->at(idx).has_up(); }
        inline bool has_down(u32 idx) const { return m_nodes->at(idx).has_down(); } 

        inline u32 get_left(u32 idx) const { return m_nodes->at(idx).has_left() ? idx - 1 : static_cast<u32>(-1); }
        inline u32 get_right(u32 idx) const { return m_nodes->at(idx).has_right() ? idx + 1 : static_cast<u32>(-1); }
        inline u32 get_up(u32 idx) const { return m_nodes->at(idx).up; }
        inline u32 get_down(u32 idx) const { return m_nodes->at(idx).down; } 

        inline cluster get_cluster(u32 idx) const { return cluster(*this, idx); }

};
