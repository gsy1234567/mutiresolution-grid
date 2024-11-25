#include "cell_vertex_layer_2d.hpp"
#include "circle.hpp"
#include <iostream>
#include <format>

void print_layer(const cell_vertex_layer_2d& layer, std::uint32_t nx = 75, std::uint32_t ny = 75, std::uint32_t x0 = 20, std::uint32_t y0 = 20) {
    for (auto i = 0U; i < layer.get_cluster_size(); ++i) {
        auto cluster = layer.get_cluster(i);
        std::cout << std::format("Cluster {}:\n", i);
        std::vector<bool> denseGrid;
        cluster.to_dense(denseGrid, nx, ny, x0, y0);
        for(auto y = ny ; y != 0U ; --y) {
            for(auto x = 0U ; x < nx ; ++x) {
                std::cout << (denseGrid[(y-1) * nx + x] ? "██" : "  ");
            }
            std::cout << std::endl;
        }
    }
}

void print_layer_obstacle(const cell_vertex_layer_2d& layer, std::uint32_t nx = 75, std::uint32_t ny = 75, std::uint32_t x0 = 20, std::uint32_t y0 = 20) {
    for (auto i = 0U; i < layer.get_cluster_size(); ++i) {
        auto cluster = layer.get_cluster(i);
        std::cout << std::format("Cluster {}:\n", i);
        std::vector<bool> denseGrid;
        cluster.obstacle_to_dense(denseGrid, nx, ny, x0, y0);
        for(auto y = ny ; y != 0U ; --y) {
            for(auto x = 0U ; x < nx ; ++x) {
                std::cout << (denseGrid[(y-1) * nx + x] ? "▓▓" : "  ");
            }
            std::cout << std::endl;
        }
    }
}



void dump_layer(std::ostream& os, const cell_vertex_layer_2d& layer) {
    for(auto i = 0U ; i < layer.get_cluster_size() ; ++i) {
        os << std::format("cluster number: {}\n", i);
        layer.get_cluster(i).to_point2d(os);
        os << "-----------------------------------------\n";

    }
}

int main() {
    auto c = circle(0.4f, 0.4f, 0.24f);
    cell_vertex_layer_2d layer0 (0.f, 0.f, 1.1f, 1.1f, 0.1f, c);
    print_layer(layer0);
    print_layer_obstacle(layer0);
    auto layer1 = layer0.refine(c);
    print_layer(layer1);
    print_layer_obstacle(layer1);
    auto layer2 = layer1.refine(c);
    print_layer(layer2);
    print_layer_obstacle(layer2);
    auto layer3 = layer2.refine(c);
    print_layer(layer3);
    print_layer_obstacle(layer3);
    dump_layer(std::cout, layer3);


    return 0;
}
