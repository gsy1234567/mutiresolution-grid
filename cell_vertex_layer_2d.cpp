#include "cell_vertex_layer_2d.hpp"

#include <array>
#include <vector>
#include <memory>
#include <cstdint>
#include <assert.h>
#include <format>
#include <cmath>
#include <list>
#include <tuple>

struct cluster_cell_vertex {
    std::shared_ptr<std::vector<nodeCV>> nodes = std::make_shared<std::vector<nodeCV>>();
    std::uint32_t begin;
    std::uint32_t size;
    std::array<float, 2> pos0; //The position of `nodes[begin]`.
    float delta;               //The distance between two neighbors.
    inline void invalid_children() {
        for(std::uint32_t off = 0 ; off < size ; ++off) {
            (*nodes)[begin + off].invalid_child();
        }
    }
};

struct Range {
    std::uint32_t begin, end;

    inline bool in(std::uint32_t i) const {
        return begin <= i && i < end;
    } 
    
    inline std::uint32_t len() const {
        return end - begin;
    }
};

class RowRanges : public std::vector<Range> {
    public:
        inline friend std::ostream& operator<<(std::ostream& os, const RowRanges& ranges) {
            for(auto [begin, end] : ranges) {
                os << std::format("[{}, {}), ", begin, end);
            }
            return os;
        }

        inline std::uint32_t number_nodes() const {
            std::uint32_t cnt = 0;
            for(auto [begin, end] : *this) {
                cnt += end - begin;
            }
            return cnt;
        }
};

class Ranges2D : public std::vector<RowRanges> {
    public:
        inline friend std::ostream& operator<<(std::ostream& os, const Ranges2D& ranges2d) {
            for(std::uint32_t r = ranges2d.size() ; r != 0 ; --r) {
                os << std::format("{}: ", r-1) << ranges2d[r-1] << std::endl;
            } 
            return os;
        }

    /**
     * \brief Calculate the number of nodes below each row.
     * \param[out] lowerCnt `lowerCnt[r]` is the number of nodes below row `r`.
     * \attention `lowerCnt`[this->size()] is valid, it equal to the number
     * all nodes.
    */
    inline void get_lower(std::vector<std::uint32_t>& lowerCnt) const {
        lowerCnt.resize(size() + 1, 0);
        for(auto r = 0U ; r != size() ; ++r) {
            lowerCnt[r + 1] = (*this)[r].number_nodes();
        }
        for(auto r = 1U ; r != size() + 1 ; ++r) {
            lowerCnt[r] += lowerCnt[r-1];
        }
    }
};

static void get_ranges2d(
    float xmin, float ymin, 
    float xmax, float ymax, 
    float delta, 
    const object<2>& obj, 
    Ranges2D& ranges2d
);

static void get_neighbor(
    const Ranges2D& ranges2d, 
    cluster_cell_vertex& nodes
);

static void get_obstacle(
    const object<2>& obj, 
    cluster_cell_vertex& nodes
);

static void bfs(
    const cluster_cell_vertex& nodes, 
    std::function<void(long, long, std::uint32_t)> option
);

struct CNode {
    std::uint32_t left, right, up, down;
    std::uint32_t parent;
    float x, y;
    static constexpr std::uint32_t invalid = static_cast<std::uint32_t>(-1);

    CNode(std::uint32_t par, float x, float y) : left(invalid), right(invalid), up(invalid), down(invalid), parent(par), x(x), y(y) {}
    CNode(float x, float y) : CNode(invalid, x, y) {}

    inline void make_edge_center_x(CNode& l, CNode& r, std::uint32_t iSelf, std::uint32_t iLeft, std::uint32_t iRight) {
        assert(left == invalid);
        assert(right == invalid);
        assert(l.right == invalid);
        assert(r.left == invalid);

        left = iLeft;
        right = iRight;
        l.right = r.left = iSelf;
    }

    inline void make_edge_center_y(CNode& u, CNode& d, std::uint32_t iSelf, std::uint32_t iUp, std::uint32_t iDown) {
        assert(up == invalid);
        assert(down == invalid);
        assert(u.down == invalid);
        assert(d.up == invalid);

        up = iUp;
        down = iDown;
        u.down = d.up = iSelf;
    }

    inline void make_face_center(
        CNode& l, CNode& r, CNode& u, CNode& d, 
        std::uint32_t iSelf, 
        std::uint32_t iLeft, std::uint32_t iRight, std::uint32_t iUp, std::uint32_t iDown
    ) {
        assert(up == invalid);
        assert(down == invalid);
        assert(left == invalid);
        assert(right == invalid);
        assert(l.right == invalid);
        assert(r.left == invalid);
        assert(u.down == invalid);
        assert(d.up == invalid);

        left = iLeft;
        right = iRight;
        up = iUp;
        down = iDown;
        l.right = r.left = u.down = d.up = iSelf;
    } 

    inline bool has_left() const { return left != invalid; }
    inline bool has_right() const { return right != invalid; }
    inline bool has_up() const { return up != invalid; }
    inline bool has_down() const { return down != invalid; }
};

static void bfs(
    const std::vector<CNode>& cnodes, 
    std::function<void(long, long, std::uint32_t, std::uint32_t)> option
);

static void to_dense(
    const std::vector<CNode>& nodes, 
    std::uint32_t nx, std::uint32_t ny, 
    std::uint32_t x0, std::uint32_t y0, 
    std::vector<std::vector<bool>>& grids
);

static void refine_to_cnode(
    cluster_cell_vertex& curCluster, 
    std::vector<CNode>& cnodes
);

static void cnodes_to_ranges2d(
    const std::vector<CNode>& cnodes, 
    std::uint32_t& prevChildren,
    cluster_cell_vertex& parCluster,
    std::vector<Ranges2D>& ranges2ds, 
    std::vector<std::pair<float, float>>& pos0
);

/**
 * \brief Initialize the first cluster.
 * \param[in] obj The object.
 * \param[out] initCluster The first cluster, this function fill the neighbor
 *  information and obstacle informtion of it.
*/
static void init_first_cluster(
    float xmin, float ymin, float xmax, float ymax, float delta, 
    const object<2>& obj,
    cluster_cell_vertex& initCluster
);

/**
 * \brief Refine the `curCluster` cluster, put the result into `nextCluster`.
 * \param[in, out] curCluster The current cluster, the neighbor information and 
 *  obstacle information must be given. This function fill the children information
 * of it.
 * \param[in, out] nextCluster The next cluster, this function fill the neighbor 
 *  information and obstacle information of it.      
*/
static void refinement(
    const object<2>& obj,
    std::vector<cluster_cell_vertex>& curCluster, 
    std::vector<cluster_cell_vertex>& nextCluster 
);

static void init_first_cluster(
    float xmin, float ymin, float xmax, float ymax, float delta, 
    const object<2>& obj,
    cluster_cell_vertex& initCluster
) {
    Ranges2D ranges2d;
    std::vector<std::uint32_t> lowerCnt;
    get_ranges2d(xmin, ymin, xmax, ymax, delta, obj, ranges2d);

    initCluster.pos0[0] = xmin;
    initCluster.pos0[1] = ymin;
    initCluster.delta = delta;
    get_neighbor(ranges2d, initCluster);
    get_obstacle(obj, initCluster);
}

static void refinement(
    const object<2>& obj,
    std::vector<cluster_cell_vertex>& curClusters, 
    std::vector<cluster_cell_vertex>& nextClusters 
) {
    std::uint32_t prevChildren = 0;
    nextClusters.clear();
    for(auto& curCluster : curClusters) {
        curCluster.invalid_children();
        //一组CNode可能会由多个不连同的子图组成
        std::vector<CNode> cnodes;
        std::vector<Ranges2D> ranges2ds;
        std::vector<std::pair<float, float>> pos0;
        std::vector<cluster_cell_vertex> nextClusters_;

        refine_to_cnode(curCluster, cnodes);
        cnodes_to_ranges2d(cnodes, prevChildren, curCluster, ranges2ds, pos0);
        cnodes.clear();
        nextClusters_.resize(ranges2ds.size());
        for(auto forestIdx = 0U ; forestIdx!= ranges2ds.size() ; ++forestIdx) {
            nextClusters_[forestIdx].pos0[0] = pos0[forestIdx].first;
            nextClusters_[forestIdx].pos0[1] = pos0[forestIdx].second;
            nextClusters_[forestIdx].delta = curCluster.delta / 2.f;
            get_neighbor(ranges2ds[forestIdx], nextClusters_[forestIdx]);
            get_obstacle(obj, nextClusters_[forestIdx]);
        }
        nextClusters.insert(nextClusters.end(), nextClusters_.begin(), nextClusters_.end());
    }
    //finally we concatenate all nodes of cluster and update begin
    std::shared_ptr<std::vector<nodeCV>> nodes = std::make_shared<std::vector<nodeCV>>();
    std::uint32_t currBegin = 0;
    for(auto& nextCluster : nextClusters) {
        nextCluster.begin = currBegin;
        for(auto& node : (*nextCluster.nodes)) {
            if(node.has_down())
                node.down += currBegin;
            if(node.has_up())
                node.up += currBegin;
        }
        assert(nextCluster.nodes->size() == nextCluster.size);
        currBegin += nextCluster.size;
        nodes->insert(nodes->end(), nextCluster.nodes->begin(), nextCluster.nodes->end());
        nextCluster.nodes = nodes;
    }
}

static void get_ranges2d(
    float xmin, float ymin, 
    float xmax, float ymax, 
    float delta, 
    const object<2>& obj, 
    Ranges2D& ranges2d
) {
    using u32 = std::uint32_t;
    using f32 = float;
    assert(xmax > xmin);
    const u32 nx = std::floor((xmax - xmin) / delta);
    assert(ymax > ymin);
    const u32 ny = std::floor((ymax - ymin) / delta);
    u32 row;
    f32 x, y;
    ranges2d.clear();
    ranges2d.resize(ny);
    for(row = 0, y = ymin ; y < ymax && row < ny; ++row, y += delta) {
        bool upDateRangeStart = false;
        u32 rangeStart = 0;
        u32 col;
        for(col = 0, x = xmin ; x < xmax && col < nx ; ++col, x += delta) {
            //whether the node is inside of the object
            bool self = false, left = false, right = false, up = false, down = false;
            bool leftDown = false, leftUp = false, rightDown = false, rightUp = false;

            if(obj.in(x, y)) self = true;
            if(x - delta >= xmin && obj.in(x - delta, y)) left = true;
            if(x + delta < xmax && obj.in(x + delta, y)) right = true;
            if(y + delta < ymax && obj.in(x, y + delta)) up = true;
            if(y - delta >= ymin && obj.in(x, y - delta)) down = true;
            if(x - delta >= xmin && y - delta >= ymin && obj.in(x - delta, y - delta)) leftDown = true;
            if(x - delta >= xmin && y + delta < ymax && obj.in(x - delta, y + delta)) leftUp = true;
            if(x + delta < xmax && y - delta >= ymin && obj.in(x + delta, y - delta)) rightDown = true;
            if(x + delta < xmax && y + delta < ymax && obj.in(x + delta, y + delta)) rightUp = true;

            //This node should be deleted.
            if(self && left && right && up && down && leftDown && leftUp && rightDown && rightUp) {
                if(!upDateRangeStart) {
                    //We found a valid range.
                    ranges2d[row].emplace_back(rangeStart, col);
                    upDateRangeStart = true;
                }
            } else {
                if(upDateRangeStart) {
                    //We start to find a new valid range.
                    rangeStart = col;
                    upDateRangeStart = false;
                }
            }
        }

        if(!upDateRangeStart) {
            //Record the last valid range.
            ranges2d[row].emplace_back(rangeStart, nx);
        }
    }
}

static void get_neighbor(
    const Ranges2D& ranges2d, 
    cluster_cell_vertex& cluster
) {
    using u32 = std::uint32_t;

    //calculate the number of nodes below each row
    std::vector<u32> lowerCnt;
    ranges2d.get_lower(lowerCnt);

    //initialize `nodes`
    cluster.nodes->clear();
    cluster.nodes->resize(lowerCnt.back());
    cluster.begin = 0;
    cluster.size = lowerCnt.back();

    //fill the neighbor information of each node
    for(u32 r = 0 ; r < ranges2d.size() ; ++r) {
        u32 downRangeIdx = 0, 
            currRangeIdx = 0, 
            upRangeIdx = 0;
        u32 downBaseOff = (r == 0) ? -1 : lowerCnt[r-1], 
            currOff = ranges2d[r][0].begin,
            upBaseOff = (r == ranges2d.size() - 1) ? -1 : lowerCnt[r+1];
        //Here `off` is index of the current node in `nodes`.nodes.
        for(u32 off = lowerCnt[r] ; off != lowerCnt[r+1] ; ++off) {

            auto [curBegin, curEnd] = ranges2d[r][currRangeIdx];

            assert(curBegin <= currOff && currOff < curEnd);
            if(r != 0) assert(!ranges2d[r-1].empty());
            if(r != ranges2d.size() - 1) assert(!ranges2d[r+1].empty());
            //Judge whether the node has the left and right neighbor
            (*cluster.nodes)[off].set_left(currOff != curBegin);
            (*cluster.nodes)[off].set_right(currOff != curEnd - 1);
            //Judge whether the node has the up and down neighbor
            if(r != ranges2d.size() - 1 && ranges2d[r+1][upRangeIdx].in(currOff)) {
                (*cluster.nodes)[off].up = upBaseOff + currOff - ranges2d[r+1][upRangeIdx].begin;
            } else {
                (*cluster.nodes)[off].up = static_cast<u32>(-1);
            }
            if(r != 0 && ranges2d[r-1][downRangeIdx].in(currOff)) {
                (*cluster.nodes)[off].down = downBaseOff + currOff - ranges2d[r-1][downRangeIdx].begin;
            } else {
                (*cluster.nodes)[off].down = static_cast<u32>(-1);
            }
            //Update `currRangeIdx` and `currOff`
            ++currOff;
            if(currOff == ranges2d[r][currRangeIdx].end) {
                ++currRangeIdx;
                if(currRangeIdx < ranges2d[r].size()) {
                    currOff = ranges2d[r][currRangeIdx].begin;
                }
            }
            //Update `upRangeIdx` and `upBaseOff`
            while((r != ranges2d.size() - 1) && (upRangeIdx < ranges2d[r+1].size() - 1) && (currOff >= ranges2d[r+1][upRangeIdx].end)) {
                upBaseOff += ranges2d[r+1][upRangeIdx].len();
                ++upRangeIdx;
            }
            //Update `downRangeIdx` and `downBaseOff`
            while((r != 0) && (downRangeIdx < ranges2d[r-1].size() - 1) && (currOff >= ranges2d[r-1][downRangeIdx].end)) {
                downBaseOff += ranges2d[r-1][downRangeIdx].len();
                ++downRangeIdx;
            }
        }
    }
}

static void get_obstacle(
    const object<2>& obj, 
    cluster_cell_vertex& cluster
) {
    bfs(cluster, [&cluster, &obj](long x, long y, std::uint32_t idx)->void{
        float fx  = cluster.pos0[0] + static_cast<float>(x) * cluster.delta;
        float fy = cluster.pos0[1] + static_cast<float>(y) * cluster.delta;
        (*cluster.nodes)[idx].set_obstacle(obj.in(fx, fy));
    });
}

static void bfs(const cluster_cell_vertex& cluster, std::function<void(long, long, std::uint32_t)> option) {
    std::vector<bool> vis (cluster.size, false);
    std::list<std::tuple<long, long, std::uint32_t>> frontier;

    vis[0] = true;
    frontier.emplace_back(0L, 0L, cluster.begin);

    while(!frontier.empty()) {
        auto [x, y, idx] = frontier.front();
        frontier.pop_front();
        assert(cluster.begin <= idx && idx < cluster.begin + cluster.size);

        option(x, y, idx);

        if((*cluster.nodes)[idx].has_left() && !vis[idx - 1]) {
            vis[idx - 1] = true;
            frontier.emplace_back(x-1, y, idx - 1);
        }
        if((*cluster.nodes)[idx].has_right() && !vis[idx + 1]) {
            vis[idx + 1] = true;
            frontier.emplace_back(x + 1, y, idx + 1);
        }
        if((*cluster.nodes)[idx].has_up() && !vis[(*cluster.nodes)[idx].up]) {
            vis[(*cluster.nodes)[idx].up] = true;
            frontier.emplace_back(x, y + 1, (*cluster.nodes)[idx].up);
        }
        if((*cluster.nodes)[idx].has_down() && !vis[(*cluster.nodes)[idx].down]) {
            vis[(*cluster.nodes)[idx].down] = true;
            frontier.emplace_back(x, y-1, (*cluster.nodes)[idx].down);
        }
    }
}

static void bfs(
    const std::vector<CNode>& cnodes, 
    std::function<void(long, long, std::uint32_t, std::uint32_t)> option
) {
    using u32 = std::uint32_t;
    std::vector<bool> vis(cnodes.size(), false);
    u32 forestIdx = 0;
    u32 unVis = 0;
    while(unVis < cnodes.size()) {
        if(vis[unVis]) {
            ++unVis;
            continue;
        }

        std::list<std::tuple<long, long, u32>> frontier;

        vis[unVis] = true;
        frontier.emplace_back(0L, 0L, unVis);

        while(!frontier.empty()) {
            auto [x, y, idx] = frontier.front();
            frontier.pop_front();

            option(x, y, idx, forestIdx);

            if(cnodes[idx].has_left() && !vis[cnodes[idx].left]) {
                vis[cnodes[idx].left] = true;
                frontier.emplace_back(x-1, y, cnodes[idx].left);
            }

            if(cnodes[idx].has_right() && !vis[cnodes[idx].right]) {
                vis[cnodes[idx].right] = true;
                frontier.emplace_back(x + 1, y, cnodes[idx].right);
            }

            if(cnodes[idx].has_up() && !vis[cnodes[idx].up]) {
                vis[cnodes[idx].up] = true;
                frontier.emplace_back(x, y + 1, cnodes[idx].up);
            }

            if(cnodes[idx].has_down() && !vis[cnodes[idx].down]) {
                vis[cnodes[idx].down] = true;
                frontier.emplace_back(x, y-1, cnodes[idx].down);
            }
        }

        ++forestIdx;
    }
}

static void to_dense(
    const std::vector<CNode>& nodes, 
    std::uint32_t nx, std::uint32_t ny, 
    std::uint32_t x0, std::uint32_t y0, 
    std::vector<std::vector<bool>>& grids
) {
    grids.clear();

    bfs(nodes, [nx, ny, x0, y0, &grids](long x, long y, std::uint32_t idx, std::uint32_t forestIdx)->void{
        if(forestIdx >= grids.size())
            grids.resize(forestIdx + 1, std::vector<bool>(nx * ny, false));

        grids[forestIdx][(y + y0) * nx + (x + x0)] = true;
    });
}

static void refine_to_cnode(
    cluster_cell_vertex& curCluster, 
    std::vector<CNode>& cnodes
) {
    bfs(curCluster, [&curCluster, &cnodes](long x, long y, std::uint32_t idx)->void{
        using u32 = std::uint32_t;
        using f32 = float;
        //iterate all bottom left corner

        bool exist01 = (*curCluster.nodes)[idx].has_right();
        std::uint32_t i01 = exist01 ? idx + 1 : static_cast<u32>(-1);
        bool exist10 = (*curCluster.nodes)[idx].has_up();
        std::uint32_t i10 = exist10 ? (*curCluster.nodes)[idx].get_up() : static_cast<u32>(-1);
        bool exist11 = exist01 && (*curCluster.nodes)[i01].has_up();
        std::uint32_t i11 = exist11 ? (*curCluster.nodes)[i01].get_up() : static_cast<u32>(-1);

        if(exist01 && exist10 && exist11 && 
            ((*curCluster.nodes)[idx].is_obstacle() || (*curCluster.nodes)[i01].is_obstacle() || (*curCluster.nodes)[i10].is_obstacle() || (*curCluster.nodes)[i11].is_obstacle()) && 
            !((*curCluster.nodes)[idx].is_obstacle() && (*curCluster.nodes)[i01].is_obstacle() && (*curCluster.nodes)[i10].is_obstacle() && (*curCluster.nodes)[i11].is_obstacle())
        ) {
            const float x00 = curCluster.pos0[0] + static_cast<f32>(x) * curCluster.delta;
            const float y00 = curCluster.pos0[1] + static_cast<f32>(y) * curCluster.delta;
            const float delta = curCluster.delta;
            bool hasChild00 = (*curCluster.nodes)[idx].has_child();
            bool hasChild01 = (*curCluster.nodes)[i01].has_child();
            bool hasChild10 = (*curCluster.nodes)[i10].has_child();
            bool hasChild11 = (*curCluster.nodes)[i11].has_child();

            //try to mount cnode to four corner parent nodes
            if(!hasChild00) {
                (*curCluster.nodes)[idx].child = cnodes.size();
                cnodes.emplace_back(idx, x00, y00);
            }

            if(!hasChild01) {
                (*curCluster.nodes)[i01].child = cnodes.size();
                cnodes.emplace_back(i01, x00 + delta, y00);
            }

            if(!hasChild10) {
                (*curCluster.nodes)[i10].child = cnodes.size();
                cnodes.emplace_back(i10, x00, y00 + delta);
            }

            if(!hasChild11) {
                (*curCluster.nodes)[i11].child = cnodes.size();
                cnodes.emplace_back(i11, x00 + delta, y00 + delta);
            }

            std::uint32_t ci00 = (*curCluster.nodes)[idx].child, 
                        ci01 = (*curCluster.nodes)[i01].child, 
                        ci10 = (*curCluster.nodes)[i10].child, 
                        ci11 = (*curCluster.nodes)[i11].child;
            std::uint32_t ciLeft, ciRight, ciUp, ciDown;

            //try to mount cnode to four edge nodes

            if(
                !(hasChild00 && hasChild10) || 
                (hasChild00 && !cnodes[ci00].has_up())
            ) {
                ciLeft = cnodes.size();
                cnodes.emplace_back(x00, y00 + delta / 2.f);
                cnodes[ciLeft].make_edge_center_y(cnodes[ci10], cnodes[ci00], ciLeft, ci10, ci00);
            } else {
                ciLeft = cnodes[ci00].up;
            }

            if(
                !(hasChild01 && hasChild11) || 
                (hasChild01 && !cnodes[ci01].has_up())
            ) {
                ciRight = cnodes.size();
                cnodes.emplace_back(x00 + delta, y00 + delta / 2.f);
                cnodes[ciRight].make_edge_center_y(cnodes[ci11], cnodes[ci01], ciRight, ci11, ci01);
            } else {
                ciRight = cnodes[ci01].up;
            }

            if(
                !(hasChild00 && hasChild01) ||
                (hasChild00 && !cnodes[ci00].has_right())
            ) {
                ciDown = cnodes.size();
                cnodes.emplace_back(x00 + delta / 2.f, y00);
                cnodes[ciDown].make_edge_center_x(cnodes[ci00], cnodes[ci01], ciDown, ci00, ci01);
            } else {
                ciDown = cnodes[ci00].right;
            }

            if(
                !(hasChild10 && hasChild11) || 
                (hasChild10 && !cnodes[ci10].has_right())
            ) {
                ciUp = cnodes.size();
                cnodes.emplace_back(x00 + delta / 2.f, y00 + delta);
                cnodes[ciUp].make_edge_center_x(cnodes[ci10], cnodes[ci11], ciUp, ci10, ci11);
            } else {
                ciUp = cnodes[ci10].right;
            }

            //mount the face center
            std::uint32_t ciCent = cnodes.size();
            cnodes.emplace_back(x00 + delta / 2.f, y00 + delta / 2.f);
            cnodes[ciCent].make_face_center(
                cnodes[ciLeft], cnodes[ciRight], cnodes[ciUp], cnodes[ciDown], 
                ciCent, ciLeft, ciRight, ciUp, ciDown
            );
        }
        
    });
}

static void cnodes_to_ranges2d(
    const std::vector<CNode>& cnodes,
    std::uint32_t& prevChildren, 
    cluster_cell_vertex& parCluster,
    std::vector<Ranges2D>& ranges2ds, 
    std::vector<std::pair<float, float>>& pos0
) {
    using i64 = std::int64_t;
    using u32 = std::uint32_t;
    std::vector<i64> ymin;
    std::vector<i64> ymax;

    //get ymin and ymax of each forest
    bfs(cnodes, [&ymin, &ymax](long x, long y, u32 idx, u32 forestIdx)->void{
        if(forestIdx >= ymin.size()) ymin.resize(forestIdx + 1, std::numeric_limits<i64>::max());
        if(forestIdx >= ymax.size()) ymax.resize(forestIdx + 1, std::numeric_limits<i64>::min());

        ymin[forestIdx] = std::min(ymin[forestIdx], y);
        ymax[forestIdx] = std::max(ymax[forestIdx], y);
    });

    std::vector<std::vector<std::vector<std::pair<long, std::uint32_t>>>> pairXParent(ymin.size());
    pos0.clear();
    pos0.resize(ymin.size(), std::make_pair(1e20f, 1e20f));

    for(auto forestIdx = 0U ; forestIdx != ymin.size() ; ++forestIdx) {
        pairXParent[forestIdx].resize(ymax[forestIdx] - ymin[forestIdx] + 1);
    }

    bfs(cnodes, [&pairXParent, &cnodes, &ymin, &pos0](long x, long y, u32 idx, u32 forestIdx)->void{
        pairXParent[forestIdx][y - ymin[forestIdx]].emplace_back(x, cnodes[idx].parent);
        if(std::abs(cnodes[idx].y - pos0[forestIdx].second) < 1e-7f) {
            pos0[forestIdx].first = std::min(pos0[forestIdx].first, cnodes[idx].x);
        } else if(cnodes[idx].y < pos0[forestIdx].second) {
            pos0[forestIdx].second = cnodes[idx].y;
            pos0[forestIdx].first = cnodes[idx].x;
        }
    });

    for(auto& forest: pairXParent) {
        for(auto& row : forest) {
            std::sort(row.begin(), row.end(), [](const auto& lhs, const auto& rhs)->bool{
                return lhs.first < rhs.first;
            });
        }
    }

    ranges2ds.resize(pairXParent.size());
    auto prevCnt = 0U;
    std::vector<u32> lowerCnt(pairXParent.size() + 1);
    for(auto forestIdx = 0U ; forestIdx != pairXParent.size() ; ++forestIdx) {
        ranges2ds[forestIdx].resize(pairXParent[forestIdx].size());
        auto xmin = std::numeric_limits<i64>::max();
        for(auto& row : pairXParent[forestIdx]) xmin = std::min(xmin, row.front().first);

        for(auto rowIdx = 0U ; rowIdx != pairXParent[forestIdx].size() ; ++rowIdx) {
            for(auto [x, parentIdx] : pairXParent[forestIdx][rowIdx]) {
                if(ranges2ds[forestIdx][rowIdx].empty())
                    ranges2ds[forestIdx][rowIdx].emplace_back(x-xmin, x-xmin);
                if(x - xmin == ranges2ds[forestIdx][rowIdx].back().end) {
                    ++ranges2ds[forestIdx][rowIdx].back().end;
                } else {
                    ranges2ds[forestIdx][rowIdx].emplace_back(x-xmin, x-xmin + 1);
                }
                if(parentIdx < parCluster.nodes->size())
                    (*parCluster.nodes)[parentIdx].child = prevChildren;
                ++prevChildren;
            }
        }

        pairXParent[forestIdx].clear();
    }
}
cell_vertex_layer_2d::cell_vertex_layer_2d(float xmin, float ymin, float xmax, float ymax, float delta, const object<2>& obj) : 
    m_nodes(), 
    m_delta(delta), 
    m_layer(0), 
    m_clusters(0) {
    cluster_cell_vertex initCluster;
    init_first_cluster(xmin, ymin, xmax, ymax, delta, obj, initCluster);
    m_nodes = initCluster.nodes;
    m_clusters.emplace_back(initCluster.begin, initCluster.size, initCluster.pos0[0], initCluster.pos0[1]);
}

cell_vertex_layer_2d cell_vertex_layer_2d::refine(const object<2>& obj) {
    std::vector<cluster_cell_vertex> currClusters;
    std::vector<cluster_cell_vertex> nextClusters;
    for(const auto& cluster : m_clusters) {
        currClusters.emplace_back(m_nodes, cluster.begin, cluster.size, std::array<float, 2>{cluster.x, cluster.y}, m_delta);
    }
    refinement(obj, currClusters, nextClusters);
    currClusters.clear();
    cell_vertex_layer_2d nextLayer;
    nextLayer.m_nodes = nextClusters[0].nodes;
    nextLayer.m_delta = m_delta / 2.f;
    nextLayer.m_layer = m_layer + 1;
    for(const auto& cluster : nextClusters) {
        nextLayer.m_clusters.emplace_back(cluster.begin, cluster.size, cluster.pos0[0], cluster.pos0[1]);
    }
    return nextLayer;
}

void cell_vertex_layer_2d::cluster::traverse(std::function<void(long, long, long)> f) const {
    std::vector<bool> vis(m_nodes.size(), false);
    //<local x, local y, global off>
    std::list<std::tuple<long, long, long>> frontier;
    vis[0] = true;
    frontier.emplace_back(0, 0, m_begin);
    while(!frontier.empty()) {
        auto [x, y, idx] = frontier.front();
        frontier.pop_front();

        f(x, y, idx);

        std::uint32_t localIdx = idx - m_begin;

        if(m_nodes[localIdx].has_left() && !vis[localIdx - 1]) {
            vis[localIdx - 1] = true;
            frontier.emplace_back(x-1, y, idx - 1);
        }

        if(m_nodes[localIdx].has_right() && !vis[localIdx + 1]) {
            vis[localIdx + 1] = true;
            frontier.emplace_back(x+1, y, idx + 1);
        }

        if(m_nodes[localIdx].has_up() && !vis[m_nodes[localIdx].get_up() - m_begin]) {
            vis[m_nodes[localIdx].get_up() - m_begin] = true;
            frontier.emplace_back(x, y+1, m_nodes[localIdx].get_up());
        }

        if(m_nodes[localIdx].has_down() && !vis[m_nodes[localIdx].get_down() - m_begin]) {
            vis[m_nodes[localIdx].get_down() - m_begin] = true;
            frontier.emplace_back(x, y-1, m_nodes[localIdx].get_down());
        }
    }
}

void cell_vertex_layer_2d::cluster::to_point2d(std::ostream& os) const {
    float x0 = m_x, y0 = m_y, delta = m_delta;
    traverse([&os, x0, y0, delta](long x, long y, [[maybe_unused]] long idx )->void{
        os << std::format("{}, {}\n", x0 + x * delta, y0 + y * delta);
    });
}

void cell_vertex_layer_2d::cluster::to_dense(std::vector<bool>& denseGrid, u32 nx, u32 ny, u32 x0, u32 y0) const {
    denseGrid.clear();
    denseGrid.resize(nx * ny, false);
    traverse([&denseGrid, nx, ny, x0, y0](long x, long y, [[maybe_unused]] long idx)->void{
        x += x0;
        y += y0;
        assert(0 <= x && x < nx);
        assert(0 <= y && y < ny);
        denseGrid[y * nx + x] = true;
    });
}

void cell_vertex_layer_2d::cluster::obstacle_to_dense(std::vector<bool>& denseGrid, u32 nx, u32 ny, u32 x0, u32 y0) const {
    denseGrid.clear();
    denseGrid.resize(nx * ny, false);
    traverse([this, &denseGrid, nx, ny, x0, y0](long x, long y, long idx)->void{
        x += x0;
        y += y0;
        assert(0 <= x && x < nx);
        assert(0 <= y && y < ny);
        if(m_nodes[idx - m_begin].is_obstacle())
            denseGrid[y * nx + x] = true;
    });
}


