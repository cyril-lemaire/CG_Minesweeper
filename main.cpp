//----------------------------------------------------------------------------//
//                                  includes                                  //
//----------------------------------------------------------------------------//

// trigger optimisation from source file
# pragma GCC optimize("Ofast", "inline", "omit-frame-pointer", "unroll-loops")
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>

using namespace std;

constexpr int W = 30;
constexpr int H = 16;

typedef char Grid[H][W];

struct Pos {
    int x, y;

    int cmp(Pos const& p) const {
        return (x - p.x || y - p.y);
    }
    bool operator==(Pos const& p) const {
        return (cmp(p) == 0);
    }
    bool operator!=(Pos const& p) const {
        return (cmp(p) != 0);
    }
    bool operator<(Pos const& p) const {
        return (cmp(p) < 0);
    }
    bool operator<=(Pos const& p) const {
        return (cmp(p) <= 0);
    }
    bool operator>(Pos const& p) const {
        return (cmp(p) > 0);
    }
    bool operator>=(Pos const& p) const {
        return (cmp(p) >= 0);
    }
};

struct Group {
    set<Pos> cells;
    // Possibilities store bombs
    vector<set<Pos>> possibilities;

    Group(set<Pos> const& cells): cells {cells}, possibilities {{}} {}

    void merge(Group const* group, set<Pos> const& common_cells) {
        map<set<Pos>, vector<set<Pos>>> my_poss_by_common;
        for (set<Pos> const& poss: possibilities) {
            set<Pos> poss_key;
            for (Pos const& cell: common_cells) {
                if (poss.find(cell) != end(poss)) {
                    poss_key.insert(cell);
                }
            }
            my_poss_by_common[poss_key].push_back(poss);
        }
        vector<set<Pos>> new_possibilities;
        for (set<Pos> const& g_poss: group->possibilities) {
            set<Pos> g_poss_key;
            for (Pos const& cell: common_cells) {
                if (g_poss.find(cell) != end(g_poss)) {
                    g_poss_key.insert(cell);
                    auto it = my_poss_by_common.find(g_poss_key);
                    if (it == end(my_poss_by_common)) {
                        continue;
                    }
                    for (set<Pos> const& g_poss: it->second) {
                        new_possibilities.emplace_back(g_poss);
                    }
                }
            }

        }
    }
};

ostream& operator<<(ostream& os, Grid const& grid) {
    for (int l = 0; l < H; ++l) {
        os << "|";
        for (int c = 0; c < W; ++c) {
            if (c > 0)
                os << "";
            os << grid[l][c];
        }
        os << "|\n";
    }
    return os;
}

void read_grid(Grid & grid) {
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            cin >> grid[y][x];
            cin.ignore();
        }
    }
}

template<typename T>
set<T> intersection(set<T> const& a, set<T> const& b) {
    set<T> res;
    for (auto const& elt: a) {
        if (b.find(elt) != end(b)) {
            res.insert(elt);
        }
    }
    return res;
}

int cell_val(Grid const& grid, Pos const& p) {
    char cell_char = grid[p.y][p.x];
    int val;

    switch (cell_char) {
        case '?':
            val = -1;
            break;
        case '.':
            val = 0;
            break;
        default:
            val = cell_char - '0';
            break;
    }
    return val;
}

template<class OutputIterator>
int unknown_around(Grid const& grid, Pos const& p, OutputIterator res) {
    int count = 0;
    Pos n;

    for (n.y = max(0, p.y - 1); n.y < min(H, p.y + 2); ++n.y) {
        for (n.x = max(0, p.x - 1); n.y < min(W, p.x + 2); ++n.x) {
            if (grid[n.y][n.x] == '?' && (n.y != p.y || n.x != p.x)) {
                *res = n;
                ++res;
                ++count;
            }
        }
    }
    return (count);
}

Group * cell_group(Grid const& grid, Pos const& p) {
    static Pos neighbours[9];
    int bombs_around = cell_val(grid, p);
    int neighbours_count = unknown_around(grid, p, begin(neighbours));

    Group * group = new Group({begin(neighbours), begin(neighbours) + neighbours_count});

    for (int ni = 0; ni < neighbours_count; ++ni) {
        vector<set<Pos>> next_possibilities;
        for (set<Pos> const& poss: group->possibilities) {
            if (static_cast<int>(poss.size()) < bombs_around) {
                // Possibility of a bomb at neighbours[ni]
                next_possibilities.emplace_back(poss);
                next_possibilities.back().insert(neighbours[ni]);
            }
            if (bombs_around - static_cast<int>(poss.size()) < neighbours_count - ni) {
                // Possibility of no bomb at neighbours[ni]
                next_possibilities.emplace_back(poss);
            }
        }
        group->possibilities = next_possibilities;
    }
    return (group);
}

void compute_cell(Grid const& grid, Pos const& p, vector<Group*> & groups) {
    vector<int> to_del;
    Group * c_group = cell_group(grid, p);

    for (size_t gi = 0; gi < groups.size(); ++gi) {
        set<Pos> common_cells = intersection(c_group->cells, groups[gi]->cells);
        if (!common_cells.empty()) {
            c_group->merge(groups[gi], common_cells);
            delete groups[gi];
            groups[gi] = nullptr;
        }
    }
    auto it = remove(begin(groups), end(groups), nullptr);
    groups.erase(it, end(groups));
    groups.push_back(c_group);
}

vector<Pos>

int main(void) {
    Grid grid;
    while (true) {
        read_grid(grid);
        cerr << grid;
    }
    return (0);
}
