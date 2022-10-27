//----------------------------------------------------------------------------//
//                                  includes                                  //
//----------------------------------------------------------------------------//

// trigger optimisation from source file
// # pragma GCC optimize("Ofast", "inline", "omit-frame-pointer", "unroll-loops")
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>


constexpr int W = 30;
constexpr int H = 16;

typedef char Grid[W * H];
typedef uint16_t Pos;

/*
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
*/

inline Pos get_pos(int x, int y) {
    if (x < 0 ||x >= W || y < 0 || y >= H) {
        throw std::runtime_error(std::string("Error! Invalid coordinates [x") + std::to_string(x) + ", y" + std::to_string(y) + "]");
    }
    return y * W + x;
}

inline Pos get_pos(std::pair<int, int> const& c) {
    return get_pos(c.first, c.second);
}

inline std::pair<int, int> get_coords(Pos const& p) {
    if (p >= W * H) {
        throw std::runtime_error(std::string("Error! Invalid Pos [") + std::to_string(p) + "]");
    }
    return {p % W, p / W};
}

std::ostream& operator<<(std::ostream& os, std::pair<int, int> const& coords) {
    return os << "(x=" << coords.first << ", y=" << coords.second << ")";
}

template<class T>
std::ostream& operator<<(std::ostream& os, std::set<T> const& s) {
    os << "{";
    for (auto it = std::begin(s); it != std::end(s); ++it) {
        if (it != std::begin(s)) {
            os << ", ";
        }
        os << *it;
    }
    return os << "}";
}

template<class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& v) {
    os << "[";
    for (auto it = std::begin(v); it != std::end(v); ++it) {
        if (it != std::begin(v)) {
            os << ", ";
        }
        os << *it;
    }
    return os << "]";
}

template<class K, class V>
std::ostream& operator<<(std::ostream& os, std::map<K, V> const& m) {
    os << "<";
    for (auto it = std::begin(m); it != std::end(m); ++it) {
        if (it != std::begin(m)) {
            os << ", ";
        }
        K const& k = it->first;
        V const& v = it->second; 
        os << k << " => " << v;
    }
    return os << ">";
}


struct Group {
    std::set<Pos> known_cells;
    std::set<Pos> unknown_cells;
    // Possibilities store bombs
    std::vector<std::set<Pos>> possibilities;

    Group(std::set<Pos> const& known_cells, std::set<Pos> const& unknown_cells): known_cells {known_cells}, unknown_cells {unknown_cells}, possibilities {{}} {}

    void merge(Group const* group, std::set<Pos> const& common_cells) {
        // std::cerr << "Merging groups with common key " << common_cells << ":\n";
        std::map<std::set<Pos>, std::vector<std::set<Pos>>> my_poss_by_common;
        for (std::set<Pos> const& poss: possibilities) {
            std::set<Pos> poss_key;
            for (Pos const& cell: common_cells) {
                if (poss.find(cell) != end(poss)) {
                    poss_key.insert(cell);
                }
            }
            my_poss_by_common[poss_key].push_back(poss);
        }
        // std::cerr << "Group " << unknown_cells << " by common key: " << my_poss_by_common << "\n";
        std::vector<std::set<Pos>> new_possibilities;
        // std::cerr << "Group " << group->unknown_cells << " by common key:\n";
        for (std::set<Pos> const& g_poss: group->possibilities) {
            std::set<Pos> g_poss_key;
            for (Pos const& cell: common_cells) {
                if (g_poss.find(cell) != end(g_poss)) {
                    g_poss_key.insert(cell);
                }
            }
            auto it = my_poss_by_common.find(g_poss_key);
            // std::cerr << "Poss " << g_poss << " (common: " << g_poss_key << ") Matches ";
            if (it == std::end(my_poss_by_common)) {
            //     std::cerr << "NONE\n";
                continue;
            }
            // std::cerr << it->second << "\n";
            for (std::set<Pos> const& my_poss: it->second) {
                new_possibilities.emplace_back(g_poss);
                new_possibilities.back().insert(std::begin(my_poss), std::end(my_poss));
            }
        }
        known_cells.insert(std::begin(group->known_cells), std::end(group->known_cells));
        unknown_cells.insert(std::begin(group->unknown_cells), std::end(group->unknown_cells));
        possibilities = new_possibilities;
    }

    // std::set<Pos> safe_cells(void) const {
    //     std::set<Pos> res {unknown_cells};
    //     for (auto const& poss: possibilities) {
    //         res.erase(std::begin(poss), std::end(poss));
    //     }
    //     return res;
    // }

    Pos safest_cell(double & safety_level) const {
        std::map<Pos, size_t> cells_safety;

        for (Pos p: unknown_cells) {
            cells_safety[p] = possibilities.size();
        }
        for (auto const& poss: possibilities) {
            for (Pos p: poss) {
                --cells_safety[p];
            }
        }
        Pos res = *std::max_element(std::begin(unknown_cells), std::end(unknown_cells), [cells_safety](Pos p1, Pos p2) {return cells_safety.at(p1) < cells_safety.at(p2);});
        safety_level = static_cast<double>(cells_safety[res]) / possibilities.size();
        return res;
    }
};

std::ostream& operator<<(std::ostream& os, Grid const& grid) {
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            if (x > 0) {
                os << " ";
            }
            os << grid[get_pos(x, y)];
        }
        os << "\n";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, Group const& group) {
    int x_min = W - 1;
    int x_max = 0;
    int y_min = H - 1;
    int y_max = 0;

    os << "Group " << group.unknown_cells << "\n";
    for (Pos const& p: group.unknown_cells) {
        std::pair<int, int> coords = get_coords(p);
        x_min = std::min(x_min, coords.first);
        x_max = std::max(x_max, coords.first);
        y_min = std::min(y_min, coords.second);
        y_max = std::max(y_max, coords.second);
    }
    for (std::set<Pos> const& possibility: group.possibilities) {
        os << "Poss " << possibility << "\n";
        for (int y = y_min; y <= y_max; ++y) {
            os << "|";
            for (int x = x_min; x <= x_max; ++x) {
                if (x > 0) {
                    os << "";
                }
                Pos p = get_pos(x, y);
                if (possibility.find(p) != std::end(possibility)) {
                    os << "x";
                } else if (group.unknown_cells.find(p) != std::end(group.unknown_cells)) {
                    os << "o";
                } else {
                    os << ".";
                }
            }
            os << "|\n";
        }
    }
    return os;
}

void read_grid(Grid & grid) {
    for (auto & c: grid) {
        std::cin >> c;
        std::cin.ignore();
    }
}

template<typename T>
std::set<T> intersection(std::set<T> const& a, std::set<T> const& b) {
    std::set<T> res;
    for (auto const& elt: a) {
        if (b.find(elt) != std::end(b)) {
            res.insert(elt);
        }
    }
    return res;
}

int cell_val(Grid const& grid, Pos const& p) {
    char cell_char = grid[p];
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

/**
 * @brief Note for optimisation purpose, p is its own neighbour.
 * Should have no practical impact as this function is always called on numbered cells.
 */
template<class OutputIterator>
OutputIterator unknown_neighbours(Grid const& grid, Pos const& p, OutputIterator const res) {
    OutputIterator it = res;
    std::pair<int, int> c = get_coords(p);
    
    for (int n_y = std::max(0, c.second - 1); n_y < std::min(H, c.second + 2); ++n_y) {
        for (int n_x = std::max(0, c.first - 1); n_x < std::min(W, c.first + 2); ++n_x) {
            if (grid[get_pos(n_x, n_y)] != '?') {
                continue;
            }
            *it = get_pos(n_x, n_y);
            ++it;
        }
    }
    return (it);
}

Group * cell_group(Grid const& grid, Pos const& p) {
    static Pos neighbours[8];
    size_t bombs_around = cell_val(grid, p);
    auto neighbours_end = unknown_neighbours(grid, p, std::begin(neighbours));
    size_t neighbours_count = std::distance(std::begin(neighbours), neighbours_end);
    Group * group = new Group({p}, {std::begin(neighbours), neighbours_end});

    for (size_t ni = 0; ni < neighbours_count; ++ni) {
        std::vector<std::set<Pos>> next_possibilities;
        for (std::set<Pos> const& poss: group->possibilities) {
            if (poss.size() < bombs_around) {
                // Possibility of a bomb at neighbours[ni]
                next_possibilities.emplace_back(poss);
                next_possibilities.back().insert(neighbours[ni]);
            } 
            if (bombs_around + ni < neighbours_count + poss.size()) {   // <=> (bombs_around - poss.size() < neighbours_count - ni) avoiding substraction with unsigned
                // Possibility of no bomb at neighbours[ni]
                next_possibilities.emplace_back(poss);
            }
        }
        group->possibilities = next_possibilities;
    }
    return (group);
}

void compute_cell(Grid const& grid, Pos const& p, std::vector<Group*> & groups) {
    std::vector<int> to_del;
    Group * c_group = cell_group(grid, p);

    // std::cerr << "cell_group(" << p << ") = " << *c_group << "\n";
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        std::set<Pos> common_cells = intersection(c_group->unknown_cells, groups[gi]->unknown_cells);
        if (!common_cells.empty()) {
            c_group->merge(groups[gi], common_cells);
            // std::cerr << "Merged => " << *(c_group) << "\n";
            delete groups[gi];
            groups[gi] = nullptr;
        }
    }
    auto it = std::remove(std::begin(groups), std::end(groups), nullptr);
    groups.erase(it, std::end(groups));
    groups.push_back(c_group);
}

int main(void) {
    std::set<Pos> computed_cells;
    std::vector<Group*> groups;
    Grid grid;
    read_grid(grid);
    std::cout << W / 2 << " " << H / 2 << std::endl;    // First move in the middle of the grid

    while (true) {
        read_grid(grid);
        // std::cerr << grid;
        // std::cerr << "Already computed cells: " << computed_cells << "\n";
        for (Pos p = 0; p < W * H; ++p) {
            if ('0' <= grid[p] && grid[p] <= '9' && computed_cells.find(p) == std::end(computed_cells)) {
                std::cerr << "Computing cell " << p << " " << get_coords(p) << "\n";
                computed_cells.insert(p);
                compute_cell(grid, p, groups);
            }
        }
        Pos guess;
        double guess_safety = 0.0;
        std::vector<Group *>::iterator guess_group_it = std::end(groups);
        // std::cerr << groups.size() << " groups.\n";
        for (auto group_it = std::begin(groups); group_it != std::end(groups); ++group_it) {
            Group const* g = *group_it;
            std::cerr << "Group #" << std::distance(std::begin(groups), group_it) + 1 << " (" << g->possibilities.size() << " possibilities), "
                    << g->possibilities.front().size() << " bombs, cells: " << g->unknown_cells << "\n";
            double cell_safety;
            Pos p = g->safest_cell(cell_safety);
            if (cell_safety > guess_safety) {
                guess = p;
                guess_safety = cell_safety;
                guess_group_it = group_it;
                if (cell_safety == 1.0) {   // This cell is perfectly safe
                    break;
                }
            }
        }
        std::pair<int, int> coords = get_coords(guess);
        std::cerr << "Guess safety = " << guess_safety * 100.0 << "%\n"; 
        std::cout << coords.first << " " << coords.second << std::endl;
        Group * g = *guess_group_it;

        std::cerr << "Forgetting cells " << g->known_cells << "\n";
        for (Pos const& p: g->known_cells) {
            computed_cells.erase(p);
        }
        groups.erase(guess_group_it);
        // delete g;
        // throw std::runtime_error("stop after first turn!");
    }
    return (0);
}
