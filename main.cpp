//----------------------------------------------------------------------------//
//                                  includes                                  //
//----------------------------------------------------------------------------//

// trigger optimisation from source file
// # pragma GCC optimize("Ofast", "inline", "omit-frame-pointer", "unroll-loops")
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <array>
#include <map>
#include <queue>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>
#include <functional>


constexpr int W = 30;
constexpr int H = 16;
constexpr int CELLS_COUNT = W * H;

constexpr std::chrono::duration FIRST_TURN_DURATION = std::chrono::seconds(1);
constexpr std::chrono::duration TURN_DURATION = std::chrono::milliseconds(50);
constexpr std::chrono::duration TIME_EPSILON = std::chrono::microseconds(1000);

std::condition_variable game_player_cv;

typedef std::array<char, CELLS_COUNT> Grid;
typedef uint16_t Pos;

inline Pos get_pos(int x, int y) {
    if (x < 0 || x >= W || y < 0 || y >= H) {
        throw std::runtime_error(std::string("Error! Invalid coordinates [x") + std::to_string(x) + ", y" + std::to_string(y) + "]");
    }
    return y * W + x;
}

inline Pos get_pos(std::pair<int, int> const& c) {
    return get_pos(c.first, c.second);
}

inline std::pair<int, int> get_coords(Pos const& p) {
    if (p >= CELLS_COUNT) {
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

int cell_val(char cell_char) {
    switch (cell_char) {
        case '?':
            return -1;
        case '.':
            return 0;
        case '1':
            [[fallthrough]];
        case '2':
            [[fallthrough]];
        case '3':
            [[fallthrough]];
        case '4':
            [[fallthrough]];
        case '5':
            [[fallthrough]];
        case '6':
            [[fallthrough]];
        case '7':
            [[fallthrough]];
        case '8':
            return cell_char - '0';
        default:
            throw std::runtime_error("Error! Invalid cell value [" + std::to_string(cell_char) + "]");
    }
}

/**
 * @brief Counts neighbour of p in grid.
 * Note: for optimisation purpose, p is considered its own neighbour.
 * Should have no practical impact as this function is always called on numbered - thus revealed - cells.
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

Group * raw_cell_group(Grid const& grid, Pos const& p) {
    static Pos neighbours[8];
    size_t bombs_around = cell_val(grid[p]);
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

Group * get_cell_group(Grid const& grid, Pos const& p, std::vector<Group*> & groups) {
    std::vector<int> to_del;
    Group * c_group = raw_cell_group(grid, p);

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
    return c_group;
}

void game_interface(Grid & grid, std::mutex & grid_mtx, std::queue<Pos> & to_compute, std::mutex to_compute_mtx, std::queue<Pos> & to_play, Pos & to_guess, std::mutex & to_play_mtx) {
    std::chrono::high_resolution_clock::time_point turn_start;
    std::chrono::high_resolution_clock::time_point turn_end;
    std::set<Pos> computed_cells;
    Pos played_pos;

    to_guess = W * H;   // This indicates the game_player that the guess was played.
    read_grid(grid);
    std::cout << W / 2 << " " << H / 2 << std::endl;    // First move in the middle of the grid

    while (true) {
        {
            std::lock_guard<std::mutex> lck(grid_mtx);  // Only this thread can write, block against concurrent read access 
            read_grid(grid);
        }
        turn_start = std::chrono::high_resolution_clock::now();
        turn_end = turn_start + TURN_DURATION;
        for (Pos p = 0; p < CELLS_COUNT; ++p) {
            if ('1' <= grid[p] && grid[p] <= '1' && computed_cells.find(p) == std::end(computed_cells)) {
                {
                    std::lock_guard<std::mutex> lck(to_compute_mtx);
                    to_compute.push(p);
                    computed_cells.insert(p);
                    game_player_cv.notify_one();
                }
            }
        }
        std::this_thread::sleep_until(turn_end - TIME_EPSILON);
        {
            std::lock_guard<std::mutex> lck(to_play_mtx);
            if (!to_play.empty()) {
                played_pos = to_play.front();
                to_play.pop();
            } else {
                played_pos = to_guess;
                to_guess = W * H;   // This indicates game_player thread that the guess was played.
            }
        }
        std::pair<int, int> coords = get_coords(played_pos);
        std::cout << coords.first << " " << coords.second << std::endl;

        // throw std::runtime_error("stop after first turn!");
    }
}

void game_player(Grid & grid, std::mutex & grid_mtx, std::queue<Pos> & to_compute, std::mutex to_compute_mtx, std::queue<Pos> & to_play, Pos & to_guess, std::mutex & to_play_mtx) {
    std::vector<Group*> groups;
    std::array<double, CELLS_COUNT> cells_safety;
    double guess_safety;
    Pos cell;
    Group * cell_group;
    std::unique_lock<std::mutex> to_compute_lck {to_compute_mtx};

    std::fill(std::begin(cells_safety), std::end(cells_safety), 0.0);
    guess_safety = 0.0;

    while (true) {
        game_player_cv.wait(to_compute_lck, [&to_compute]{return !to_compute.empty();});
        if (to_guess < 0 || CELLS_COUNT <= to_guess) {
            guess_safety = 0.0;
            for (Group const* g: groups) {
                for (Pos p: g->unknown_cells) {
                    if (guess_safety < cells_safety[p] && cells_safety[p] < 1.0) {
                        std::lock_guard<std::mutex> lck(to_play_mtx);
                        to_guess = p;
                        guess_safety = cells_safety[p];
                    }
                }
            }
        }
        {
            std::lock_guard<std::mutex> lck(to_compute_mtx);
            cell = to_compute.front();
            to_compute.pop();
        }
        {
            std::lock_guard<std::mutex> lck(grid_mtx);
            cell_group = get_cell_group(grid, cell, groups);
        }
        for (Pos p: cell_group->unknown_cells) {
            cells_safety[p] = 1.0;
        }
        for (auto const& poss: cell_group->possibilities) {
            double const step = 1.0 / poss.size();
            for (Pos p: poss) {
                cells_safety[p] -= step;
            }
        }
        for (Pos const& p: cell_group->unknown_cells) {
            if (cells_safety[p] == 1.0) {    // Safe cell
                {
                    std::lock_guard<std::mutex> lck(to_play_mtx);
                    to_play.push(p);
                }
            } else if (cells_safety[p] > guess_safety) {
                {
                    std::lock_guard<std::mutex> lck(to_play_mtx);
                    to_guess = p;
                    guess_safety = cells_safety[p];
                }
            }

        }
    }
}

int main(void) {
    Grid grid;
    std::mutex grid_mtx;
    std::queue<Pos> to_compute;
    Pos to_guess;
    std::mutex to_compute_mtx;
    std::queue<Pos> to_play;
    std::mutex to_play_mtx;

    std::thread my_game_interface(game_interface, std::ref(grid), std::ref(grid_mtx), std::ref(to_compute), std::ref(to_guess), std::ref(to_compute_mtx), std::ref(to_play), std::ref(to_play_mtx));
    std::thread my_player(game_player, std::ref(grid), std::ref(grid_mtx), std::ref(to_compute), std::ref(to_guess), std::ref(to_compute_mtx), std::ref(to_play), std::ref(to_play_mtx));

    my_game_interface.join();
    my_player.join();
    return (0);
}
