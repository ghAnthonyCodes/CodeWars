#include <vector>
#include <set>
#include <cmath>
#include <tuple>
#include <numeric>

using vi = std::vector<int>;
using vvi = std::vector<vi>;
using vsi = std::vector<std::set<int>>;
using vvsi = std::vector<vsi>;

int N = 7;
int total = 0;
vi tb_clues;
vi lr_clues;
vi rl_clues;
vi bt_clues;
std::vector<std::pair<int, int>> parse_order;

void print_grid(vvi &grid) {
  printf("~~~~~~~~~~~~~~~~~~~~~~\n");
  for (int r=0; r<N; r++) {
    for (int c=0; c<N; c++) 
      printf("%2d ", grid[r][c]); 
    printf("\n");
  }
}

int lr_check(vvi &grid, int r) {
  int lr_max = 0;
  int lr = 0;
  for (int i=0; i<N; i++) {
    if (grid[r][i] > lr_max) {
      lr_max = grid[r][i];
      lr++;
    } 
  }
  return lr;
}

int rl_check(vvi &grid, int r) {
  int rl_max = 0;
  int rl = 0;
  for (int i=0; i<N; i++) {
    if (grid[r][N-i-1] > rl_max) {
      rl_max = grid[r][N-i-1];
      rl++;
    } 
  }
  return rl;
}

int tb_check(vvi &grid, int c) {
  int tb_max = 0;
  int tb = 0;
  for (int i=0; i<N; i++) {
    if (grid[i][c] > tb_max) {
      tb_max = grid[i][c];
      tb++;
    } 
  }
  return tb;
}

int bt_check(vvi &grid, int c) {
  int bt_max = 0;
  int bt = 0;
  for (int i=0; i<N; i++) {
    if (grid[N-i-1][c] > bt_max) {
      bt_max = grid[N-i-1][c];
      bt++;
    } 
  }
  return bt;
}

bool dive(int parse_idx, vvi &grid, vsi &rows, vsi &cols) {
  total++;
  
  if (parse_idx == 49 || total > 10000000)
    return true;
  
  int r = parse_order[parse_idx].first;
  int c = parse_order[parse_idx].second;
  
  vi v_intersection;
  std::set_intersection(rows[r].begin(), rows[r].end(), cols[c].begin(), cols[c].end(), std::back_inserter(v_intersection));
  for (auto v : v_intersection) {
    
    // Propose grid[r][c] = v
    grid[r][c] = v;
    
    // If we have filled a row
    int row_sum = std::accumulate(grid[r].begin(), grid[r].end(), 0);
    if (row_sum == 28) {
      if (lr_clues[r] != 0 && lr_check(grid, r) != lr_clues[r]) {
        grid[r][c] = 0;
        continue;
      }
      if (rl_clues[r] != 0 && rl_check(grid, r) != rl_clues[r]) {
        grid[r][c] = 0;
        continue;
      }
    }
    
    // If we have filled a col
    int col_sum = 0;
    for (int x=0; x<N; x++)
      col_sum += grid[x][c];
    
    if (col_sum == 28) {
      if (tb_clues[c] != 0 && tb_check(grid, c) != tb_clues[c]) {
        grid[r][c] = 0;
        continue;
      }
      if (bt_clues[c] != 0 && bt_check(grid, c) != bt_clues[c]) {
        grid[r][c] = 0;
        continue;
      }
    }
    
    rows[r].erase(v);
    cols[c].erase(v);
    if (dive(parse_idx + 1, grid, rows, cols))
      return true;
    grid[r][c] = 0;
    rows[r].insert(v);
    cols[c].insert(v);
  }
  return false;
}

std::vector<std::vector<int>> SolvePuzzle(const std::vector<int> &clues) {

  vsi rows = vsi(N, std::set<int>({ 1, 2, 3, 4, 5, 6, 7 }));
  vsi cols = vsi(N, std::set<int>({ 1, 2, 3, 4, 5, 6, 7 }));
  tb_clues = vi(clues.begin(), clues.begin() + 7);
  rl_clues = vi(clues.begin() + 7, clues.begin() + 14);
  bt_clues = vi(clues.begin() + 14, clues.begin() + 21);
  lr_clues = vi(clues.begin() + 21, clues.begin() + 28);
  std::reverse(bt_clues.begin(), bt_clues.end());
  std::reverse(lr_clues.begin(), lr_clues.end());
  
  // Figure out the parse order
  std::vector<std::tuple<char, int, int>> G;
  for (int i=0; i<N; i++) {
    int r_sq = (int)std::pow(rl_clues[i], 2) + std::pow(lr_clues[i], 2);
    int c_sq = (int)std::pow(tb_clues[i], 2) + std::pow(bt_clues[i], 2);
    G.push_back({ 'R', r_sq, i });
    G.push_back({ 'C', c_sq, i });
  }
  
  std::sort(G.begin(), G.end(), [](const std::tuple<char, int, int> &A, const std::tuple<char, int, int> &B) {
    return std::get<1>(B) < std::get<1>(A); 
  });
  
  // Create parse order
  vvi grid_parse_order = vvi(7, vi(7, -1));
  parse_order = std::vector<std::pair<int, int>>();
  total = 0;
  int parse_idx = 0;
  for (auto &t : G) {
    // If next is row 
    if (std::get<0>(t) == 'R') {
      for (int i=0; i<N; i++) {
        int r = std::get<2>(t);
        if (grid_parse_order[r][i] == -1) {
          grid_parse_order[r][i] = parse_idx++;
          parse_order.push_back({ r, i });
        } 
      }
    } else {
      for (int i=0; i<N; i++) {
        int c = std::get<2>(t);
        if (grid_parse_order[i][c] == -1) {
          grid_parse_order[i][c] = parse_idx++;
          parse_order.push_back({ i, c });
        } 
      }
    }
  }
  
  // How long does it take to generate every solution of a nxn?
  vvi grid = vvi(N, vi(N, 0));
  dive(0, grid, rows, cols);
  print_grid(grid);
  printf("Total: %d\n", total);
  
  return grid;
}
