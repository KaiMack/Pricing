#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>
#include <numeric> 
#include <stack>

struct Cell {int row, col;};

class Maze{
    public: 
    Maze(int rows, int cols) : rows(rows), cols(cols), grid(2*rows+1, std::vector<char>(2 * cols + 1, '#')){
        visited.resize(rows, std::vector<bool>(cols, false));
        generate_maze();
    }
    
    void print_maze(){
        for (const auto& row : grid){
            for (char cell: row)
            std::cout << cell; 
        std::cout << '\n';
        }
    }
    private:
    int rows, cols;
    std::vector<std::vector<char>> grid;
    std::vector<std::vector<bool>> visited;
    std::mt19937 rng{std::random_device{}()};
    const std::vector<std::pair<int, int>> directions{{-1,0}, {1,0},{0,-1},{0,1}};
    void generate_maze(){
        std::stack<Cell> stack;
        stack.push({0,0});
        visited[0][0] = true;
        grid[1][1] = ' ';

        while (!stack.empty()){
            Cell current = stack.top();
            std::vector<Cell> neighbors;
            for (const auto& [dr, dc]: directions){
                int nr = current.row + dr;
                int nc = current.col + dc;
                if (is_valid(nr, nc))
                neighbors.push_back({nr, nc});
            }
            if (!neighbors.empty()){
                std::shuffle(neighbors.begin(), neighbors.end(), rng);
                Cell next = neighbors.front();

                int wall_r = current.row + next.row + 1;
                int wall_c = current.col + next.col + 1;
                grid[wall_r][wall_c] = ' ';
                grid[2*next.row+1][2*next.col+1] = ' ';
                visited[next.row][next.col] = true;
                stack.push(next);
            } else{
                stack.pop();
            }
        }
    }
    bool is_valid(int r, int c){
        return r >= 0 && c >= 0 && r < rows && c < cols && !visited[r][c];
    }
};

int main(){
    int rows = 10, cols = 20;
    Maze maze(rows,cols);
    maze.print_maze();
    return 0; 
}

