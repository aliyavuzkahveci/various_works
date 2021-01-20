
#ifndef _CHESS_KNIGHT_PROBLEM_H_
#define _CHESS_KNIGHT_PROBLEM_H_

#include <iostream>
#include <queue>

namespace chess_knight_problem {

/*
https://www.geeksforgeeks.org/minimum-steps-reach-target-knight/
Given a square chessboard of N x N size, the position of Knight and position of a target is given. 
We need to find out the minimum steps a Knight will take to reach the target position.

Approach:
This problem can be seen as shortest path in unweighted graph.
Therefore we use BFS to solve this problem.
We try all 8 possible positions where a Knight can reach from its position.
If reachable position is not already visited and is inside the board, we push this state into queue with distance 1 more than its parent state.
Finally we return distance of target position, when it gets pop out from queue.

Below code implements BFS for searching through cells, where each cell contains its coordinate and distance from starting node.
In worst case, below code visits all cells of board, making worst-case time complexity as O(N^2)
*/

// structure for storing a cell's data 
struct cell { 
    int x, y; 
    int dis;

    cell() {} 

    cell(int x, int y, int dis) 
        : x(x), y(y), dis(dis) 
    {}
};


// Utility method returns true if (x, y) lies inside Board 
bool isInside(int x, int y, int N) {
    if (x >= 1 && x <= N && y >= 1 && y <= N) {
        return true; 
    } else {
        return false; 
    }
}

// Method returns minimum step 
// to reach target position 
int minStepToReachTarget(int knightPos[], int targetPos[], int N) {
    // x and y direction, where a knight can move 
    // (relative to the current position)
    int dx[] = { -2, -1, 1, 2, -2, -1, 1, 2 };
    int dy[] = { -1, -2, -2, -1, 1, 2, 2, 1 };
    auto targetSize = 8;
  
    // queue for storing states of knight in board
    std::queue<cell> q;
  
    // push starting position of knight with 0 distance
    q.push(cell(knightPos[0], knightPos[1], 0));
  
    bool visit[N + 1][N + 1]; // 1-indexed search will be performed!

    // initially, make all cell unvisited
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            visit[i][j] = false;
        }
    }

    // visit starting state
    visit[knightPos[0]][knightPos[1]] = true;

    cell t;
    int x, y;

    // loop untill we have one element in queue
    while (!q.empty()) {
        t = q.front();
        q.pop();

        // if current cell is equal to target cell,return its distance
        if (t.x == targetPos[0] && t.y == targetPos[1]) {
            return t.dis;
        }

        // loop for all reachable states
        for (int i = 0; i < targetSize; i++) {
            x = t.x + dx[i];
            y = t.y + dy[i];

            // If reachable state is not yet visited and inside board, push that state into queue
            if (isInside(x, y, N) && !visit[x][y]) {
                visit[x][y] = true;
                q.push(cell(x, y, t.dis + 1));
            }
        }
    }

    return -1;
}

#pragma region Operations

void minStepToReachTarget_trials() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Chess Knight Problem operations..." << std::endl;
    
    int N = 30;
    int knightPos[] = {1, 1};
    int targetPos[] = {30, 30};
    auto result = minStepToReachTarget(knightPos, targetPos, N);
    std::cout << "Knight Position: (" << knightPos[0] << "," << knightPos[1] << ") & Target: (" << targetPos[0] << "," << targetPos[1] << ")" << std::endl;
    std::cout << "Minimum Number of Steps to Reach Target: " << result << std::endl;

    std::cout << "Ending Chess Knight Problem operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

#pragma endregion Operations

}

#endif // _CHESS_KNIGHT_PROBLEM_H_
