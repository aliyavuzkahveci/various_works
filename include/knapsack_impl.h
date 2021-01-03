
#ifndef _KNAPSACK_IMPL_H_
#define _KNAPSACK_IMPL_H_

#include <vector>
#include <iostream>
#include <algorithm>

namespace knapsack {

/* Method-1 START*/
// recursive approach
// Time complexity: O(2^n)
// Space complexity: O(1)
int knapsack(int W, std::vector<int>& weight, std::vector<int>& value, int size) {
    if(W == 0 || size == 0) {
        return 0;
    }

    // if the last item's weight is more than the knapsack capacity W
    // we shall exclude it since it cannot be in the optimal solution
    if(weight[size-1] > W) {
        return knapsack(W, weight, value, size-1);
    }

    // check the two cases:
    // 1. include the last item
    // 2. exclude the last item
    // choose the maximum of two
    else {
        auto last_included = value[size-1] + knapsack(W-weight[size-1], weight, value, size-1);
        auto last_excluded = knapsack(W, weight, value, size-1);
        return std::max(last_included, last_excluded);
    }

}
/* Method-1 END*/

/* Method-2 START*/
// Dynamic programming (Bottom-Up) based solution
// Time complexity: O(N*M)
// Space complexity: O(N*M) => N: size, M: knapsack capacity W
int knapsack_2(int W, std::vector<int>& weight, std::vector<int>& value, int size) {
    // initialize capacity table with 0s
    std::vector<std::vector<int>> K(size+1, std::vector(W+1, 0)); // dimensions: size+1 x W+1

    // build the capacity table in a bottom-up manner
    for(int i=0; i<=size; ++i) {
        for(int j=0; j<=W; ++j) {
            if(i==0 || j==0) {
                K[i][j] = 0;
            } else if(weight[i-1] <= j) { // compare the two cases and choose the maximum
                auto included = value[i] + K[i-1][j-weight[i-1]];
                K[i][j] = std::max(included, K[i-1][j]);
            } else { // do not include 
                K[i][j] = K[i-1][j];
            }
        }
    }

    return K[size][W];
}
/* Method-2 END*/

/* Method-3 START*/
// Uses Memorization Technique (an extension of recursive technique)
// Dynamic programming (Top-Down)

int knapsackRecursive(int W, std::vector<int>& weight, std::vector<int>& value, int index, std::vector<std::vector<int>>& K) {
    if(index < 0) { // base condition
        return 0;
    } else if(K[index][W] != -1) { // if it was alread calculated, return immediately!
        return K[index][W];
    } else if(weight[index] > W) { // first calculate & store then return the value
        K[index][W] = knapsackRecursive(W, weight, value, index-1, K);
        return K[index][W];
    } else { // choose the maximum value in case the weight of the item is less than W
        auto included = value[index] + knapsackRecursive(W-weight[index], weight, value, index-1, K);
        auto excluded = knapsackRecursive(W, weight, value, index-1, K);
        // store the value first
        K[index][W] = std::max(included, excluded);
        // return the value
        return K[index][W];
    }
}

// Time complexity: O(N*M)
// Space complexity: O(N*M) => N: size, M: knapsack capacity W
int knapsack_3(int W, std::vector<int>& weight, std::vector<int>& value, int size) {
    // initialize state map with -1
    std::vector<std::vector<int>> K(size, std::vector<int>(W+1, -1));
    
    return knapsackRecursive(W, weight, value, size-1, K);
}

/* Method-3 END*/

void knapsack_trials() {
    auto W = 50;
    std::vector<int> value = {60, 100, 120};
    std::vector<int> weight = {10, 20, 30};
    auto size = 3;

    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Knapsack trials for parameters W=50, weight={10, 20, 30}, value={60, 100, 120}, size=3" << std::endl;
    std::cout << "Method1: knapsack() = " << knapsack(W, weight, value, size) << std::endl;
    std::cout << "Method2: knapsack_2() = " << knapsack_2(W, weight, value, size) << std::endl;
    std::cout << "Method3: knapsack_3() = " << knapsack_3(W, weight, value, size) << std::endl;
    std::cout << "=======================================" << std::endl;
}

} // knapsack

#endif // _KNAPSACK_IMPL_H_
