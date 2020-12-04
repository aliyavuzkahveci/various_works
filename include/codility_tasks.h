
#ifndef _CODILITY_TASKS_H_
#define _CODILITY_TASKS_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>

#include "bitwise_operations.h"

namespace codility_tasks {

void intToBinary(int n, std::vector<int>& binaryList) {
    binaryList.push_back(n % 2);

    if(n/2 == 0) {
        return;
    }

    intToBinary(n/2, binaryList);
}

int binaryGap(int n) {
    std::vector<int> binaryList;
    intToBinary(n, binaryList);

    bool started = false;
    int counter = 0;
    int maxBinaryGap = 0;
    for(const auto iter : binaryList) {
    //for(auto iter = binaryList.begin(); iter != binaryList.end(); ++iter) {
        //if(*iter == 1) {
        if(iter == 1) {
            started = true;
            if(counter > maxBinaryGap) {
                maxBinaryGap = counter;
            }
            counter = 0;
            continue;
        }

        if(started) {
            ++counter;
        }
    }

    return maxBinaryGap;
}

void rotatingArray(std::vector<int> &A, int numOfRotate) {
    if(A.empty()) { // to protect segmentation fault!!!
        return;
    }

    for(int j=0; j<numOfRotate; ++j) {
        //shift to right by 1
        int protectedValue = A[A.size()-1]; // segmentation fault in case of empty!!!
        for(size_t i=0; i<A.size(); ++i) {
            bitwise_operations::swapXOR<int>(protectedValue, A[i]);
        }
    }
}

int oddOccurrenceInArray(std::vector<int> &A) {
    std::unordered_map<int, int> occurenceMap;
    for(auto const & iter : A) {
        occurenceMap[iter]++;
    }

    for(auto const & iter : occurenceMap) {
        if(iter.second % 2 != 0) {
            return iter.first;
        }
    }

    return -1; // should not reach here!!!
}

int missingElement(std::vector<int> &A) { // [1,....,(N+1)]
    int N = A.size();
    // the sum for 1.....N => N*(N+1)/2
    long seriesSum = (long)N*(N+1)/2;
    //since range is [1...(N+1)]
    seriesSum += (N+1);

    for(auto const & iter : A) {
        seriesSum -= iter;
    }

    return (int)seriesSum; // the leftover is the missing element!
}

int tapeEquilibrium(std::vector<int> &A) {
    int preSum = A[0];
    int postSum = 0;
    for(size_t i=1; i<A.size(); ++i) {
        postSum += A[i];
    }

    auto minDiff = std::abs(postSum - preSum);
    for(size_t i=1; i<A.size()-1; ++i) {
        preSum += A[i];
        postSum -= A[i];
        auto tempDiff = std::abs(postSum - preSum);
        if(tempDiff < minDiff) {
            minDiff = tempDiff;
        }
    }

    return minDiff;
}

int frogRiverOne(int X, std::vector<int> &A) {
    std::unordered_set<int> pathSet;

    for(size_t i=0; i<A.size(); ++i) {
        if(A[i] <= X) {
            pathSet.insert(A[i]); // adds each distinct element only once!

            if(pathSet.size() == (size_t)X) { // path accross the river established!
                return i;
            }
        }
    }

    return -1; // if frog does not make it!
}

std::vector<int> maxCounters(int N, std::vector<int> &A) {
    std::vector<int> counters(N, 0); // initialize the list
    int base = 0;
    int maxValue = 0;

    for(size_t i = 0; i<A.size(); ++i) {
        if(A[i] == N+1) {
            base = maxValue;
        } else { // 1 <= iter <= N
            counters[A[i] - 1] = std::max(base, counters[A[i]-1]) + 1;
            maxValue = std::max(maxValue, counters[A[i]-1]);
        }
    }

    std::for_each(counters.begin(), counters.end(), [&](auto & iter){
        if(iter < base) {
            iter = base;
        }
    });

	/* original solution got 77% => failed on performance!!!
    vector<int> counters(N, 0); // initialize the list
    int maxValue = 0;

    for(auto const & iter : A) {
        if(iter == N+1) {
            std::for_each(counters.begin(), counters.end(), [&](auto & iter){
                iter = maxValue;
            });
        } else { // 1 <= iter <= N
            ++counters[iter-1]; // 0-based indexing
            
            if(maxValue < counters[iter-1]) {
                maxValue = counters[iter-1];
            }
        }
    }
    */

    return counters;
}

}

#endif // _CODILITY_TASKS_H_
