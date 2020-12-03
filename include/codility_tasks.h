
#ifndef _CODILITY_TASKS_H_
#define _CODILITY_TASKS_H_

#include <vector>
#include <unordered_map>
#include <cmath>

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

}

#endif // _CODILITY_TASKS_H_
