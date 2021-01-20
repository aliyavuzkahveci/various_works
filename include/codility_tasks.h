
#ifndef _CODILITY_TASKS_H_
#define _CODILITY_TASKS_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <stack>
#include <climits>
#include <set>

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

int missingInteger(std::vector<int> &A) {
    int smallest = 1;

    std::sort(A.begin(), A.end());

    for(auto const & iter : A) {
        if(smallest == iter) {
            ++smallest;
        }
    }

    return smallest;
}

int permutationCheck(std::vector<int> &A) {
    // return 1 => permutation OK, 0 => permutation FAIL
    /*
    // result = 75% (83% correctness, 66% performance)
    int N = A.size();
    long total = (long) N * (N + 1) / 2;

    for(auto const & iter : A) {
        total -= iter;
    }

    return total == 0 ? 1 : 0;
    */
    int maxNumPossible = A.size() + 1;
    // initialize list with 0's
    std::vector<int> controlList(maxNumPossible, 0);

    for(auto const & iter : A) {
        if(iter > maxNumPossible) {
            // number shall not exceed maximum number which is A.size() + 1
            return 0;
        }

        ++controlList[iter];
        if(controlList[iter] != 1) {
            // number encountered in the list before!
            return 0;
        }
    }

    return 1;
}

int countDiv(int A, int B, int K) {
    int inclusive = (A%K) == 0 ? 1 : 0;

    return (B/K) - (A/K) + inclusive;
}

std::vector<int> genomicRangeQuery(std::string &S, std::vector<int> &P, std::vector<int> &Q) {
    /*
    // Total Score: %62: Correctness: 100% - Performance: 0%
    //start control with the minimal impact factor i.e. A
    // A: 1, C: 2, G: 3, T: 4
    std::vector<char> nucleotideList = {'A', 'C', 'G', 'T'};

    std::vector<int> impactFactorList;
    for(size_t i=0; i<P.size(); ++i) {
        int start = P[i];
        int end = Q[i];

        auto subStr = S.substr(start, end-start+1);

        for(size_t i=0; i<nucleotideList.size(); ++i) {
            if(subStr.find(nucleotideList[i]) != std::string::npos) {
                impactFactorList.push_back(i+1);
                break;
            }
        }
    }

    return impactFactorList;
    */

   // Total Score: 100%
   //holding the cumulative sums for each nucleotide
    std::vector<int> cum_sum_A(S.size());
    std::vector<int> cum_sum_C(S.size());
    std::vector<int> cum_sum_G(S.size());
    //vector<int> cum_sum_T(S.size()); // not really necessary since it falls into ELSE case!

    int A_counter = 0;
    int C_counter = 0;
    int G_counter = 0;
    //int T_counter = 0; // not really necessary since it falls into ELSE case!

    // calculate the cumulative sums
    for(size_t i=0; i<S.size(); ++i) {
        if(S[i] == 'A') {
            ++A_counter;
        } else if(S[i] == 'C') {
            ++C_counter;
        } else if(S[i] == 'G') {
            ++G_counter;
        }/* else if(S[i] == 'T') { //not really necessary since it falls into ELSE case!
            ++T_counter;
        }*/

        // cumulative is the key to solve this problem
        cum_sum_A[i] = A_counter;
        cum_sum_C[i] = C_counter;
        cum_sum_G[i] = G_counter;
        //cum_sum_T.push_back(T_counter); // not really necessary since it falls into ELSE case!
    }

    std::vector<int> result(P.size());
    for(size_t i=0; i<P.size(); ++i) {
        int A_at_start = (S[P[i]] == 'A') ? 1 : 0;
        int C_at_start = (S[P[i]] == 'C') ? 1 : 0;
        int G_at_start = (S[P[i]] == 'G') ? 1 : 0;
        //int T_at_start = S[i] == 'T' ? 1 : 0; // not really necessary since it falls into ELSE case!

        if(cum_sum_A[Q[i]] - cum_sum_A[P[i]] + A_at_start > 0) {
            result[i] = 1; // A exists in the sub range!
        } else if(cum_sum_C[Q[i]] - cum_sum_C[P[i]] + C_at_start > 0) {
            result[i] = 2; // C exists in the sub range!
        } else if(cum_sum_G[Q[i]] - cum_sum_G[P[i]] + G_at_start > 0) {
            result[i] = 3; // G exists in the sub range!
        } else { // if none of the above nucleotides exists in   the range
            result[i] = 4; // T exists in the sub range
        }
    }

    return result;
}

int minAvgTwoSlice(std::vector<int> &A) {
    if(A.size() < 3) {
        return 0; // the min slice length must be 2
    }

    int min_2_slice_index = 0;
    int min_3_slice_index = 0;
    int min_2_slice_sum = A[0] + A[1];
    int min_3_slice_sum = min_2_slice_sum + A[2];

    for(size_t i=2; i<A.size(); ++i) {
        int temp = A[i-1] + A[i];
        if(temp < min_2_slice_sum) {
            min_2_slice_sum = temp;
            min_2_slice_index = i-1;
        }

        temp += A[i-2];
        if(temp < min_3_slice_sum) {
            min_3_slice_sum = temp;
            min_3_slice_index = i-2;
        }
    }

    // instead of dividing into 2 & 3
    // we are multplying to their smallest common multiple which is 6
    min_2_slice_sum *= 3;
    min_3_slice_sum *= 2;

    if(min_2_slice_sum == min_3_slice_sum) {
        // return the smallest index in case their averages are same!
        if(min_3_slice_index < min_2_slice_index) return min_3_slice_index;
        else return min_2_slice_index;
    } else if (min_3_slice_sum < min_2_slice_sum) {
        return min_3_slice_index;
    } else {
        return min_2_slice_index;
    }

/*	// Total Score: 90% - Correctness: 80%, Performance: 100%
    double minAvg = ((double)(A[0] + A[1])) / 2;
    double tempAvg;
    int minSliceIndex = 0;
    for(size_t i=0; i<A.size()-2; ++i) {
        // check for 2-element slices
        tempAvg =((double)(A[i] + A[i+1])) / 2;
        if(tempAvg < minAvg) {
            minAvg = tempAvg;
            minSliceIndex = i;
        }

        // check for 3-element slices
        tempAvg = ((double)(A[i] + A[i+1] + A[i+2])) / 3;
        if(tempAvg < minAvg) {
            minAvg = tempAvg;
            minSliceIndex = i;
        }
    }

    //controlling the slice with last two elements
    tempAvg = ((double)(A[A.size()-1] + A.size()-2)) / 2;
    if(tempAvg < minAvg) {
        minSliceIndex = A.size() - 2;
    }

    return minSliceIndex;
*/
}

int passingCars(std::vector<int> &A) {
    int passingCarCounter = 0;

    if(A.size() == 1) {
        return passingCarCounter;
    }

    int numOf0s = 0;
    for(size_t i=0; i<A.size(); ++i) {
        if(A[i] == 0) {
            ++numOf0s;
        } else {
            passingCarCounter += numOf0s;
            if(passingCarCounter > 1000000000) {
                return -1;
            }
        }
    }

    return passingCarCounter;
	/*
	// Total Score: 60% - Correctness: 100% & Performance: 20%
    int passingCarCounter = 0;

    if(A.size() == 1) {
        return passingCarCounter;
    }

    for(size_t i=0; i<A.size()-1; ++i) {
        if(A[i] != 0) continue;
        for(size_t j=i+1; j<A.size(); ++j) {
            // if enter here, A[i] = 0 !
            passingCarCounter += (A[i] ^ A[j]);
            
            if(passingCarCounter > 1000000000) {
                return -1;
            }
        }
    }

    return passingCarCounter;
	*/
}

int distinct(std::vector<int> &A) {
    std::unordered_set<int> container;

    for(auto const & iter : A) {
        container.insert(iter);
    }

    return container.size();
}

int maxProductOfThree(std::vector<int> &A) {
    std::sort(A.begin(), A.end());
    //we are given the size will be at least 3!
    auto size = A.size();
    auto maxProduct = A[size-1] * A[size-2] * A[size-3];

    //there are might be at least 2 negative numbers
    // which their product becomes positive!!
    if(A[0] < 0 && A[1] < 0) {
        auto temp = A[0] * A[1] * A[size-1];
        maxProduct = (temp > maxProduct) ? temp : maxProduct;
    }

    return maxProduct;
}

int numberOfDiscIntersections(std::vector<int> &A) {
    if(A.size() < 2) {
        return 0;
    }

    std::vector<std::pair<long, long>> sortedDiscRanges;

    for(size_t i=0; i<A.size(); ++i) {
        long start = (long)i - A[i];
        long end = (long)i + A[i];
        sortedDiscRanges.push_back({start, end});
    }

    std::sort(sortedDiscRanges.begin(), sortedDiscRanges.end(), 
    [&](std::pair<long, long> p1, std::pair<long, long> p2){
        return p1.first < p2.first;
    });

    int intersectNum = 0;
    for(size_t i=0; i<sortedDiscRanges.size(); ++i) {
        auto endPoint = sortedDiscRanges[i].second;

        for(size_t j=i+1; j<sortedDiscRanges.size() && sortedDiscRanges[j].first <=endPoint; ++j) {
            ++intersectNum;

            if(intersectNum > 10000000) {
                return -1;
            }
        }
    }

    return intersectNum;
	/*
    // Total Score: 81% - Correctness: 100% & Performance: 62%
    if(A.size() < 2) {
        return 0;
    }
	
    int intersectNum = 0;

    for(size_t i=0; i<A.size()-1; ++i) {
        long startOfDisc = (long)i - A[i];
        long endOfDisc = (long)i + A[i];
        for(size_t j=i+1; j<A.size(); ++j) {
            long start = (long)j - A[j];
            long end = (long)j + A[j];

            if(end < startOfDisc || start > endOfDisc) {
                continue; // no intersection
            } else {
                ++intersectNum;

                // since we count A-B as well as B-A
                if(intersectNum > 10000000) {
                    return -1; // exceeded the max disc intersection
                }
            }
        }
    }

    return intersectNum;
	*/
}

bool triangularityCheck(int A, int B, int C) {
    long sum1 = (long)A + B;
    long sum2 = (long)A + C;
    long sum3 = (long)B + C;
    if(sum1 > C && sum2 > B && sum3 > A) {
        return true; // triangular!
    } else {
        return false; // NOT triangular!
    }
}

int triangle(std::vector<int> &A) {
    auto size = A.size();
    if(size < 3) {
        return 0; // list should have at least 3 elements!
    }

    std::sort(A.begin(), A.end());
    auto midPoint = size / 2;

    if((triangularityCheck(A[midPoint-1], A[midPoint], A[midPoint+1])) || 
       (triangularityCheck(A[0], A[1], A[2])) || 
       (triangularityCheck(A[size-1], A[size-2], A[size-3]))) {
        return 1;
    } else {
        return 0;
    }
}

int brackets(std::string &S) {
    //return 1 => properly nested
    std::stack<char> properlyStack;

    for(size_t i=0; i<S.size(); ++i) {
        if(S[i] == '(' || S[i] == '[' || S[i] == '{') {
            properlyStack.push(S[i]);
        } else if(S[i] == ')' || S[i] == ']' || S[i] == '}') {
            auto checkChar = (S[i] == ')' ? '(' : (S[i] == ']' ? '[' : '{'));
            if(properlyStack.top() == checkChar) {
                properlyStack.pop(); // remove from the stack
                continue;
            } else {
                return 0; // improperly nested!
            }
        }
    }

    return properlyStack.empty() ? 1 : 0;
}

int fish(std::vector<int> &A, std::vector<int> &B) {
    std::stack<int> upFish; // 0: flowing upstream
    std::stack<int> downFish; // 1: flowing downstream
    // meet criteria: P < Q && B[P] = 1 && B[Q] = 0
    //B[i-1] = 1 && B[i] = 0

    // sizes of A and B are the same!!!
    for(size_t i=0; i<A.size(); ++i) {
        if(B[i] == 0) { // upstream
            // upstream fish will eat all the downstream SMALLER ones!!!
            while(downFish.size() && downFish.top() < A[i]) {
                downFish.pop(); // fish flowing downstream dies!
            }

            if(downFish.empty()) {
                //std::cout << "adding " << A[i] << " to UPSTREAM stack" << std::endl;
                upFish.push(A[i]);
            }// else {
                // fish at A[i] dies!
                // fish flowing downstream lives!
            //}
        } else { // B[i] = 1 {downstream}
           downFish.push(A[i]);
        }
    }

    // either one or both of team could be empty!
    return upFish.size() + downFish.size();
}

int nested(std::string &S) {
    // S contains only '(' or ')'
    // return 1 => properly nested
    if(S.empty()) {
        return 1;
    }

    std::stack<char> nestedParaStack;
    for(size_t i=0; i<S.size(); ++i) {
        if(S[i] == '(') {
            nestedParaStack.push(S[i]);
        } else if(nestedParaStack.empty()) {
            return 0; // not nested properly!
        } else {
            nestedParaStack.pop(); // remove one '('
        }
    }

    // stack should be empty at the end to be properly nested!
    return nestedParaStack.empty() ? 1 : 0;
}

int stoneWall(std::vector<int> &H) {	
	std::stack<int> s;
    int numOfBlocks = 0;
    for(auto const & iter : H) {
        while (s.size() && s.top() > iter) {
            s.pop();
        }
        
        if (s.size() && iter == s.top()) {
            continue;
        } else {
            s.push(iter);
            ++numOfBlocks;
        }
    }

    // Total Score: 35% - Correctness: 80% & Performance: 11%
    /*int numOfBlocks = 0;
    std::queue<int> leftOvers;
    int currentMinHeight = INT_MAX;

    for(auto const & iter : H) {
        if(iter < currentMinHeight) {
            currentMinHeight = iter;
            ++numOfBlocks;

            //std::cout << "add block for " << iter << std::endl;

            int localMax = 0;
            int localMin = INT_MAX;
            while(leftOvers.size()) { //continue until empty!
                if(leftOvers.front() < localMin) {
                    localMin = leftOvers.front();
                }

                if(leftOvers.front() > localMax) {
                    //std::cout << "WHILE add block for " << leftOvers.front() << std::endl;
                    localMax = leftOvers.front();
                    ++numOfBlocks;
                } else if(leftOvers.front() < localMax && leftOvers.front() > localMin) {
                    //std::cout << "WHILE add block for " << leftOvers.front() << std::endl;
                    ++numOfBlocks;
                //} else {
                //    std::cout << "WHILE NOT add block " << leftOvers.front() << std::endl;
                }
                leftOvers.pop();
            }

        } else {
            //std::cout << "NOT add block " << iter << std::endl;
            auto leftOver = iter - currentMinHeight;
            if(leftOver > 0) {
                //std::cout << "add to Queue " << leftOver <<  " for " << iter << std::endl;
                leftOvers.push(leftOver);
            //} else {
            //    std::cout << "Skipping both for " << iter << std::endl;
            }
        }
    }

    //std::cout << "numOfBlocks: " << numOfBlocks << std::endl;
    //std::cout << "leftOvers.size(): " << leftOvers.size() << std::endl;

    int localMax = 0;
    int localMin = INT_MAX;
    while(leftOvers.size()) { //continue until empty!
        if(leftOvers.front() < localMin) {
            localMin = leftOvers.front();
        }

        if(leftOvers.front() > localMax) {
            //std::cout << "WHILE add block for " << leftOvers.front() << std::endl;
            localMax = leftOvers.front();
            ++numOfBlocks;
        } else if(leftOvers.front() < localMax && leftOvers.front() > localMin) {
            //std::cout << "WHILE add block for " << leftOvers.front() << std::endl;
            ++numOfBlocks;
        //} else {
        //    std::cout << "WHILE NOT add block " << leftOvers.front() << std::endl;
        }
        leftOvers.pop();
    }*/

    return numOfBlocks;
}

int dominator(std::vector<int> &A) {
    std::unordered_map<int, int> occurrenceMap;

    int leastDomCount = A.size() / 2 + 1;

    for(size_t i=0; i<A.size(); ++i) {
        occurrenceMap[A[i]]++;
        if(occurrenceMap[A[i]] >= leastDomCount) {
            return i; //there could be only one!
        }
    }

    return -1; // if reach here, no dominator!
}

int equiLeader(std::vector<int> &A) {
    int missing = -1000000001; // range: [-1,000,000,000..1,000,000,000]
    int numOfEquiLeaders = 0;
    int lastLeader = missing;
    int length = A.size();

    if(length < 2) {
        return numOfEquiLeaders; // 0
    }

    std::vector<int> subArrayLeaders(length, missing);
    std::unordered_map<int, int> occurrenceMap;

    for (int i = 0; i < length; ++i)
    {
        ++occurrenceMap[A[i]];

        if (occurrenceMap[A[i]] > (i + 1) / 2) {
            lastLeader = subArrayLeaders[i] = A[i];
        } else if (occurrenceMap[lastLeader] > (i + 1) / 2) {
            subArrayLeaders[i] = lastLeader;
        }
    }

    occurrenceMap.clear();
    lastLeader = missing;
    for (int i = length - 1; i > 0; --i)
    {
        ++occurrenceMap[A[i]];

        // starting from back, comparing leaders from both sub arrays
        if (occurrenceMap[A[i]] > (length - i) / 2 && subArrayLeaders[i - 1] == A[i]) {
            lastLeader = A[i];
            ++numOfEquiLeaders;
        } else if (occurrenceMap[lastLeader] > (length - i) / 2 && subArrayLeaders[i - 1] == lastLeader) {
            ++numOfEquiLeaders;
        }
    }

    return numOfEquiLeaders;
	/*
	// Total Score: 88% - Correctness: 100% - Performance: 75%
    int lowLeader, highLeader, numOfEquiLeaders = 0;
    std::unordered_map<int, int> lowOccurrenceMap, highOccurrenceMap;

    ++lowOccurrenceMap[A[0]];
    lowLeader = A[0];

    highLeader = 1000000001; // simply out of available range

    // get high leader
    int minNumOfOccurrence = ((A.size()-1) / 2) + 1;
    for(size_t i=1; i<A.size(); ++i) {
        ++highOccurrenceMap[A[i]];
        if(highOccurrenceMap[A[i]] >= minNumOfOccurrence) {
            highLeader = A[i];
        }
    }

    if(lowLeader == highLeader) {
        ++numOfEquiLeaders;
    }

    // iterate through vector by inserting into low and removing from high
    for(size_t i=1; i<A.size()-1; ++i) {
        ++lowOccurrenceMap[A[i]];
        --highOccurrenceMap[A[i]];
        if(highOccurrenceMap[A[i]] == 0) {
            highOccurrenceMap.erase(A[i]);
        }

        // low part's size = i+1
        // high part's size = A.size()- i-1
        size_t lowSize = i+1;
        size_t highSize = A.size()-i-1;

        //number of distinct values in the amp is an indication 
        //for a leader in the array
        if(lowOccurrenceMap.size() <= lowSize && 
           highOccurrenceMap.size() <= highSize) {
            //find leaders
            int minOccurForLow = (lowSize / 2) + 1;
            int minOccurForHigh = (highSize / 2) + 1;

            lowLeader = -1;
            highLeader = -2;

            for(auto const &iter : lowOccurrenceMap) {
                if(iter.second >= minOccurForLow) {
                    lowLeader = iter.first;
                    break;
                }
            }

            for(auto const &iter : highOccurrenceMap) {
                if(iter.second >= minOccurForHigh) {
                    highLeader = iter.first;
                    break;
                }
            }

            if(lowLeader == highLeader) {
                ++numOfEquiLeaders;
            }
        }
    }

    return numOfEquiLeaders;
	*/
}

int maxDoubleSliceSum(std::vector<int> &A) {
    // 3 <= A.size() <= 100,000
    // range: [-10,000, 10,000]
    // 0 <= X < Y < Z < N
    auto length = A.size();
    std::vector<int> firstSliceSum(length, 0);
    std::vector<int> secondSliceSum(length, 0);

    int max_sum = 0;
    int current_sum = 0;
    int min_sum = 0;

    // by definition, FIRST and last elements cannot be involved!
    // compute BEST sum among values assume Y=i !
    for(size_t i=1; i<length-1; ++i) {
        firstSliceSum[i] = std::max(0, current_sum - min_sum);
        current_sum += A[i];
        min_sum = std::min(current_sum, min_sum);
    }

    current_sum = 0;
    min_sum = 0;
    // by definition, first and LAST elements cannot be involved!
    // compute BEST sum among values assume Y=i !
    for(size_t i=length-2; i>0; --i) {
        secondSliceSum[i] = std::max(0, current_sum - min_sum);
        current_sum += A[i];
        min_sum = std::min(current_sum, min_sum);
    }

    // find the maximal sum of double slices according to y position!
    for(size_t i=1; i<length-1; ++i) {
        max_sum = std::max(max_sum, firstSliceSum[i] + secondSliceSum[i]);
    }

    return max_sum;
}

int maxProfit(std::vector<int> &A) {
    // N: [0, 400,000]
    // 0 <=  A[i] <= 200,000

    if(std::is_sorted(A.rbegin(), A.rend())) {
        // share descreases each day!
        return 0;
    }

    int min = A[0];
    int max = A[0];
    int max_profit = 0;
    for(size_t i=0; i<A.size(); ++i) {
        if(A[i] < min) {
            min = max = A[i];
        } else if(A[i] > max) {
            max = A[i];
        }

        auto profit = max - min; 
        max_profit = std::max(max_profit, profit);
    }

    return max_profit;
}

int maxSliceSum(std::vector<int> &A) {
    // N: [1, 1,000,000]
    // -2,147,483,648 < A[i] < 2,147,483,648
    long max_slice_ending = (long)0;
    long max_sum = (long)INT_MIN;
    for(auto const & iter : A) {
        // max sum of slice ending at iter (might start anywhere)!
        max_slice_ending = std::max(max_slice_ending + iter, (long)iter);

        // keep the maximum sum!
        max_sum = std::max(max_sum, max_slice_ending);
    }   

    return max_sum;
}

int countFactors(int N) {
    // N = D * M => 24 = 6 * 4
    // count until root of N, add 2 for each divisor
    // since the other divisor is bigger than root of N
    // which we won't visit
    int numOfDivisor = 0;
    int i;
    auto sq_root = std::sqrt(N); // don't have to count multiply on each loop
    for(i=1; i<sq_root; ++i) {
    //for(i=1; i*i<N; ++i) { // condition: i is less than root of N
        if(N % i == 0) {
            numOfDivisor += 2; // because N = i * j => count also j
        }
    }
    if(i*i == N) { // N is an exact square of a number
        numOfDivisor += 1; // because N = i * i => don't count the same number twice
    }

    return numOfDivisor;
}

bool placementTrial(int candidateFlagCount, std::vector<int> peakList) {
    if(candidateFlagCount == 1) {
        // it is always possible to put 1 flag!
        return true;
    }

    auto prevPeakPos = peakList[0];
    auto remainingFlags = candidateFlagCount - 1; // first flag into the first peak
    for(size_t i=1; i<peakList.size(); ++i) {
        auto distance = peakList[i] - prevPeakPos;
        if(distance >= candidateFlagCount) {
            // we can put the flag to this peak!
            --remainingFlags;

            if(remainingFlags == 0) {
                // all the flags are put!
                return true;
            }

            prevPeakPos = peakList[i];
        }
    }

    return false; // failed to put given number of flags to the peaks!
}

int flags(std::vector<int> &A) {
    auto size = A.size();
    if(size < 3) { // there cannot be a peak!
        return 0;
    }

    //the max number of the peeks is always less than N / 2.
    //Indeed, the max number of the peeks, (N - 3) / 2 + 1. 
    std::vector<int> peakList(size/2, 0);
    int numOfPeaks = 0;  

    for(size_t i=1; i<size-1; ++i) {
        if(A[i] > A[i-1] && A[i] > A[i+1]) {
            peakList[numOfPeaks] = i; // store peak positions!
            ++numOfPeaks;
        }
    }

    if(numOfPeaks == 0) { //no peak => no flag 
        return 0;
    }

    //we use the bisection algorithm to find how many flags can be set!
    int maxFlags   = 0;
    int minFlags = 1; //as we have at least one peek, min can be 1.

    //we need to round up the sqrt(size-2) to get upper limit!
    auto temp1 = sqrt(size-2);
    auto temp2 = (int)temp1;

    // extra logic to cover size=3, since we don't want it to be 0
    auto candidateMax = 1; // there must be at least 1 flag!!!
    if(temp2 > 1) {
        candidateMax = (temp1 == (double)temp2) ? (temp2 - 1) : (temp2 + 1);
    }

    // setting the maximum flags possible
    maxFlags = (numOfPeaks < candidateMax) ? numOfPeaks : candidateMax;

    //cout << "maxFlags: " << maxFlags << endl;

    int maxFlagsAllowed = 0;
    while(minFlags <= maxFlags) {
        auto avg = (minFlags + maxFlags) / 2;

        //cout << "maxFlags: " << maxFlags << endl;
        //cout << "minFlags: " << minFlags << endl;
        //cout << "avg: " << avg << endl;

        // check if we can put avg # flags to the peaks!
        auto result = placementTrial(avg, peakList);

        //cout << "result:" << (result ? "True" : "False") << endl;

        if(result) { // successfully put flags to the peaks
            // try a higher number for the next round!
            minFlags = avg + 1;
            maxFlagsAllowed = avg;
        } else { // failed to put flags to the peaks
            // try a smaller number for the next round!
            maxFlags = avg - 1;
        }
    }

    return maxFlagsAllowed;
}

int minPerimeterRectangle(int N) {
    // we have to find the factors of N = A*B
    std::set<int> perimeters; // to preserve ascending order!
    double sqrt = std::sqrt(N);
    int sqrt_casted = (int)sqrt;

    for(int i=1; i<=sqrt_casted; ++i) {
        if(N % i == 0) {
            int otherSide = N / i;
            perimeters.insert(i + otherSide);
        }
    }

    // perimeter = 2 * (A + B)
    return 2 * (*perimeters.begin());
}

int peaks(std::vector<int> &A) {
    auto size = A.size();
    if(size < 3) {
        // there cannot be a peak in such a list!!!
        return 0;
    }
    // since blocks should contain equal number of elements
    // we have to count the dividers of the size of the array.

    std::set<int> dividerList; // we need to keep the list sorted!
    auto sqrt = (size_t)std::sqrt(size);
    for(size_t i=1; i<=sqrt; ++i) {
        if(size % i == 0) {
            dividerList.insert(i);
            dividerList.insert(size / i);
        }
    }
    /////////////////////////////////////////////////////////////////////

    // counting the peak indexes!
    std::vector<int> peakList(size, 0);
    size_t numOfPeaks = 0;
    // A[0] && A[size-1] cannot be peak by convention!!!
    for(size_t i=1; i<size-1; ++i) {
        if(A[i] > A[i-1] && A[i] > A[i+1]) {
            peakList[i] = 1;
            ++numOfPeaks;
        }
    }

    if(numOfPeaks < 2) {
        // if no peak, then result is 0
        // if peak count is 1, there could only be 1 block.
        return numOfPeaks;
    }
    /////////////////////////////////////////////////////////////////////

    // calculate the cumulative sum so that we shall know
    // how many peaks there are  until the given index!
    std::vector<int> peakCountList(size, 0);
    for(size_t i=1; i<size; ++i) {
        peakCountList[i] += peakCountList[i-1] + peakList[i];
    }
    /////////////////////////////////////////////////////////////////////

    // start with the smallest block size i.e. dividerList.begin()
    for(auto const & blockSize : dividerList) {
        bool noPeak = false;
        for(size_t i=0; i<size; i+=blockSize) {
            auto peakCountUntilThisBlock = (i==0) ? 0 : peakCountList[i-1];
            auto peakCountInThisBlock = peakCountList[i+blockSize-1] - peakCountUntilThisBlock;
            if(peakCountInThisBlock == 0) {
                // list cannot be divided with this block size!
                noPeak = true;
                break;
            }
        }

        if(!noPeak) {
            // all the blocks has at least 1 peak in it!
            return size / blockSize; // block count!!!
        }
    }

    return 0; // list could not be divided into blocks!
}

std::vector<int> countNonDivisible(std::vector<int> &A) {
    auto size = A.size();

    auto iter = std::max_element(A.begin(), A.end());
    auto existListSize = *iter + 1; // since 0-indexed :)
    std::vector<int> existList(existListSize, 0); // initially all 0

    for(size_t i=0; i<size; ++i) {
        ++existList[A[i]]; // also holds multiple existence!
    }

    // initially set maximum 
    std::vector<int> nonDivCountList(size, size);
    for(size_t i=0; i<size; ++i) {
        auto element = A[i];
        for(int j=1; (j*j)<=element; ++j) {
            if(element % j == 0) { // A[i] can be divided by j
                // remove dividible elements (includes also multiple existences)
                nonDivCountList[i] -= existList[j];

                // also subtract the other divider
                auto otherDivider = element / j;
                if(otherDivider != j) { // don't count same j again
                    nonDivCountList[i] -= existList[otherDivider];
                }
            }
        }
    }

    return nonDivCountList;
}

std::vector<int> countSemiprimes(int N, std::vector<int> &P, std::vector<int> &Q) {
    // N: [1, 50000] => maximum number in the two lists!
    // M: [1, 30000]
    // P[i] <= Q[i]
    std::vector<int> semiPrimeList(N+1, 0); // N+1: 0-based indexing
    std::vector<int> primesDivisors(N+1, 0);

    for(int i=2; (i*i)<=N; ++i) {
        if(semiPrimeList[i] == 0) { // still considered as prime!
            for(int k=(i*i); k<=N; k+=i) {
                semiPrimeList[k] = i; // set the divider 
            }
        }
    }
/*
    cout << "semiPrimeList:";
    for(int i=0; i<N+1; i++) {
        cout << " " << semiPrimeList[i];
    }
    cout << endl;
*/

    int counter=0;
    std::vector<int> semiPrimeCounterList(N+1, 0);
    for(int i=2; i<=N; ++i) {
        if(semiPrimeList[i] != 0 && semiPrimeList[i / semiPrimeList[i]] == 0) {
            ++counter;
        }
        semiPrimeCounterList[i] = counter;
    }
/*
    cout << "semiPrimeList:";
    for(int i=0; i<N+1; i++) {
        cout << " " << semiPrimeCounterList[i];
    }
    cout << endl;
*/
    auto size = Q.size();
    std::vector<int> result(size, 0);
    for(size_t i=0; i<size; ++i) {
        result[i] = semiPrimeCounterList[Q[i]] - semiPrimeCounterList[P[i]-1];
    }

    return result;
}

int gcd(int a, int b) { // greatest common divisor
    if(a == b) {
        return a;
    } else if(a%2 == 0 && b%2 == 0) {
        return 2 * gcd(a/2, b/2);
    } else if(a%2 == 0) {
        return gcd(a/2, b);
    } else if(b%2 == 0) {
        return gcd(a, b/2);
    } else if(a>b) {
        return gcd(a-b, b);
    } else /*if(b>a)*/{
        return gcd(a, b-a);
    }
}

int chocolatesByNumbers(int N, int M) {
    return N / gcd(N, M);
}

int commonPrimeDivisors(std::vector<int> &A, std::vector<int> &B) {
    auto size = A.size();
    auto returnCounter = 0;

    for(size_t i=0; i<size; i++) {
        auto a = A[i];
        auto b = B[i];
        if(a == b) {
            // if the numbers are the same increase directly!
            ++returnCounter;
            continue;
        }

        auto divisor = gcd(a, b);
        while(divisor != 1) {
            auto a_div = divisor;
            while(a % a_div == 0 && a_div != 1) {
                a /= a_div;

                a_div = gcd(a, a_div);
            }

            auto b_div = divisor;
            while(b % b_div == 0 && b_div != 1) {
                b /= b_div;

                b_div = gcd(b, b_div);
            }

            divisor = gcd(a, b);
        }

        if(a == b) {
            ++returnCounter;
        }
    }

    return returnCounter;
}

int fib_frog(std::vector<int>& A) {
    int size = A.size();

    // calculate the fibonacci numbers 
    // as the last member is less than or equal to size of array A
    std::vector<int> fiboList;
    fiboList.push_back(0); // for 0
    fiboList.push_back(1); // for 1
    auto lastFiboNum = 1; // for 1
    for(int i=2; lastFiboNum<=size; ++i) {
        lastFiboNum = fiboList[i-1] + fiboList[i-2];
        fiboList.push_back(lastFiboNum);
        //cout << "adding " << lastFiboNum << endl;
    }

    std::set<int> positions = {size};
    for(int jumps = 1; jumps <=(size+1); ++jumps) {
        //cout << "jumps: " << jumps << endl;
        std::set<int> newPositions;
        for(auto const & position : positions) {
            //cout << "position: " << position << endl;
            for(auto const & fiboNum : fiboList) {
                //cout << "fiboNum: " << fiboNum << endl;
                
                // return current jumps counter if reach to the start point
                if(position - fiboNum == -1) { // pos - (fiboNum-1) == 0
                    return jumps;
                }

                auto prevPos = position - fiboNum;
                //cout << "prevPos: " << prevPos << endl;

                // we won't calculate the bigger jumps!
                if(prevPos < 0) {
                    break;
                }

                if(prevPos < size && A[prevPos] == 1) {
                    newPositions.insert(prevPos);
                    //cout << "insert to newPositions: " << prevPos << endl;
                }
            }
        }

        //cout << "newPositions.size(): " << newPositions.size() << endl;
        
        // failed to cross the river!
        if(newPositions.size() == 0) {
            return -1;
        }

        positions = newPositions;
    }

    // if reach here, tried all the options but failed to cross the river!
    return -1;
}

std::vector<int> ladder(std::vector<int> &A, std::vector<int> &B) {
    // L: the maximum number in array A
    auto size = A.size();
    
    auto L = A[0]; // L: limit
    for(size_t i=1; i<size; ++i) {
        if(A[i] > L) {
            L = A[i];
        }
    }

    // calculate the fibonacci numbers up until L+2
    std::vector<unsigned long long> fibList(L+2, 0);
    //fibList[0] = 0; // for 0
    fibList[1] = 1; // for 1
    for(int i=2; i<L+2; ++i) {
        fibList[i] = fibList[i-1] + fibList[i-2];
    }

    std::vector<int> result(size, 0);
    for(size_t i=0; i<size; ++i) {
        // number of different ways to climb is
        // the next fibonacci number of A[i]
        auto numOfWays = fibList[A[i] + 1];

        // we need to apply modulo to down-cast to int
        result[i] = numOfWays % (int)std::pow(2, B[i]);
    }

    return result;
}

/*
There are N empty glasses with a capacity of 1, 2, ..., N liters (there is exactly one glass of each unique capacity).
You want to pour exactly K liters of water into glasses.
 Each glass may be either full or empty (a glass cannot be partially filled).
 What is the minimum number of glasses that you need to contain K liters of water?
 Write a function: class Solution { public int solution(int N, int K); } that, given two integers N and K, returns the minimum number of glasses that are needed to contain exactly K liters of water.
 If it is not possible to pour exactly K liters of water into glasses then the function should return -1.

 Examples: 1. Given N = 5 and K = 8, the function should return 2. There are five glasses of capacity 1, 2, 3, 4 and 5. You can use two glasses with capacity 3 and 5 to hold 8 liters of water.
 2. Given N = 4 and K = 10, the function should return 4. You must use all the glasses to contain 10 liters of water.
 3. Given N = 1 and K = 2, the function should return -1. There is only one glass with capacity 1, so you cannot pour 2 liters of water.
 4. Given N = 10 and K = 5, the function should return 1. You can use the glass with capacity
 5. Write an efficient algorithm for the following assumptions:
 N is an integer within the range [1..1,000,000]; K is an integer within the range [1..1,000,000,000].
*/

int pourWaterIntoGlasses(int N, int K){
	// N: number of glasses
	// K: liters of water

    if(N >= K) {
		return 1;
	}
    
	// 1 + 2 + 3 + .... + N = (N*(N+1))/2
	long sum = ((long)N*(N+1))/2;
	
    if(sum < K){ // not possible to pour all the water!
        return -1;
    } else if(sum == K){ // it takes all the glasses to pour all the water 
        return N;
    }

	int leftOver = K;
	int biggestGlass = N;
	int usedGlasses = 0;
	while(leftOver) {
		++usedGlasses;
		if(leftOver > biggestGlass) {
			leftOver -= biggestGlass;
			--biggestGlass;
		} else {
			break;
		}
	}
	
	return usedGlasses;

    //Alternative solution
/*
    long longN = N;
	long longK = K;
	long max = longN * (longN + 1) / 2;
	int numGlasses = -1;
	if (longK <= max)
	{
		numGlasses = (int) (longN - (int) ((1 + std::sqrt(1 + 4 * (-2 * longK + longN * longN + longN))) / 2) + 1);
	} // if

	return numGlasses;
*/
}

#pragma region Operations

template<typename T>
void printArray(const std::vector<T> & tList) {
    std::cout << "Printing List: ";
    for(auto & iter : tList) std::cout << iter << " ";
    std::cout << std::endl;
}

void binaryGap_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Binary Gap operations..." << std::endl;
    
    std::cout << "6 -> " << binaryGap(6) << std::endl;
    std::cout << "328 -> " << binaryGap(328) << std::endl;
    std::cout << "1162 -> " << binaryGap(1162) << std::endl;
    std::cout << "51712 -> " << binaryGap(51712) << std::endl;
    std::cout << "1610612737 -> " << binaryGap(1610612737) << std::endl;

    std::cout << "Ending Binary Gap operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void rotatingArray_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Rotating Array operations..." << std::endl;
    
    std::vector<int> arry = {3, 8, 9, 7, 6};

    std::cout << "BEFORE" << std::endl;
    printArray<int>(arry);

    rotatingArray(arry, 3);

    std::cout << "AFTER" << std::endl;
    printArray<int>(arry);

    std::cout << "Ending Rotating Array operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void oddOccurrenceInArray_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Odd Occurrence in Array operations..." << std::endl;

    std::vector<int> arry = {9, 3, 9, 3, 9, 7, 9};
    printArray<int>(arry);
    std::cout << "Odd Occurrence in Array: " << oddOccurrenceInArray(arry) << std::endl;

    std::cout << "Ending Odd Occurrence in Array operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void missingElement_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Missing Element operations..." << std::endl;

    std::vector<int> arry = {2, 3, 1, 5}; // 4 is missing!
    printArray<int>(arry);
    std::cout << "Missing Element in Array: " << missingElement(arry) << std::endl;

    std::cout << "Ending Missing Element operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void tapeEquilibrium_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Tape Equilibrium operations..." << std::endl;

    std::vector<int> arry = {3, 1, 2, 4, 3};
    printArray<int>(arry);
    std::cout << "Minimum Difference of Parts in Array: " << tapeEquilibrium(arry) << std::endl;

    std::cout << "Ending Tape Equilibrium operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void frogRiverOne_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Frog River One operations..." << std::endl;

    std::vector<int> arry = {1, 3, 1, 4, 2, 3, 5, 4};
    int X = 5;
    printArray<int>(arry);
    std::cout << "Position when path established: " << frogRiverOne(X, arry) << std::endl;

    std::cout << "Ending Frog River One operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void maxCounters_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Max Counters operations..." << std::endl;

    std::vector<int> arry = {3, 4, 4, 6, 1, 4, 4};
    int X = 5;
    printArray<int>(arry);
    auto result = maxCounters(X, arry);
    std::cout << "After applying MaxCounters:" << std::endl;
    printArray(result);

    std::cout << "Ending Max Counters operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void missingInteger_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Missing Integer operations..." << std::endl;

    std::vector<int> arry = {1, 3, 6, 4, 1, 2};
    printArray<int>(arry);
    std::cout << "Missing Integer from the Array: " << missingInteger(arry) << std::endl;

    arry = {-1, -3};
    printArray<int>(arry);
    std::cout << "Missing Integer from the Array: " << missingInteger(arry) << std::endl;

    arry = {1, 2, 3};
    printArray<int>(arry);
    std::cout << "Missing Integer from the Array: " << missingInteger(arry) << std::endl;

    std::cout << "Ending Missing Integer operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void permutationCheck_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Permutation Check operations..." << std::endl;

    std::vector<int> arry = {4, 1, 3, 2};
    printArray<int>(arry);
    std::cout << "Permutation Check Result: " << permutationCheck(arry) << std::endl;

    arry = {4, 1, 3};
    printArray<int>(arry);
    std::cout << "Permutation Check Result: " << permutationCheck(arry) << std::endl;

    std::cout << "Ending Permutation Check operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void countDiv_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Count Div operations..." << std::endl;

    int A = 6;
    int B = 11;
    int K = 2;
    std::cout << "Count Div for (" << "A: " << A << " B: " << B << " K: " << K << ") Result: " << countDiv(A, B, K) << std::endl;

    A = 0;
    B = 2000000000;
    K = 1;
    std::cout << "Count Div for (" << "A: " << A << " B: " << B << " K: " << K << ") Result: " << countDiv(A, B, K) << std::endl;

    std::cout << "Ending Count Div operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void genomicRangeQuery_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Genomic Range Query operations..." << std::endl;

    std::string S("CAGCCTA");
    std::vector<int> P = {2, 5, 0};
    std::vector<int> Q = {4, 5, 6};

    auto result = genomicRangeQuery(S, P, Q);
    printArray(result);

    std::cout << "Ending Genomic Range Query operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void minAvgTwoSlice_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Min Average Two Slice operations..." << std::endl;

    std::vector<int> arry = {4, 2, 2, 5, 1, 5, 8};

    printArray(arry);
    std::cout << "Starting index of minimal average slice: " << minAvgTwoSlice(arry) << std::endl;

    std::cout << "Ending Min Average Two Slice operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void passingCars_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Passing Cars operations..." << std::endl;

    std::vector<int> arry = {0, 1, 0, 1, 1};

    printArray(arry);
    std::cout << "Number of passing cars: " << passingCars(arry) << std::endl;

    std::cout << "Ending Passing Cars operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void distinct_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Distinct operations..." << std::endl;

    std::vector<int> arry = {2, 1, 1, 2, 3, 1};

    printArray(arry);
    std::cout << "Number of distinct values: " << distinct(arry) << std::endl;

    std::cout << "Ending Distinct operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void maxProductOfThree_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Max Product of Three operations..." << std::endl;

    std::vector<int> arry = {-3, 1, 2, -2, 5, 6};
    printArray(arry);
    std::cout << "Max product of three: " << maxProductOfThree(arry) << std::endl;

    arry = {1, 5, 3, 4, -6, -5, 6};
    printArray(arry);
    std::cout << "Max product of three: " << maxProductOfThree(arry) << std::endl;

    std::cout << "Ending Max Product of Three operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void numberOfDiscIntersections_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Number of Disc Intersections operations..." << std::endl;

    std::vector<int> arry = {1, 5, 2, 1, 4, 0};
    printArray(arry);
    std::cout << "Number of Disc Intersections " << numberOfDiscIntersections(arry) << std::endl;

    std::cout << "Ending Number of Disc Intersections operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void triangle_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Triangle operations..." << std::endl;

    std::vector<int> arry = {10, 2, 5, 1, 8, 20};
    printArray(arry);
    std::cout << "Triangle exists: " << (triangle(arry)==1 ? "Yes" : "No") << std::endl;

    arry = {10, 50, 5, 1};
    printArray(arry);
    std::cout << "Triangle exists: " << (triangle(arry)==1 ? "Yes" : "No") << std::endl;

    std::cout << "Ending Triangle operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void brackets_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Brackets operations..." << std::endl;

    std::string S("{[()()]}");
    std::cout << "String: " << S << std::endl;
    std::cout << "Properly Nested: " << (brackets(S)==1 ? "Yes" : "No") << std::endl;

    S = "([)()]";
    std::cout << "String: " << S << std::endl;
    std::cout << "Properly Nested: " << (brackets(S)==1 ? "Yes" : "No") << std::endl;

    std::cout << "Ending Brackets operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void fish_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Fish operations..." << std::endl;

    std::vector<int> size = {4, 3, 2, 1, 5};
    std::vector<int> direction = {0, 1, 0, 0, 0};
    std::cout << "(Fish Size)";
    printArray(size);
    std::cout << "(Fish Direction)";
    printArray(direction);
    std::cout << "Living fish count: " << fish(size, direction) << std::endl;

    std::cout << "Ending Fish operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void nested_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Nested operations..." << std::endl;

    std::string S("(()(())())");
    std::cout << "String: " << S << std::endl;
    std::cout << "Properly Nested: " << (nested(S)==1 ? "Yes" : "No") << std::endl;

    S = "())";
    std::cout << "String: " << S << std::endl;
    std::cout << "Properly Nested: " << (nested(S)==1 ? "Yes" : "No") << std::endl;

    std::cout << "Ending Nested operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void stoneWall_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Stone Wall operations..." << std::endl;

    std::vector<int> heights = {8, 8, 5, 7, 9, 8, 7, 4, 8};
    std::cout << "(Wall Heights)";
    printArray(heights);
    std::cout << "Minimum Wall Count: " << stoneWall(heights) << std::endl;

    heights = {1, 2, 3, 3, 2, 1};
    std::cout << "(Wall Heights)";
    printArray(heights);
    std::cout << "Minimum Wall Count: " << stoneWall(heights) << std::endl;

    heights = {2, 5, 1, 4, 6, 7, 9, 10, 1};
    std::cout << "(Wall Heights)";
    printArray(heights);
    std::cout << "Minimum Wall Count: " << stoneWall(heights) << std::endl;

    heights = {1, 1000000000, 1};
    std::cout << "(Wall Heights)";
    printArray(heights);
    std::cout << "Minimum Wall Count: " << stoneWall(heights) << std::endl;

    std::cout << "Ending Stone Wall operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void dominator_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Dominator operations..." << std::endl;

    std::vector<int> arry = {3, 4, 3, 2, 3, -1, 3, 3};
    printArray(arry);
    std::cout << "Dominator Index: " << dominator(arry) << std::endl;

    arry = {7, 7, 7, 7, 7, 7, 3, 2, 5};
    printArray(arry);
    std::cout << "Dominator Index: " << dominator(arry) << std::endl;

    std::cout << "Ending Dominator operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void equiLeader_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting EquiLeader operations..." << std::endl;

    std::vector<int> arry = {4, 3, 4, 4, 4, 2};
    printArray(arry);
    std::cout << "EquiLeader count: " << equiLeader(arry) << std::endl;

    std::cout << "Ending EquiLeader operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void maxDoubleSliceSum_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Max Double Slice Sum operations..." << std::endl;

    std::vector<int> arry = {3, 2, 6, -1, 4, 5, -1, 2};
    printArray(arry);
    std::cout << "Max Double Slice Sum: " << maxDoubleSliceSum(arry) << std::endl;

    std::cout << "Ending Max Double Slice Sum operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void maxProfit_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Max Profit operations..." << std::endl;

    std::vector<int> arry = {23'171, 21'011, 21'123, 21'366, 21'013, 21'367};
    printArray(arry);
    std::cout << "Max Profit: " << maxProfit(arry) << std::endl;

    std::cout << "Ending Max Profit operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void maxSliceSum_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Max Slice Sum operations..." << std::endl;

    std::vector<int> arry = {3, 2, -6, 4, 0};
    printArray(arry);
    std::cout << "Max Profit: " << maxSliceSum(arry) << std::endl;

    arry = {-10};
    printArray(arry);
    std::cout << "Max Profit: " << maxSliceSum(arry) << std::endl;

    arry = {-2, 1};
    printArray(arry);
    std::cout << "Max Profit: " << maxSliceSum(arry) << std::endl;

    std::cout << "Ending Max Slice Sum operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void countFactors_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Count Factors operations..." << std::endl;

    auto N = 24;
    std::cout << "Number of Factors for " << N << " : " << countFactors(N) << std::endl;

    N = 2'147'483'647;
    std::cout << "Number of Factors for " << N << " : " << countFactors(N) << std::endl;

    std::cout << "Ending Count Factors operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void flags_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Flags operations..." << std::endl;

    std::vector<int> arry = {1, 5, 3, 4, 3, 4, 1, 2, 3, 4, 6, 2};
    printArray(arry);
    std::cout << "Max Number of Flags: " << flags(arry) << std::endl;

    arry = {3, 2, 1};
    printArray(arry);
    std::cout << "Max Number of Flags: " << flags(arry) << std::endl;

    arry = {1, 3, 2};
    printArray(arry);
    std::cout << "Max Number of Flags: " << flags(arry) << std::endl;

    std::cout << "Ending Flags operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void minPerimeterRectangle_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Min Perimeter Rectangle operations..." << std::endl;

    auto N = 30;
    std::cout << "Min Perimeter for " << N << " : " << minPerimeterRectangle(N) << std::endl;

    N = 1;
    std::cout << "Min Perimeter for " << N << " : " << minPerimeterRectangle(N) << std::endl;

    N = 1'000'000'000;
    std::cout << "Min Perimeter for " << N << " : " << minPerimeterRectangle(N) << std::endl;

    std::cout << "Ending Min Perimeter Rectangle operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void peaks_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Peaks operations..." << std::endl;

    std::vector<int> arry = {1, 2, 3, 4, 3, 4, 1, 2, 3, 4, 6, 2};
    printArray(arry);
    std::cout << "Max Number of Blocks: " << peaks(arry) << std::endl;

    std::cout << "Ending Peaks operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void countNonDivisible_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Count Non-Divisible operations..." << std::endl;

    std::vector<int> arry = {3, 1, 2, 3, 6};
    printArray(arry);
    auto result = countNonDivisible(arry);
    std::cout << "Amount of Non-Divisors" << std::endl;
    printArray(result);

    std::cout << "Ending Count Non-Divisible operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void countSemiprimes_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Count Semiprimes operations..." << std::endl;

    auto N = 26;
    std::vector<int> arryP = {1, 4, 16};
    std::vector<int> arryQ = {26, 10, 20};
    std::cout << "Parameters => N: " << N << std::endl << "P: ";
    printArray(arryP);
    std::cout << "Q: ";
    printArray(arryQ);

    auto result = countSemiprimes(N, arryP, arryQ);
    std::cout << "Amount of Semiprimes" << std::endl;
    printArray(result);

    std::cout << "Ending Count Semiprimes operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void chocolatesByNumbers_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Chocolates by Numbers operations..." << std::endl;

    //N = 10 & M = 4 => 5 (You will eat the following chocolates: 0, 4, 8, 2, 6)
    auto N = 10;
    auto M = 4;
    auto numOfChocolates = chocolatesByNumbers(10, 4);
    
    std::cout << "Parameters => N: " << N << std::endl << "M: " << M << std::endl;
    std::cout << "Amount of chocolates eaten: " << numOfChocolates << std::endl;

    std::cout << "Ending Chocolates by Numbers operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void commonPrimeDivisors_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Common Prime Divisors operations..." << std::endl;

    // ([15, 10, 3], [75,30, 5]) -> 1
    std::vector<int> A = {15, 10, 3};
    std::vector<int> B = {75, 30, 5};
    std::cout << "A: ";
    printArray(A);
    std::cout << "B: ";
    printArray(B);

    auto result = commonPrimeDivisors(A, B);
    std::cout << "Number of Pairs having the same prime divisors: " << result << std::endl;

    // ([1], [1]) -> 1
    A = {1};
    B = {1};
    std::cout << "A: ";
    printArray(A);
    std::cout << "B: ";
    printArray(B);

    result = commonPrimeDivisors(A, B);
    std::cout << "Number of Pairs having the same prime divisors: " << result << std::endl;

    // ([2, 1, 2], [1, 2, 2]) -> 1
    A = {2, 1, 2};
    B = {1, 2, 2};
    std::cout << "A: ";
    printArray(A);
    std::cout << "B: ";
    printArray(B);

    result = commonPrimeDivisors(A, B);
    std::cout << "Number of Pairs having the same prime divisors: " << result << std::endl;

    std::cout << "Ending Common Prime Divisors operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void fibFrog_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting FibFrog operations..." << std::endl;

    // {0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0} => 3 (F[5]=5, F(3)=2, F(5)=5)
    std::vector<int> A = {0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0};
    std::cout << "A: ";
    printArray(A);
    
    auto result = fib_frog(A);
    std::cout << "Minimum number of jumps required: " << result << std::endl;

    // {0, 0, 0, 1, 0} => -1
    A = {0, 0, 0, 1, 0};
    std::cout << "A: ";
    printArray(A);
    
    result = fib_frog(A);
    std::cout << "Minimum number of jumps required: " << result << std::endl;

    // {} => 1
    A = {};
    std::cout << "A: ";
    printArray(A);
    
    result = fib_frog(A);
    std::cout << "Minimum number of jumps required: " << result << std::endl;

    std::cout << "Ending FibFrog operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void ladder_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Ladders operations..." << std::endl;

    // A = {4, 4, 5, 5, 1} & B = {3, 2, 4, 3, 1} => {5, 1, 8, 0, 1}
    std::vector<int> A = {4, 4, 5, 5, 1};
    std::vector<int> B = {3, 2, 4, 3, 1};
    std::cout << "A: ";
    printArray(A);
    std::cout << "B: ";
    printArray(B);

    auto result = ladder(A, B);
    std::cout << "Number of different ways to climb to the top of the ladder: " << std::endl;
    printArray(result);

    std::cout << "Ending Ladders operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void pourWaterIntoGlasses_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Pour Water into Glasses operations..." << std::endl;

    auto N = 5;
    auto K = 8;
    auto result = pourWaterIntoGlasses(N, K); // 2
    std::cout << "Number of glasses to pour all water: " << result << " for number of glasses: " << N << " & liters of water: " << K << std::endl;

    N = 4;
    K = 10;
    result = pourWaterIntoGlasses(N, K); // 4
    std::cout << "Number of glasses to pour all water: " << result << " for number of glasses: " << N << " & liters of water: " << K << std::endl;

    N = 1;
    K = 2;
    result = pourWaterIntoGlasses(N, K); // -1
    std::cout << "Number of glasses to pour all water: " << result << " for number of glasses: " << N << " & liters of water: " << K << std::endl;

    N = 10;
    K = 5;
    result = pourWaterIntoGlasses(N, K); // 1
    std::cout << "Number of glasses to pour all water: " << result << " for number of glasses: " << N << " & liters of water: " << K << std::endl;

    N = 1'000'000;
    K = 1'000'000'000;
    result = pourWaterIntoGlasses(N, K); // 1001
    std::cout << "Number of glasses to pour all water: " << result << " for number of glasses: " << N << " & liters of water: " << K << std::endl;

    std::cout << "Ending Pour Water into Glasses operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

#pragma endregion Operations

void codility_operations() {
    binaryGap_operations();

    rotatingArray_operations();

    oddOccurrenceInArray_operations();

    missingElement_operations();

    tapeEquilibrium_operations();

    frogRiverOne_operations();

    maxCounters_operations();

    missingInteger_operations();

    permutationCheck_operations();

    countDiv_operations();

    genomicRangeQuery_operations();

    minAvgTwoSlice_operations();

    passingCars_operations();

    distinct_operations();

    maxProductOfThree_operations();

    numberOfDiscIntersections_operations();

    triangle_operations();

    brackets_operations();

    fish_operations();

    nested_operations();

    stoneWall_operations();

    dominator_operations();

    equiLeader_operations();

    maxDoubleSliceSum_operations();

    maxProfit_operations();

    maxSliceSum_operations();

    countFactors_operations();

    flags_operations();

    minPerimeterRectangle_operations();

    peaks_operations();

    countNonDivisible_operations();

    countSemiprimes_operations();

    chocolatesByNumbers_operations();

    commonPrimeDivisors_operations();

    fibFrog_operations();

    ladder_operations();

    pourWaterIntoGlasses_operations();
}

}

#endif // _CODILITY_TASKS_H_
