
#ifndef _CODILITY_TASKS_H_
#define _CODILITY_TASKS_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <stack>

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

}

#endif // _CODILITY_TASKS_H_
