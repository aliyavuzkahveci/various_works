
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

}

#endif // _CODILITY_TASKS_H_
