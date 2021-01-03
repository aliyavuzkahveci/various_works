
#include "bitwise_operations.h"
#include "online_assessments.h"
#include "codility_tasks.h"
#include "codesignal_tasks.h"
#include "trial_codes.h"
#include "fibonacci_impl.h"
#include "knapsack_impl.h"

#include <math.h>

using namespace bitwise_operations;
using namespace online_assessments;
using namespace codility_tasks;
using namespace codesignal_tasks;
using namespace trial_codes;
using namespace fibonacci;
using namespace knapsack;

using namespace std::string_literals;

template<typename T>
void printArray(const std::vector<T> & tList) {
    std::cout << "Printing List: ";
    for(auto & iter : tList) std::cout << iter << " ";
    std::cout << std::endl;
}

void apply_bitwise_operations()
{
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Bitwise operations..." << std::endl;

    unsigned int result = set_bit(256, 5);
    std::cout << "set_bit result=" << result << std::endl
              << std::endl;

    result = clear_bit(255, 3);
    std::cout << "clear_bit result=" << result << std::endl
              << std::endl;

    result = flip_bit(1233, 7);
    std::cout << "flip_bit result=" << result << std::endl
              << std::endl;

    bool result_ = is_bit_set(233, 3);
    std::cout << "is_bit_set result=" << std::boolalpha << result_ << std::endl
              << std::endl;

    result = modify_bit(233, 3, 1);
    std::cout << "modify_bit result=" << result << std::endl
              << std::endl;

    result_ = check_if_even(233);
    std::cout << "check_if_even result=" << std::boolalpha << result_ << std::endl
              << std::endl;

    result_ = check_if_even(62);
    std::cout << "check_if_even result=" << std::boolalpha << result_ << std::endl
              << std::endl;

    result_ = check_if_power_of_two(128);
    std::cout << "check_if_even result=" << std::boolalpha << result_ << std::endl
              << std::endl;

    result_ = check_if_power_of_two(62);
    std::cout << "check_if_even result=" << std::boolalpha << result_ << std::endl
              << std::endl;

    result = count_num_of_ones(255);
    std::cout << "count_num_of_ones result=" << result << std::endl
              << std::endl;

    result = count_num_of_ones(127);
    std::cout << "count_num_of_ones result=" << result << std::endl
              << std::endl;

    result = count_num_of_ones(128);
    std::cout << "count_num_of_ones result=" << result << std::endl
              << std::endl;

    result = count_num_of_ones(0);
    std::cout << "count_num_of_ones result=" << result << std::endl
              << std::endl;

    result = count_different_bits(128, 127);
    std::cout << "count_different_bits result=" << result << std::endl
              << std::endl;

    result = count_different_bits(128, 128);
    std::cout << "count_different_bits result=" << result << std::endl
              << std::endl;

    result = count_different_bits(157, 33);
    std::cout << "count_different_bits result=" << result << std::endl
              << std::endl;

    try
    {
        unsigned long int value = 43951385;
        result = roundUpToNextPowerOfTwo(value);
        std::cout << "value=" << value << " is " << std::log2(value) << " power of 2" << std::endl;
        std::cout << "result=" << result << " is " << std::log2(result) << " power of 2" << std::endl;
        std::cout << std::endl;

        value = 5823502467;
        result = roundUpToNextPowerOfTwo(value);
        std::cout << "value=" << value << " is " << std::log2(value) << " power of 2" << std::endl;
        std::cout << "result=" << result << " is " << std::log2(result) << " power of 2" << std::endl;
        std::cout << std::endl;
    }
    catch (std::exception &ex)
    {
        std::cout << "Exception occurred: " << ex.what() << std::endl;
    }

    int value1 = 345;
    int value2 = 678;

    swapXOR<int>(value1, value2);
    std::cout << std::endl;
    swapWithTemp<int>(value1, value2);

    std::cout << std::endl;
    absoluteValue(-415263748);

    std::cout << std::endl;
    auto x = 98;
    auto y = 34;
    std::cout << "Sum of " << x << " and " << y << " is " << add(x, y) << std::endl;

    std::cout << std::endl;
    std::cout << "Subtract of " << x << " and " << y << " is " << subtract(x, y) << std::endl;

    std::cout << std::endl;
    std::cout << x << "==" << y << ": " << std::boolalpha << compare(x, y) << std::endl;
    x = 34;
    std::cout << x << "==" << y << ": " << std::boolalpha << compare(x, y) << std::endl;

    std::cout << "Ending Bitwise operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void isbalanced_operations()
{
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Is Balanced operations..." << std::endl;

    std::string name = "What is your name? ";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    name = "([])[]()";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    name = "((([([])]))())";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    name = "][";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    name = "([]]()";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    name = "";
    std::cout << "string: " << name << " => " << "isBalanced = " << std::boolalpha << isBalanced(name) << std::endl;

    std::cout << "Ending Is Balanced operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void ispalindrome_operations()
{
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Is Palindrome operations..." << std::endl;

    std::string name = "Deleveled";
    std::cout << "string: " << name << " => " << "isPalindrome = " << std::boolalpha << isPalindrome(name) << std::endl;

    std::cout << "Ending Is Palindrome operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void eightBallProblem_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting 8 Ball Problem operations..." << std::endl;
    std::vector<size_t> ballList = {5, 5, 5, 5, 5, 5, 5, 5};

    auto defectiveBallIndex = 1;
    ballList[defectiveBallIndex] = 4;
    auto foundIndex = eightBallProblem(ballList);
    std::cout << "8 Ball Problem: defective ball index:" << defectiveBallIndex << " vs found index: " << foundIndex << std::endl;
    ballList[defectiveBallIndex] = 5;

    defectiveBallIndex = 5;
    ballList[defectiveBallIndex] = 4;
    foundIndex = eightBallProblem(ballList);
    std::cout << "8 Ball Problem: defective ball index:" << defectiveBallIndex << " vs found index: " << foundIndex << std::endl;
    ballList[defectiveBallIndex] = 5;

    defectiveBallIndex = 7;
    ballList[defectiveBallIndex] = 4;
    foundIndex = eightBallProblem(ballList);
    std::cout << "8 Ball Problem: defective ball index:" << defectiveBallIndex << " vs found index: " << foundIndex << std::endl;
    ballList[defectiveBallIndex] = 5;

    std::cout << "Ending 8 Ball Problem operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void mutateTheArray_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Mutate The Array operations..." << std::endl;
    
    std::vector<int> list = {4, 0, 1, -2, 3};
    printArray(list);
    auto mutatedArray = mutateTheArray(5, list);
    std::cout << "Mutated Array" << std::endl;
    printArray(mutatedArray);

    std::cout << "Ending Mutate The Array operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void alternatingSort_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Alternating Sort operations..." << std::endl;
    
    std::vector<int> list = {1, 3, 5, 6, 4, 2};
    printArray(list);
    auto sorted = alternatingSort(list);
    std::cout << "Result: " << (sorted ? "array is sorted in strictly ascending order" : "array is NOT sorted in strictly ascending order") << std::endl;

    std::cout << "Ending Alternating Sort operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void countTinyPairs_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Count Tiny Pairs operations..." << std::endl;

    std::vector<int> a = {16, 1, 4, 2, 14};
    std::vector<int> b = {7, 11, 2, 0, 15};
    auto k=743;

    printArray(a);
    printArray(b);
    std::cout << "k: " << k << std::endl;
    
    auto numOfTinyPairs = countTinyPairs(a, b, k);
    std::cout << "Number of Tiny Pairs: " << numOfTinyPairs << std::endl;    

    std::cout << "Ending Count Tiny Pairs operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void meanGroups_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Mean Groups operations..." << std::endl;
    
    std::vector<std::vector<int>> a = {{3, 3, 4, 2}, {4, 4}, {4, 0, 3, 3}, {2, 3}, {3, 3, 3}};
    std::cout << "Input List:" << std::endl;
    for(auto const & iter : a) {
        printArray(iter);
    }
    auto meanGroupedList = meanGroups(a);
    std::cout << "Grouped List:" << std::endl;
    for(auto const & iter : meanGroupedList) {
        printArray(iter);
    }

    std::cout << "Ending Mean Groups operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void mergeStrings_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Merge Strings operations..." << std::endl;
    
    std::string s1("dce");
    std::string s2("cccbd");
    std::cout << "s1: " << s1 << " - s2: " << s2 << std::endl;
    auto merged = mergeStrings(s1, s2);
    std::cout << "Merged String: " << merged << std::endl;

    std::cout << "Ending Merge Strings operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void hashMap_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting HashMap operations..." << std::endl;
    
    std::vector<std::string> queryTypeList = {"insert"s, "insert"s, "addToValue"s, "addToKey"s, "get"s};
    std::vector<std::vector<int>> queryList = {{1, 2}, {2, 3}, {2}, {1}, {3}};
    printArray(queryTypeList);
    auto sumOfGets = hashMap(queryTypeList, queryList);
    std::cout << "Sum of All Gets: " << sumOfGets << std::endl;

    std::cout << "Ending HashMap operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void digitAnagrams_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Digit Anagrams operations..." << std::endl;
    
    std::vector<int> list = {345, 657, 134, 879, 987, 534, 435};
    printArray(list);
    auto numOfAnagrams = digitAnagrams(list);
    std::cout << "Number of digit anagrams:" << numOfAnagrams << std::endl;

    std::cout << "Ending Digit Anagrams operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

void concatenationsSum_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Concatenations Sum operations..." << std::endl;
    
    std::vector<int> a = {10, 2};
    printArray(a);
    auto concatSum = concatenationsSum(a);
    std::cout << "Sum of All Concatenations:" << concatSum << std::endl;

    a = {8};
    printArray(a);
    concatSum = concatenationsSum(a);
    std::cout << "Sum of All Concatenations:" << concatSum << std::endl;

    a = {1, 2, 3};
    printArray(a);
    concatSum = concatenationsSum(a);
    std::cout << "Sum of All Concatenations:" << concatSum << std::endl;

    std::cout << "Ending Concatenations Sum operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
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


void online_assessments_operations() {
    isbalanced_operations();

    ispalindrome_operations();

    eightBallProblem_operations();
}

void codesignal_operations() {
    mutateTheArray_operations();

    alternatingSort_operations();

    countTinyPairs_operations();

    mergeStrings_operations();

    meanGroups_operations();

    hashMap_operations();

    digitAnagrams_operations();

    concatenationsSum_operations();
}

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
}

int main()
{
    curly_braces_initializer();

    apply_bitwise_operations();

    online_assessments_operations();

    codesignal_operations();

    codility_operations();

    fibonacci_trials();

    knapsack_trials();

    return 0;
}