
#include "bitwise_operations.h"
#include "online_assessments.h"
#include "codility_tasks.h"

#include <algorithm>
#include <math.h>

using namespace bitwise_operations;
using namespace online_assessments;
using namespace codility_tasks;

template<typename T>
void printArray(std::vector<T> & tList) {
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

int main()
{
    apply_bitwise_operations();

    isbalanced_operations();

    ispalindrome_operations();

    eightBallProblem_operations();

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

    return 0;
}