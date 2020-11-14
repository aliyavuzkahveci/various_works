
#include "bitwise_operations.h"
#include "online_assessments.h"

#include <algorithm>
#include <math.h>

using namespace bitwise_operations;
using namespace online_assessments;

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

int main()
{
    apply_bitwise_operations();

    isbalanced_operations();

    ispalindrome_operations();

    return 0;
}