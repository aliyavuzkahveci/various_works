#ifndef _BITWISE_OPERATIONS_
#define _BITWISE_OPERATIONS_

#include <iostream>
#include <bitset>
#include <math.h>

namespace bitwise_operations
{
unsigned int set_bit(unsigned int value, unsigned int position)
{
    unsigned int mask = 1 << position; // 00000001 => (position=2) 00000100
    unsigned int result = value | mask;

    std::cout << "value:    " << std::bitset<8>(value) << std::endl;
    std::cout << "position: " << std::bitset<8>(position) << std::endl;
    std::cout << "mask:     " << std::bitset<8>(mask) << std::endl;
    std::cout << "result:   " << std::bitset<8>(result) << std::endl;

    return result;
}

unsigned int clear_bit(unsigned int value, unsigned int position)
{
    unsigned int mask = 1 << position; // 00000001 => (position=2) 00000100
    unsigned int result = value & ~mask;

    std::cout << "value:    " << std::bitset<8>(value) << std::endl;
    std::cout << "position: " << std::bitset<8>(position) << std::endl;
    std::cout << "mask:     " << std::bitset<8>(mask) << std::endl;
    std::cout << "~mask:    " << std::bitset<8>(~mask) << std::endl;
    std::cout << "result:   " << std::bitset<8>(result) << std::endl;

    return result;
}

unsigned int flip_bit(unsigned int value, unsigned int position)
{
    unsigned mask = 1 << position; // 00000001 => (position=2) 00000100
    unsigned int result = value ^ mask;

    std::cout << "value:    " << std::bitset<8>(value) << std::endl;
    std::cout << "position: " << std::bitset<8>(position) << std::endl;
    std::cout << "mask:     " << std::bitset<8>(mask) << std::endl;
    std::cout << "result:   " << std::bitset<8>(result) << std::endl;

    return result;
}

bool is_bit_set(unsigned int value, unsigned int position)
{
    unsigned int shifted = value >> position;
    unsigned int result = 1 & shifted;

    std::cout << "value:    " << std::bitset<8>(value) << std::endl;
    std::cout << "position: " << std::bitset<8>(position) << std::endl;
    std::cout << "shifted:  " << std::bitset<8>(shifted) << std::endl;
    std::cout << "result:   " << std::bitset<8>(result) << std::endl;
    std::cout << "result:   " << std::boolalpha << bool(result) << std::endl;

    return result;
}

unsigned int modify_bit(unsigned int value, unsigned int position, int state)
{
    unsigned int mask = 1 << position;
    unsigned int result = (value & ~mask) | (-state & mask); // -state => 1's complement!!!

    std::cout << "value:     " << std::bitset<8>(value) << std::endl;
    std::cout << "position:  " << std::bitset<8>(position) << std::endl;
    std::cout << "state:     " << std::bitset<8>(state) << std::endl;
    std::cout << "mask:      " << std::bitset<8>(mask) << std::endl;
    std::cout << "~mask:     " << std::bitset<8>(~mask) << std::endl;
    std::cout << "-state:    " << std::bitset<8>(-state) << std::endl;
    std::cout << "result:    " << std::bitset<8>(result) << std::endl;

    return result;
}

bool check_if_even(unsigned int value)
{
    // 0110 & 0001 = 0000
    bool result = (value & 1) == 0;

    std::cout << "result:    " << std::boolalpha << result << std::endl;

    return result;
}

bool check_if_power_of_two(unsigned int value)
{
    bool result = (value & (value - 1)) == 0;

    std::cout << "result:    " << std::boolalpha << result << std::endl;

    return result;
}

unsigned int count_num_of_ones(char value) noexcept
{
    int numOfOnes = value;
    char mask1 = 0b01010101;
    char mask2 = 0b00110011;
    char mask3 = 0b00001111;

    char shifted = numOfOnes >> 1;
    numOfOnes = (numOfOnes & mask1) + (shifted & mask1);

    shifted = numOfOnes >> 2;
    numOfOnes = (numOfOnes & mask2) + (shifted & mask2);

    shifted = numOfOnes >> 4;
    numOfOnes = (numOfOnes & mask3) + (shifted & mask3);

    std::cout << "number of 1's in " << std::bitset<8>(value) << " is: " << numOfOnes << std::endl;

    return numOfOnes;
}

unsigned int count_different_bits(char value1, char value2) noexcept
{
    // ony different bits will provide result as 1 via XOR operation
    char output = value1 ^ value2;

    unsigned int diffBits = count_num_of_ones(output);

    std::cout << "value1          : " << std::bitset<8>(value1) << std::endl;
    std::cout << "value2          : " << std::bitset<8>(value2) << std::endl;
    std::cout << "value1 ^ value2 : " << std::bitset<8>(output) << std::endl;
    std::cout << "num of diff bit : " << diffBits << std::endl;

    return diffBits;
}

unsigned int roundUpToNextPowerOfTwo(unsigned long int value)
{
    unsigned int twoTo31 = 2147483648; // 2^31

    if (value == 0 || value > twoTo31)
    {
        throw std::out_of_range("Given value is outside of the allowed range!");
    }

    std::cout << "--BEFORE OPERATION--" << std::endl;
    std::cout << "value=" << value << std::endl;

    value--;
    value |= value >> 1;  // handle 2 bit numbers
    value |= value >> 2;  // handle 4 bit numbers
    value |= value >> 4;  // handle 8 bit numbers
    value |= value >> 8;  // handle 16 bit numbers
    value |= value >> 16; // handle 32 bit numbers
    value++;

    std::cout << "--AFTER OPERATION--" << std::endl;
    std::cout << "value=" << value << std::endl;

    return value;
}

template <typename T>
void swapXOR(T &value1, T &value2) noexcept
{
    //std::cout << "--Before SWAP--" << std::endl;
    //std::cout << "value1: " << value1 << " value2: " << value2 << std::endl;

    value1 = value1 ^ value2;
    value2 = value2 ^ value1;
    value1 = value1 ^ value2;

    //std::cout << "--After SWAP--" << std::endl;
    //std::cout << "value1: " << value1 << " value2: " << value2 << std::endl;
}

template <typename T>
void swapWithTemp(T &value1, T &value2)
{
    //std::cout << "--Before SWAP--" << std::endl;
    //std::cout << "value1: " << value1 << " value2: " << value2 << std::endl;

    auto tempVal = value1;
    value1 = value2;
    value2 = tempVal;

    //std::cout << "--Before SWAP--" << std::endl;
    //std::cout << "value1: " << value1 << " value2: " << value2 << std::endl;
}

int absoluteValue(int value) noexcept
{
    const int bit31 = value >> 31; // valid for 32bit machines !!!
    int result = (value ^ bit31) - bit31;

    std::cout << "absolute value of " << value << " is " << result << std::endl;

    return result;
}

// adding two numbers without arithmetic operations.
int add(int x, int y) {
    // continue until no carry bit left
    while(y != 0) {
        // binary AND would provide common 1s in x and y
        int carry = x & y;

        // binary XOR would provide sum of x and y except the carry bits!
        x = x ^ y;

        // shifting carry by 1 gives the necessary sum
        y = carry << 1;
    }

    return x;
}

// subtracting a number from another one without arithmetic operations.
int subtract(int x, int y) {
    // continue until no borrow left 
    while(y != 0) {
        // binary AND of y with negated x would provide 1s in x with the same position of 0s in y
        int borrow = (~x) & y;

        // subtraction of bits in x and y where at least one of the bits is not set
        x = x ^ y;

        // shifting borrow by 1 gives the necessary results of subtracting y from x
        y = borrow << 1;
    }

    return x;
}

// comparison of two numbers without arithmetic operations.
bool compare(int x, int y) {
    // based on XOR operation
    // 0 XOR 1 = 1 XOR 0 = 1
    // 0 XOR 0 = 1 XOR 1 = 0
    return (x ^ y) == 0;
}

void print_binary(unsigned int n) {
    if(n > 1) {
        print_binary(n >> 1); // shift by 1 and recursively call itself!
    }

    std::cout << (n & 1);
}

#pragma region Operations

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

    std::cout << std::endl;
    std::cout << x << " in binary format: ";
    print_binary(x);
    std::cout << std::endl;

    std::cout << "Ending Bitwise operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}


#pragma endregion Operations


} // namespace bitwise_operations

#endif // _BITWISE_OPERATIONS_H_
