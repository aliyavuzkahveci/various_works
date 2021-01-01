
#ifndef _FIBONACCI_IMPL_H_
#define _FIBONACCI_IMPL_H_

#include <cmath>
#include <iostream>

namespace fibonacci {

// fibonacci numbers: 
// 0 1 2 3 4 5 6 7  8  9 ...
// 0 1 1 2 3 5 8 13 21 34 ...

/* Method-1 START */
// f(n) = f(n-1) + f(n-2)
int fibonacci(int n) {
    if(n == 0 || n == 1) {
        return n;
    } else {
        return fibonacci(n-1) + fibonacci(n-2);
    }
}
/* Method-1 END*/


/* Method-2 START*/
// dynamic programming
int fibonacci_2(int n) {
    int a = 0, b = 1, c;
    for(int i=2; i<=n; ++i) {
        c = a + b;
        a = b;
        b = c;
    }

    return b;
}
/* Method-2 END*/

/* Method-3 START*/
// matrix multiplication 
/*
(1 1        f(n+1) f(n)
 1 0) ^ n = f(n)   f(n-1)
*/
// multiply matrix and writes the result back into F
void multiply(int F[2][2], int M[2][2]) {
    auto x_00 = F[0][0] * M[0][0] +
                F[0][1] * M[1][0];
    
    auto x_01 = F[0][0] * M[0][1] +
                F[0][1] * M[1][1];
    
    auto x_10 = F[1][0] * M[0][0] +
                F[1][1] * M[1][0];
    
    auto x_11 = F[1][0] * M[0][1] +
                F[1][1] * M[1][1];
    
    F[0][0] = x_00;
    F[0][1] = x_01;
    F[1][0] = x_10;
    F[1][1] = x_11;
}

// raise to the number n and puts the result back into F
void power(int F[2][2], int n) {
    int M[2][2] = {{1, 1}, {1, 0}};

    for(size_t i=2; i<=n; ++i) {
        multiply(F, M);
    }
}
// Time complexity: O(n)
// Space complexity: O(1)
int fibonacci_3(int n) {
    if(n == 0) {
        return n;
    }

    int F[2][2] = {{1, 1}, {1, 0}};
    power(F, n-1);

    return F[0][0];
}
/* Method-3 END*/

/* Method-3_Optimized START*/
void power_optimized(int F[2][2], int n) {
    if(n == 0 || n == 1) {
        return;
    }

    int M[2][2] = {{1, 1}, {1, 0}};

    power(F, n/2);
    multiply(F, F);
    if(n%2 == 1) {
        multiply(F, M);
    }
}
// Time complexity: O(logn)
// Space complexity: O(1)
int fibonacci_3_optimized(int n) {
    if(n == 0) {
        return n;
    }

    int F[2][2] = {{1, 1}, {1, 0}};
    power_optimized(F, n-1);

    return F[0][0];
}
/* Method-3_Optimized END*/

/* Method-4 START*/
// f(n) = (((1 + sqrt(5)) / 2)^n) / sqrt(5)

// Time complexity: O(1)
// Space complexity: O(1)
int fibonacci_4(int n) {
    if(n == 0 || n == 1) {
        return n;
    }

    double phi = (1 + std::sqrt(5)) / 2;
    int fib_n = std::round(std::pow(phi, n) / std::sqrt(5));

    return fib_n;
}
/* Method-4 END*/

void fibonacci_trials() {
    auto n = 9;
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Method1: fibonacci(9) = " << fibonacci(9) << std::endl;
    std::cout << "Method2: fibonacci_2(9) = " << fibonacci_2(9) << std::endl;
    std::cout << "Method3: fibonacci_3(9) = " << fibonacci_3(9) << std::endl;
    std::cout << "Method3: fibonacci_3_optimized(9) = " << fibonacci_3_optimized(9) << std::endl;
    std::cout << "Method4: fibonacci_4(9) = " << fibonacci_4(9) << std::endl;
    std::cout << "=======================================" << std::endl;
}

} // fibonacci

#endif // _FIBONACCI_IMPL_H_
