
#ifndef _ONLINE_ASSESSMENTS_H_
#define _ONLINE_ASSESSMENTS_H_

#include <stack>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace online_assessments
{

bool isBalanced(std::string &checkStr)
{
    std::stack<char> cStack;

    for (const auto &ch : checkStr)
    {
        if (ch == '(' || ch == '[')
        {
            cStack.push(ch);
        }
        else if (cStack.empty())
        {
            return false;
        }
        else if (ch == ')' && cStack.top() != '(')
        {
            return false;
        }
        else if (ch == ']' && cStack.top() != '[')
        {
            return false;
        }
        else
        {
            cStack.pop();
        }
    }

    return cStack.empty();
}

bool isPalindrome(const std::string &word)
{
    if (word.length() < 2) // covers empty and 1-char strings! 
        return true;

    int endPos = word.length() - 1;
    for (int i = 0; i < (word.length() / 2); ++i)
    {
        if (std::isalpha(word[i]) && std::isalpha(word[endPos - i]))
        {
            if ((word[i] == word[endPos - i]) || std::abs(word[i] - word[endPos - i]) == 32)
                continue;
            else
                return false;
        }
        else
            return false;
    }
    return true;
}

size_t eightBallProblem(const std::vector<size_t>& ballList) {
    // assume the size of list is 8!!!
    // First group the balls into 3 as 
    // 1st group: index 0, 1, 2
    // 2nd group: index 3, 4, 5
    // 3rd group: index 6, 7
    auto group1_weight = ballList[0] + ballList[1] + ballList[2];
    auto group2_weight = ballList[3] + ballList[4] + ballList[5];
    auto group3_weight = ballList[6] + ballList[7];
    //Step1: compare the weights of group 1 and group 2
    if(group1_weight == group2_weight) {
        // defective ball is in group 3 (index 6 or 7)
        if(ballList[6] < ballList[7]) {
            return 6;
        } else {
            return 7;
        }
    } else if(group1_weight < group2_weight) {
        // defective ball is in group 1
        if(ballList[0] < ballList[1]) {
            return 0;
        } else if(ballList[0] > ballList[1]) {
            return 1;
        } else {
            return 2;
        }
    } else { // if(group1_weight > group2_weight) {
        // defective ball is in group 2
        if(ballList[3] < ballList[4]) {
            return 3;
        } else if(ballList[3] > ballList[4]) {
            return 4;
        } else {
            return 5;
        }
    }
}

int countAllPalindromeSubstrings(std::string& str) {
    // counts and prints substrings with length more than 1
    // excludes single chars :)
    std::cout << "****************************************" << std::endl;
    auto size = str.size();
    if(size == 1) {
        std::cout << "Palindrome Sub-Strings:" << std::endl << str;
        return 1;
    }

    // can be solved by dynamic programming
    // initialize palindrome table
    std::vector<std::string> paliSubstrings;
    
    std::vector<std::vector<bool>> paliList(size, std::vector<bool>(size, false));
    std::vector<std::vector<int>> paliCountList(size, std::vector<int>(size, 0));

    // palindrome of 1-char sub-strings
    for(size_t i=0; i<size; ++i) {
        paliList[i][i] = true;
    }

    // palindrome of substrings with length 2
    for(size_t i=0; i<size-1; ++i) {
        if(str[i] == str[i+1]) {
            paliList[i][i+1] = true;
            paliCountList[i][i+1] = 1;
            paliSubstrings.push_back(str.substr(i, 2));
        }
    }

    // palindrome check for sub-strings with length more than 2
    // this looping is similar to Matrix Chain Multiplication
    // we will start the gap from 2 and increase the sub-string from both ends 1 by 1
    for(size_t gap=2; gap<size; ++gap) {
        for(size_t i=0; i<size-gap; ++i) {
            // i: start index & j: end index
            auto j = i + gap;

            // if current string is palindrome
            if(str[i] == str[j] && paliList[i+1][j-1]) {
                paliList[i][j] = true;
            }

            // simple formula: result = 1 + (x + y) - z
            // 1: add only if the current substring is palindrome!!!
            // x: previous rest palindrome substring count => dp[i][j-1]
            // y: after rest palindrome substring count => dp[i+1][j]
            // z: common palindrome substring count => dp[i+1][j-1]
            paliCountList[i][j] = paliCountList[i][j-1] + paliCountList[i+1][j] - paliCountList[i+1][j-1];
            if(paliList[i][j]) {
                paliSubstrings.push_back(str.substr(i, (j-i+1))); // since i & j are indices, we need to add 1 to include jth position!
                ++paliCountList[i][j]; // add 1 as itself!
            }
        }
    }

    std::cout << "Palindrome Sub-String Count for \"" << str << "\" is " << paliCountList[0][size-1] << std::endl;
    std::cout << "Palindrome Sub-Strings:" << std::endl;
    for(auto const & iter : paliSubstrings) {
        std::cout << iter << std::endl;
    }

    return paliCountList[0][size-1];
}

} // namespace online_assessments

#endif // _ONLINE_ASSESSMENTS_H_
