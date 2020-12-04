
#ifndef _ONLINE_ASSESSMENTS_H_
#define _ONLINE_ASSESSMENTS_H_

#include <stack>
#include <string>
#include <vector>

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
        if (word.length() == 0)
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

} // namespace online_assessments

#endif // _ONLINE_ASSESSMENTS_H_
