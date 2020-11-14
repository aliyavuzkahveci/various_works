
#ifndef _ONLINE_ASSESSMENTS_H_
#define _ONLINE_ASSESSMENTS_H_

#include <stack>
#include <string>

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

} // namespace online_assessments

#endif // _ONLINE_ASSESSMENTS_H_
