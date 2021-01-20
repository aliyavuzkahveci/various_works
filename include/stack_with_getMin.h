
#ifndef _STACK_WITH_GETMIN_H_
#define _STACK_WITH_GETMIN_H_

#include <stack>
#include <iostream>

namespace custom_stack {

// Custom Stack implementation with getMin() in O(1) complexity
// Time complexity: O(1)
// Space complexity: O(1)
template <typename T>
class CustomStack {

public:
    explicit CustomStack(int capacity = 100) : m_capacity(capacity), m_currentSize(0), m_min(T()) {
    }

    virtual ~CustomStack() = default;

    T getMin() {
        return m_min;
    }

    T top() {
        if(m_stack.empty()) {
            std::cout << "stack is empty! Returning NULL!" << std::endl;
            return T();
        }

        auto top = m_stack.top();
        return (top < m_min) ? m_min : top; // top element could the minimum! 
    }

    void push(const T& element) {
        if(m_currentSize == m_capacity) {
            std::cout << "Stack is full! No new element can be added!" << std::endl;
            return;
        } else if(m_stack.empty()) { // adding the first element
            m_min = element;
            m_stack.push(element);
        } else if(element < m_min) { // need to rearrange the min
            m_stack.push(2*element - m_min); // store the previous min element info to return when the new min is removed! 
            m_min = element;
        } else {
            m_stack.push(element);
        }
        ++m_currentSize;
    
    }

    void pop() {
        if(m_stack.empty()) {
            std::cout << "stack is empty! No removal can be performed!" << std::endl;
            return;
        }

        auto top = m_stack.top();
        m_stack.pop();
        --m_currentSize;

        if(top < m_min) { // removed element was the minimum
            m_min = 2*m_min - top; // restore the previous min element as the current minimum!
        }

        if(m_stack.empty()) {
            m_min = T();
        }
    }

private:
    std::stack<T> m_stack;
    T m_min;
    int m_capacity;
    int m_currentSize;

};

void custom_stack_trials() {
    CustomStack<int> customStack;
    customStack.push(3);
    customStack.push(2);
    customStack.push(5);
    customStack.push(9);
    customStack.push(1);
    customStack.push(21);
    customStack.push(0);

    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "CustomStack: Min Element: " << customStack.getMin() << " Top Element: " << customStack.top() << std::endl;
    customStack.pop();
    std::cout << "CustomStack: Min Element: " << customStack.getMin() << " Top Element: " << customStack.top() << std::endl;
    customStack.pop();
    customStack.pop();
    customStack.pop();
    customStack.pop();
    std::cout << "CustomStack: Min Element: " << customStack.getMin() << " Top Element: " << customStack.top() << std::endl;
    std::cout << "=======================================" << std::endl;
}

}

#endif // _STACK_WITH_GETMIN_H_