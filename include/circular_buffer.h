
#ifndef _CIRCULAR_BUFFER_H_
#define _CIRCULAR_BUFFER_H_

#include <memory>
#include <mutex>
#include <iostream>

namespace circular_buffer {

template <class T>
class CircularBuffer {

public:
    explicit CircularBuffer(int capacity) {
        m_buffer = std::make_unique<T[]>(capacity);
        m_capacity = capacity;
        m_headPos = 0;
        m_tailPos = 0;
        m_isFull = false;
    }

    // will overwrite the oldest value in the buffer!
    void put(const T& element) {
        std::lock_guard<std::mutex> lock(m_mutex);

        std::cout << "Adding " << element << std::endl;

        m_buffer[m_headPos] = element;

        // overwriting the oldest value in the buffer!
        if(m_isFull) {
            m_tailPos = (m_tailPos + 1) % m_capacity;
        }

        m_headPos = (m_headPos + 1) % m_capacity;

        // setting if the buffer is full
        m_isFull = m_headPos == m_tailPos;
    }

    T get() {
        if(empty()) {
            // returning teh default value of T type!
            return T();
        }

        auto tail = m_buffer[m_tailPos];

        m_tailPos = (m_tailPos + 1) % m_capacity;
        m_isFull = false;

        return tail;
    }

    virtual ~CircularBuffer() = default;

    size_t capacity() const {
        return m_capacity;
    }

    bool empty() const {
        return !m_isFull && m_headPos == m_tailPos;
    }

    bool full() const {
        return m_isFull;
    }

    size_t size() const {
        if(m_isFull) {
            return m_capacity;
        }

        if(m_headPos >= m_tailPos) {
            return (m_headPos - m_tailPos);
        } else {
            return (m_capacity + m_headPos - m_tailPos);
        }
    }

private:
    std::mutex m_mutex;
    std::unique_ptr<T[]> m_buffer;
    size_t m_capacity;
    size_t m_headPos;
    size_t m_tailPos;
    bool m_isFull;

    CircularBuffer() = delete; // deleting default constructor
    CircularBuffer(const CircularBuffer& rhs) = delete; // deleting copy constructor
    CircularBuffer& operator=(const CircularBuffer& rhs) = delete; // deleting assignment operator
    CircularBuffer(CircularBuffer&& rhs) = delete; // deleting move constructor
    CircularBuffer& operator=(CircularBuffer&& rhs) = delete; // deleting move assignment operator

};

#pragma region Operations

void circularBuffer_operations() {
    std::cout << std::endl << "=======================================" << std::endl;
    std::cout << "Starting Circular Buffer operations..." << std::endl;

    CircularBuffer<int> cb(10);
    std::cout << "CB capacity: " << cb.capacity() << std::endl;
    std::cout << "CB is empty?: " << cb.empty() << std::endl;
    std::cout << "CB is full?: " << cb.full() << std::endl;

    cb.put(1);
    std::cout << "CB is empty?: " << cb.empty() << std::endl; 
    cb.put(2);
    cb.put(3);
    cb.put(4);
    cb.put(5);
    cb.put(6);
    cb.put(7);
    cb.put(8);
    cb.put(9);
    cb.put(10);
    std::cout << "CB is full?: " << cb.full() << std::endl;
    cb.put(11);
    cb.put(12);

    std::cout << "Removing elements from the buffer..." << std::endl;
    while(cb.size()) {
        std::cout << "CB get " << cb.get() << std::endl;
    }
    std::cout << "CB is empty?: " << cb.empty() << std::endl;


    std::cout << "Ending Circular operations..." << std::endl;
    std::cout << "=======================================" << std::endl;
}

#pragma endregion Operations

}

#endif // _CIRCULAR_BUFFER_H_
