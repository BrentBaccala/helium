
/* based on a GitHub Gist from https://gist.github.com/thelinked/6997598 */

/* The original Gist implemented a locking queue (FIFO); this is a locking stack (LIFO).
 *
 * Brent Baccala added the set_num_workers() method.  If the stack is empty, and the
 * number of threads waiting equals the value passed to set_num_workers(), then
 * all the calls to waitAndPop return false.  For their normal behavior, they return true.
 *
 * The idea is to detect when the stack is empty and all the workers are waiting on it.
 */

#pragma once
#include <stack>
#include <mutex>
#include <condition_variable>

template<typename T>
class LockingStack
{
public:
    void set_num_workers(int num_workers)
    {
        total_workers = num_workers;
    }

    void push(T const& _data)
    {
        {
            std::lock_guard<std::mutex> lock(guard);
            stack.push(_data);
        }
        signal.notify_one();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(guard);
        return stack.empty();
    }

    bool tryPop(T& _value)
    {
        std::lock_guard<std::mutex> lock(guard);
        if (stack.empty())
        {
            return false;
        }

        _value = stack.top();
        stack.pop();
        return true;
    }

    bool waitAndPop(T& _value)
    {
        std::unique_lock<std::mutex> lock(guard);

	waiting_workers ++;

	if ((waiting_workers == total_workers) && stack.empty())
        {
	    /* We've exhausted all the work in the stack and all the workers are idle.  We're done. */
	    signal.notify_all();
	    return false;
        }

        while (stack.empty())
        {
            signal.wait(lock);
	    if ((waiting_workers == total_workers) && stack.empty()) {
	      return false;
	    }
        }

	waiting_workers --;

        _value = stack.top();
        stack.pop();
	return true;
    }

    std::stack<T>::size_type size(void)
    {
      std::lock_guard<std::mutex> lock(guard);
      return stack.size();
    }
private:
    std::stack<T> stack;
    mutable std::mutex guard;
    std::condition_variable signal;
    int total_workers;
    int waiting_workers;
};
