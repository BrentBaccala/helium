
/* GitHub Gist from https://gist.github.com/thelinked/6997598 */

/* Brent Baccala added the set_num_workers() method.  If the queue is empty, and the
 * number of threads waiting equals the value passed to set_num_workers(), then
 * all the calls to waitAndPop return false.  For their normal behavior, they return true.
 *
 * The idea is to detect when the queue is empty and all the workers are waiting on it.
 */

#pragma once
#include <queue>
#include <mutex>
#include <functional>
#include <condition_variable>

template<class T, class Compare = std::less<T>>
class LockingPriorityQueue
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
            queue.push(_data);
        }
        signal.notify_one();
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(guard);
        return queue.empty();
    }

    bool tryPop(T& _value)
    {
        std::lock_guard<std::mutex> lock(guard);
        if (queue.empty())
        {
            return false;
        }

        _value = queue.top();
        queue.pop();
        return true;
    }

    bool waitAndPop(T& _value)
    {
        std::unique_lock<std::mutex> lock(guard);

	waiting_workers ++;

	if ((waiting_workers == total_workers) && queue.empty())
        {
	    /* We've exhausted all the work in the queue and all the workers are idle.  We're done. */
	    signal.notify_all();
	    return false;
        }

        while (queue.empty())
        {
            signal.wait(lock);
	    if ((waiting_workers == total_workers) && queue.empty()) {
	      return false;
	    }
        }

	waiting_workers --;

        _value = queue.top();
        queue.pop();
	return true;
    }

    std::queue<T>::size_type size(void)
    {
      std::lock_guard<std::mutex> lock(guard);
      return queue.size();
    }
private:
    std::priority_queue<T, std::vector<T>, Compare> queue;
    mutable std::mutex guard;
    std::condition_variable signal;
    int total_workers;
    int waiting_workers;
};
