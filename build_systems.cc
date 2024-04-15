/*
 * build_systems - optimized version of build_systems() from Sage/Python helium.sage
 *
 * BUILD: g++ -std=c++2a -o build_systems build_systems.cc
 *
 * GOAL: Given a set of polynomials (that form an ideal), we want to factor them all and
 * form a set of ideals, all formed from irreducible polynomials.
 *
 * PREPROCESS: Sage factors the polynomials and forms a list of factors.
 *
 * INPUT: a list of bit strings.  Each input bit string corresponds to a polynomial, and each
 * bit corresponds to a factor.
 *
 * OUTPUT: a list of bit strings.  Each output bit strings corresponds to an ideal,
 * and the 1 bits are the factors in that ideal.
 *
 * All of the bit strings are of the same length (the number of factors).  We store
 * them in arrays of unsigned int (or some integer type specified as a typedef).
 *
 * ALGORITHM:
 *
 * Start with a bit string of all zeros (an empty ideal), and an empty list of finished bit strings,
 * but note that we might remove a bit string from the list of finished bit strings,
 * so they're not completed "finished".
 *
 * We track the number of the polynomial we're working on.  First check to see if there's
 * any overlap between the polynomial's bit string and the working bit string (bitwise AND).
 * If so, this polynomial is satisfied; move on to the next one.
 *
 * Otherwise, we loop over the bits in polynomial's bit string.  For each one, we form
 * a new bit string by adding that one bit to the working bit string.  This will produce
 * a set of bit strings that we keep in a list.  Alternately, we track which bit we're
 * working on, since we need to backtrack.
 *
 * See https://stackoverflow.com/questions/68798413/blsi-instruction-isolate-lowest-set-bit
 *
 * After adding a bit to the working bit string, we now want to check to see if that
 * bit string is a superset of an existing finished bit string (bitwise AND and compare).
 * If so, we can discard it and move on.
 *
 * Otherwise, move on to the next polynomial.  Once we've hit the end of the polynomials,
 * we're about to add the working bit string to the finished bit strings.  First check to
 * see if any of the existing bit strings are supersets of the working bit string, if
 * so, remove them from the list of finished bit strings.  Add the working bit string
 * to the list of finished bit strings.
 *
 * Backtrack to the last point where we had multiple bits that we could have added to
 * the working bit string, and move on to the next one.
 *
 * For each point in the stack of backtrack points:
 *    - working bit string
 *    - polynomial number to process next
 *
 * If we don't include the bit number in the backtrack data, then we need to add a
 * set of backtrack points for each polynomial with multiple factors.  Doing it this
 * way probably makes it easier to implement multi-threading.
 *
 * Keep the stack of backtrack points as a queue.
 *
 * POSIX pthreads don't implement a queue.  Could just protect it with a mutex.
 * We want something a bit more complicated (maybe a condition variable) to wait
 * on the queue having data in it.
 *
 * Any running thread could potentially add more stuff to the queue, so a queue
 * empty condition does not suffice for termination.
 *
 * Ordering of the bit strings is unspecified.  Should we put all of the irreducible
 * polynomials first, because they will always be in the output?  The next set should
 * be any polynomials with those irreducibles as factors, because they will always
 * be satisfied.  After that, what, all the remaining polynomials with two factors,
 * just because we started with those with single factors?  This sorting should
 * be done in the pre-processing step; not going to bother with it in this code.
 *
 * POLYNOMIAL BIT STRINGS:
 *    - polys[i] are the complete bit strings for the i'th polynomial
 *    - first_bit[i] is the first one bit in the i'th polynomial
 *    - remaining_bits[i] is a vector of bit strings, one bit strings for each remaining bit split out
 *
 * FINISHED BIT STRINGS:
 *    - linked list of bit strings
 *    - allow multiple readers if there's no writer, because we're constantly checking
 *      this data structure in an attempt to prune our tree
 *    - if we remove finished bit strings, make sure this doesn't break readers
 *      partway through their read loop.  If we remove a finished bit string,
 *      we're going to add a subset of that string, so add it to the end of
 *      the list and make sure readers check the list in order.  What happens
 *      if we remove a bit string that a reader is holding a pointer to?
 *      Simplest might be to not allow any readers to be looping through the
 *      list while a writer holds the write lock.
 *    - We're checking (the read operation) for supersets of finished bit strings,
 *      so it might make sense to order the list by number of bits, then
 *      terminate the search once we've hit the first finished bit string with
 *      more bits than the bit string we're checking.
 *    - MRSW (multiple readers single writer) lock
 *    - PROTECTED with pthread_rwlock_t
 *
 * BACKTRACK QUEUE:
 *    - linked list of bit string and polynomial number to process next
 *    - threads wait on the queue until either it's got something in it, or all threads are waiting on
 *      an empty queue, in which case we're done, the threads exit, and we join them
 */

#include <vector>
#include <mutex>
#include <shared_mutex>
#include "LockingQueue.hpp"

unsigned int bitstring_len = 0;

class BitString
{
public:
  BitString operator&(const BitString& rhs) const
  {}
  BitString operator|(const BitString& rhs) const
  {}
  BitString operator|=(const BitString& rhs)
  {}
  bool is_superset_of(const BitString& rhs) const
  {}
  operator bool() const
  {}
};

struct BacktrackPoint
{
  BitString bitstring;
  unsigned int next_polynomial;
};

std::vector<BitString> polys;
std::vector<BitString> first_bit;
std::vector<std::vector<BitString>> remaining_bits;

LockingQueue<BacktrackPoint> backtrack_queue;

class FinishedBitStrings
{
public:
  std::vector<BitString> finished_bitstrings;
  std::shared_mutex mutex;

  bool contain_a_superset_of(const BitString& bitstring)
  {
    std::shared_lock<std::shared_mutex> lock(mutex);

    for (const BitString& fbs: finished_bitstrings) {
      if (fbs.is_superset_of(bitstring))
	return true;
    }
    return false;
  }

  bool contain_a_subset_of(const BitString& bitstring)
  {
    std::shared_lock<std::shared_mutex> lock(mutex);

    for (const BitString& fbs: finished_bitstrings) {
      if (bitstring.is_superset_of(fbs))
	return true;
    }
    return false;
  }

  void add(const BitString &bitstring)
  {
    std::unique_lock<std::shared_mutex> lock(mutex);

    /* remove any supersets */
    std::erase_if(finished_bitstrings, [&](const BitString& fbs) { return fbs.is_superset_of(bitstring); });

    finished_bitstrings.push_back(bitstring);
  }
};

FinishedBitStrings finished_bitstrings;

void task(void)
{
  BacktrackPoint current_work;

  while (true) {
    backtrack_queue.waitAndPop(current_work);

    while (current_work.next_polynomial < polys.size()) {
      /* If there's any overlap here, the polynomial is already satisfied, move on */
      if (! (current_work.bitstring & polys[current_work.next_polynomial])) {
	/* Extend backtrack queue if there's more than one factor(bit) in the polynomial */
	for (BitString next_bit: remaining_bits[current_work.next_polynomial]) {
	  BacktrackPoint new_work;
	  new_work.bitstring = current_work.bitstring | next_bit;
	  new_work.next_polynomial = current_work.next_polynomial + 1;
	  /* Check first if this is a superset of an existing bit string */
	  backtrack_queue.push(new_work);
	}
	/* Check first if this is a superset of an existing bit string */
	current_work.bitstring |= first_bit[current_work.next_polynomial];
	current_work.next_polynomial ++;

	if (finished_bitstrings.contain_a_subset_of(current_work.bitstring)) {}
      }
    }

    finished_bitstrings.add(current_work.bitstring);
  }
}

int main(int argc, char ** argv)
{
}
