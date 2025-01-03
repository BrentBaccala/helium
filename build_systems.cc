/*
 * build_systems - optimized version of build_systems() from Sage/Python helium.sage
 *
 * BUILD: g++ -std=c++2a -march=native -O3 -o build_systems build_systems.cc -lpthread
 *
 * USAGE: ./build_systems NUM_THREADS < input-bit-strings > output-bit-strings
 *
 * GOAL: Given a set of polynomials (that form an ideal), we want to factor them all and
 * form a set of ideals, all formed from irreducible polynomials.
 *
 * PREPROCESS: Sage factors the polynomials and forms a list of factors.
 *
 * INPUT: a list of bit strings.  Each input bit string corresponds to a polynomial, and each
 * bit corresponds to a factor.
 *
 * OUTPUT: a list of bit strings.  Each output bit string corresponds to an ideal,
 * and the 1 bits are the factors in that ideal.
 *
 * See the print_build_systems() and load_systems() routines in helium.sage for generating
 * the input strings to this program, and parsing the output strings.
 *
 * All of the bit strings are of the same length (the number of factors).  We store
 * them in arrays of unsigned int (or some integer type specified as a typedef).
 *
 * ALGORITHM:
 *
 * Start with a bit string of all zeros (an empty ideal), a bit string of all ones (the allowed factors),
 * and an empty list of finished bit strings, but note that we might remove a bit string from the list of
 * finished bit strings, so they're not completely "finished".  (Or maybe not, as we don't seem to be
 * removing anything from the finished bit strings)
 *
 * We track the number of the polynomial we're working on.  First check to see if there's
 * any overlap between the polynomial's bit string and the working bit string (bitwise AND).
 * If so, this polynomial is satisfied; move on to the next one.
 *
 * Otherwise, we loop over the bits in polynomial's bit string that are also in the allowed
 * factors (logical AND).  If there are none, then we discard this case and backtrack.  For each allowed one bit, we form
 * a new bit string by adding that one bit to the working bit string.  This will produce
 * a set of bit strings that we keep in a list.  We remove the remaining bits
 * from the allowed bits, since there is no reason to consider them later, as these
 * factors will be considered when we backtrack to this point.
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
 *    - working ideal bit string
 *    - allowed bits bit string
 *    - polynomial number to process next
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
 *    - inverted_polys[i] are the bit-wise complements of polys[i]
 *    - expanded_polys[i] is a vector of bit strings, one vector for each polynomial, each bit string with a single bit set
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

#include <iostream>
#include <bitset>
#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include "LockingQueue.hpp"

// unsigned int bitstring_len = 0;

class BitString
{
public:
  typedef unsigned char data_type;
  std::vector<data_type> bitstring;
  typedef std::vector<data_type>::size_type size_type;
  unsigned int len;

  /* Bitstrings are represented MSB on the left, with the initial data_type in bitstring
   * using only as many bits as needed (LSBs) and the remaining data_type's fully populated.
   */

  BitString() {len = 0;}
  BitString(unsigned int l) : bitstring((l+8*sizeof(data_type)-1)/(8*sizeof(data_type))), len(l) {}
  BitString(std::string str)
  {
    len = str.length();
    auto bitstring_size = (len+8*sizeof(data_type)-1)/(8*sizeof(data_type));
    bitstring.resize(bitstring_size);
    for (int i=0; i < len; i++) {
      int val = (str[i] == '0' ? 0 : 1);
      int j = len - i - 1;
      int index = bitstring_size - j/(8*sizeof(data_type)) - 1;
      int offset = j%(8*sizeof(data_type));
      bitstring[index] |= val << offset;
    }
  }

  friend std::ostream& operator<<(std::ostream& stream, BitString bs);

  /* Bitwise NOT */

  BitString operator~() const
  {
    BitString *rv = new BitString(len);
    for (size_type i=0; i<bitstring.size(); i++) {
      rv->bitstring[i] = ~ bitstring[i];
    }
    return *rv;
  }

  /* Bitwise AND returning a BitString */

  BitString operator&(const BitString& rhs) const
  {
    BitString *rv = new BitString(len);
    for (size_type i=0; i<bitstring.size(); i++) {
      rv->bitstring[i] = bitstring[i] & rhs.bitstring[i];
    }
    return *rv;
  }

  /* Bitwise AND returning a boolean */

  bool operator&&(const BitString& rhs) const
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      if (bitstring[i] & rhs.bitstring[i])
	return true;
    }
    return false;
  }

  /* I don't want to copy return values, so use logical_or_assign instead. */
#if 0
  BitString operator|(const BitString& rhs) const
  {
    BitString *rv = new BitString(std::max(len, rhs.len));
    for (size_type i=0; i<std::max(bitstring.size(), rhs.bitstring.size()); i++) {
      if (i >= bitstring.size()) {
	rv->bitstring[i] = rhs.bitstring[i];
      } else if (i >= rhs.bitstring.size()) {
	rv->bitstring[i] = bitstring[i];
      } else {
	rv->bitstring[i] = bitstring[i] | rhs.bitstring[i];
      }
    }
    return *rv;
  }
#else
  BitString operator|(const BitString& rhs) const = delete;
#endif

  void logical_or_assign(BitString& lhs, const BitString& rhs) const
  {
    lhs.bitstring.resize(std::max(bitstring.size(), rhs.bitstring.size()));
    lhs.len = std::max(len, rhs.len);
    for (size_type i=0; i<std::max(bitstring.size(), rhs.bitstring.size()); i++) {
      if (i >= bitstring.size()) {
	lhs.bitstring[i] = rhs.bitstring[i];
      } else if (i >= rhs.bitstring.size()) {
	lhs.bitstring[i] = bitstring[i];
      } else {
	lhs.bitstring[i] = bitstring[i] | rhs.bitstring[i];
      }
    }
  }

  BitString operator|=(const BitString& rhs)
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] |= rhs.bitstring[i];
    }
    return *this;
  }

  BitString operator^=(const BitString& rhs)
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] ^= rhs.bitstring[i];
    }
    return *this;
  }

  bool is_superset_of(const BitString& rhs) const
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      // std::cerr << "superset test " << bitstring[i] << " " << rhs.bitstring[i] << "\n";
      if ((bitstring[i] & rhs.bitstring[i]) != rhs.bitstring[i]) {
	//std::cerr << " bitstring[i] & rhs.bitstring[i] " << (bitstring[i] & rhs.bitstring[i]) << " result " << ((bitstring[i] & rhs.bitstring[i]) != rhs.bitstring[i]) << "\n";
	return false;
      }
    }
    // std::cerr << *this << " is superset of " << rhs << "\n";
    return true;
  }
  operator bool() const
  {
    for (const auto i: bitstring) {
      if (i) return true;
    }
    return false;
  }
  BitString rightmost_set_bit(void) const
  {
    /* The return value has to be copied, but this routine is only called during initialization. */
    BitString result;
    result.bitstring.resize(bitstring.size());
    result.len = len;
    for (size_type i=0; i<bitstring.size(); i++) {
      data_type rightmost_set_bit = bitstring[i] & (-bitstring[i]);
      if (rightmost_set_bit) {
	result.bitstring[i] = rightmost_set_bit;
	return result;
      }
    }
    return result;
  }
};

std::ostream& operator<<(std::ostream& stream, BitString bs)
{
  for (int i=0; i < bs.bitstring.size(); i++) {
    auto str = std::bitset<sizeof(BitString::data_type)*8>(bs.bitstring[i]).to_string();
    if (i == 0) {
      auto bits_in_this_str = bs.len % (8*sizeof(BitString::data_type));
      if (bits_in_this_str > 0) {
	stream << str.substr(str.length() - bits_in_this_str, bits_in_this_str);
      } else {
	stream << str;
      }
    } else {
      stream << str;
    }
  }
  return stream;
}

struct BacktrackPoint
{
  BitString bitstring;
  BitString allowed_bits;
  unsigned int next_polynomial;
};

std::vector<BitString> polys;
std::vector<BitString> inverted_polys;
std::vector<std::vector<BitString>> expanded_polys;

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

    /* A final extra check for subsets, to avoid problems from the race condition mentioned below */
    for (const BitString& fbs: finished_bitstrings) {
      if (bitstring.is_superset_of(fbs))
	return;
    }

    // std::cerr << "add " << bitstring << "\n";
    /* remove any supersets */
    std::erase_if(finished_bitstrings, [&](const BitString& fbs) { return fbs.is_superset_of(bitstring); });

    finished_bitstrings.push_back(bitstring);
  }
};

FinishedBitStrings finished_bitstrings;

void task(void)
{
  /* "next work" is what this thread will do after "current work".  "extra work" has to be
   * pushed on the queue for later processing, either by this thread or by another.
   */
  BacktrackPoint current_work;
  BitString next_work_bitstring;
  BitString next_work_allowed_bits;
  BacktrackPoint extra_work;

  while (true) {
    get_next_work_from_queue:

    /* waitAndPop returns true if we've got work; false if all the work is done */
    if (! backtrack_queue.waitAndPop(current_work)) return;

    /* I put this here because while the work was on the queue, we could have added new finished bitstrings */
    /* There's probably a race condition here, as a new finished bitstring could appear at any time! */
    if (finished_bitstrings.contain_a_subset_of(current_work.bitstring)) {
      std::cerr << "detected superset\n";
      continue;
    }

    while (current_work.next_polynomial < polys.size()) {
      /* If there's any overlap here, the polynomial is already satisfied, move on */
      // std::cerr << "working on " << current_work.bitstring << "\n";
      if (! (current_work.bitstring && polys[current_work.next_polynomial])) {
	bool have_next_work = false;
	/* Extend backtrack queue if there's more than one factor(bit) in the polynomial */
	for (BitString next_bit: expanded_polys[current_work.next_polynomial]) {
	  // std::cerr << "adding " << next_bit << " " << current_work.allowed_bits << "\n";
	  /* skip this bit if it's not in the allowed_bits set */
	  if (next_bit && current_work.allowed_bits) {
	    if (! have_next_work) {
	      /* Like this, but avoids a copy */
	      /* next_work_bitstring = current_work.bitstring | next_bit; */
	      current_work.bitstring.logical_or_assign(next_work_bitstring, next_bit);
	      /* remove all of the current polynomial's bits from the allowed bits for future work */
	      next_work_allowed_bits = current_work.allowed_bits & inverted_polys[current_work.next_polynomial];
	      // std::cerr << "inverted_poly " << inverted_polys[current_work.next_polynomial] << "\n";
	      // std::cerr << "next_work_allowed_bits " << next_work_allowed_bits << "\n";
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      if (finished_bitstrings.contain_a_subset_of(next_work_bitstring)) {
		std::cerr << "detected superset\n";
		continue;
	      }
	      have_next_work = true;
	    } else {
	      /* Like this, but avoids a copy */
	      /* extra_work.bitstring = current_work.bitstring | next_bit; */
	      current_work.bitstring.logical_or_assign(extra_work.bitstring, next_bit);
	      extra_work.allowed_bits = current_work.allowed_bits & inverted_polys[current_work.next_polynomial];
	      extra_work.next_polynomial = current_work.next_polynomial + 1;
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      if (finished_bitstrings.contain_a_subset_of(extra_work.bitstring)) {
		std::cerr << "detected superset\n";
		continue;
	      }
	      backtrack_queue.push(extra_work);
	    }
	  }
	}
	if (! have_next_work) {
	  goto get_next_work_from_queue;
	}
	current_work.bitstring = next_work_bitstring;
	current_work.allowed_bits = next_work_allowed_bits;
      }
      current_work.next_polynomial ++;
    }

    finished_bitstrings.add(current_work.bitstring);
  }
}

int main(int argc, char ** argv)
{
  int nthreads = 1;
  int bitstring_len = 0;

  if (argc > 1) {
    nthreads = std::stoi(argv[1]);
  }

  std::string bitstring;
  while (std::getline(std::cin, bitstring)) {
    /* BitString's constructor can take a string */
    polys.emplace_back(bitstring);

    BitString bs(bitstring);
    if (bitstring_len == 0) {
      bitstring_len = bs.len;
    } else if (bitstring_len != bs.len) {
      std::cerr << "bitstring lengths aren't consistent\n";
      exit(1);
    }

    inverted_polys.emplace_back(~bs);

    /* Put an empty list at the back of expanded_polys and fill it with the individual bits */
    expanded_polys.emplace_back();
    BitString rmsb;
    do {
      rmsb = bs.rightmost_set_bit();
      expanded_polys.back().emplace_back(rmsb);
      bs ^= rmsb;
    } while (bs);
  }

  BacktrackPoint initial_work;
  initial_work.allowed_bits = ~ BitString(bitstring_len);
  // std::cerr << "allowed_bits " << initial_work.allowed_bits << "\n";
  initial_work.next_polynomial = 0;
  backtrack_queue.push(initial_work);
  backtrack_queue.set_num_workers(nthreads);

  std::vector<std::thread> threads;

  for (int i=0; i < nthreads; i++) {
    threads.emplace_back(task);
  }

  for (int i=0; i < nthreads; i++) {
    threads[i].join();
  }

  for (auto p:finished_bitstrings.finished_bitstrings) {
    std::cout << p << "\n";
  }

  return 0;
}
