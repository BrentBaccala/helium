/*
 * cnf2dnf - conjunctive normal form (CNF) to disjunctive normal form (DNF) converter
 *
 * BUILD: g++ -std=c++2a -march=native -O3 -o cnf2dnf cnf2dnf.cc -lpthread
 *
 * USAGE: ./cnf2dnf NUM_THREADS < input-bit-strings > output-bit-strings
 *
 * GOAL: Given an ideal defined by a set of polynomial generators, we want to factor the ideal
 * into an intersection of ideals.  We factor all of the polynomial generators, treat
 * each unique irreducible factor as a logical variable, and convert the resulting
 * logical system from CNF to DNF.  This will map to a set of ideals whose intersection
 * forms the original ideal.  This isn't yet a complete irreducible decomposition,
 * as each resulting ideal would have to be processed by something like the GTZ
 * algorithm to find its complete decomposition.
 *
 * PREPROCESS: Sage factors the polynomials and forms a list of factors.
 *
 * INPUT: a list of bit strings.  Each input bit string describes a cover, a list of logic
 * variables, any one of which can be true for the cover to be satisfied.  All of the
 * covers must be satisfied to satisfy the entire system.
 *
 * OUTPUT: a list of bit strings.  Each output bit string describes a product term, a
 * list of logic variables, all of which must be true for the product term to be satisfied.
 * Any of the product terms can be satisfied to satisfy the entire system.
 *
 * See the cnf2def_external() routine in helium.sage for generating
 * the input strings to this program, and parsing the output strings.
 *
 * All of the bit strings are of the same length (the number of factors).  We store
 * them in arrays of unsigned int (or some integer type specified as a typedef).
 *
 * ALGORITHM:
 *
 * CNF to DNF conversion is, in general, NP-complete.  We make heavy use of heuristics
 * that have proven themselves valuable for the specific case we're optimized for
 * (factoring ideals by factoring polynomials).
 *
 * We start by looking for single-bit covers, which are common because often a polynomial
 * does not factor at all.  All of these factors are required in the output, and any
 * polynomials that include those factors in their factorizations will be satisfied, too,
 * although this rarely happens for any complex ideal.
 *
 * We partition the remaining factors into covers, which are sets of factors.  Polynomials
 * not excluded in the first step either have all of their factors in a cover, or have
 * none of them in that cover.  The factors in the cover are expanded until this
 * condition is met.  Often there will only be a single cover for the entire system,
 * but sometimes it splits into two, and that's a big win, because they can be processed
 * independently and then combined together at the end.
 *
 * Within each cover, we look for a polynomials with two factors, one of which doesn't
 * appear in any other polynomial, which I call an "outlier".  This is quite a common case,
 * and can be optimized, since the presence of the outlier in the output systems depends
 * only on whether the other factor (the "attachment point") is present or not.
 *
 * Then, given a cover with all of its outliers excluded, we start with a bit string of all zeros
 * (an empty ideal) and an empty list of finished bit strings, but note that we might remove
 * a bit string from the list of finished bit strings, so they're not completely "finished".
 *
 * We move through the entire system, tracking the number of the polynomial we're working on.
 * First check to see if there's any overlap between the polynomial's bit string and the
 * working bit string (bitwise AND). If so, this polynomial is satisfied; move on to the next one.
 *
 * Otherwise, we loop over the bits in polynomial's bit string.
 * If there are none, then we discard this case and backtrack.  For each one bit, we form
 * a new bit string by adding that one bit to the working bit string.  This will produce
 * a set of bit strings that we keep in a list.
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
 *    - which cover we're working on
 *    - working ideal bit string
 *    - polynomial number to process next
 *    - some other things (list of finished polynomials, set of allowed_bits)
 *      that could be collapsed into the Cover, but just haven't been
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
 * POLYNOMIAL BIT STRINGS:
 *    - polys[i] are the complete bit strings for the i'th polynomial
 *    - inverted_polys[i] are the bit-wise complements of polys[i]
 *    - expanded_polys[i] is a vector of bit strings, one vector for each polynomial, each bit string with a single bit set
 *    - expanded_polys_bits[i] is a vector of ints, one vector for each polynomial, each a list of bit positions that are one in the polynomial
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
 *    - priority queue of bit string, polynomial number and index of factor within polynomial to process next
 *    - sorting function is standard numeric on the polynomial number; the largest polynomial number is always
 *      on top to be popped off.  It's "almost" a FILO stack, because each different thread will push
 *      items in increasing order until it backtracks onto the last one it pushed.  But collectively, the
 *      different threads might not push in increasing order, so we use a priority queue to keep the
 *      highest numbered polynomials on the top of the queue.  I don't think this is strictly needed
 *      for anything other than having a consistent organization of the backtrack queue to facilitate
 *      printing progress messages.
 *    - threads wait on the queue until either it's got something in it, or all threads are waiting on
 *      an empty queue, in which case we're done, the threads exit, and we join them
 */

#include <iostream>
#include <bitset>
#include <atomic>
#include <vector>
#include <list>
#include <map>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <csignal>
#include <cassert>
#include <unistd.h>
#include "LockingPriorityQueue.hpp"

// unsigned int bitstring_len = 0;

bool verbose = false;

// #define SLOW 1

struct {
  std::atomic<unsigned int> candidate_solutions;
  std::vector<std::atomic<unsigned int>> * subset_checks_made;
  std::vector<std::atomic<unsigned int>> * subset_checks_succeeded;
  std::vector<std::atomic<unsigned int>> * last_count_pushed;
  std::vector<std::atomic<unsigned int>> * last_count_popped;
} statistics;

class BitString
{
public:
  typedef unsigned long long int data_type;
  typedef std::vector<data_type>::size_type size_type;
  /* len - the number of bits in this BitString */
  unsigned int len;
  /* count() computes the number of one bits in a BitString.  It's used often enough
   * that we cache its value here (-1 if no value is cached).  "mutable" because
   * we need count() to be declared "const" so that it can be used in a std::multiset's
   * compare function, yet "cached_count" will be set in count() if it hasn't
   * been computed already.
   */
  mutable int cached_count;
  std::vector<data_type> bitstring;

  /* Bitstrings are represented MSB on the left, with the initial data_type in bitstring
   * using only as many bits as needed (LSBs) and the remaining data_type's fully populated.
   */

  BitString() {len = 0; cached_count = 0;}
  BitString(unsigned int l) : bitstring((l+8*sizeof(data_type)-1)/(8*sizeof(data_type))), len(l), cached_count(0) {}
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
      bitstring[index] |= data_type(val) << offset;
    }
    cached_count = -1;
  }

  friend std::ostream& operator<<(std::ostream& stream, BitString bs);

  void set_bit(int bit) {
    int index_in_bitstring = bitstring.size() - bit/(8*sizeof(BitString::data_type)) - 1;
    int bit_in_integer = bit%(8*sizeof(BitString::data_type));
    bitstring[index_in_bitstring] |= (BitString::data_type(1) << bit_in_integer);
  }

  void clear_bit(int bit) {
    int index_in_bitstring = bitstring.size() - bit/(8*sizeof(BitString::data_type)) - 1;
    int bit_in_integer = bit%(8*sizeof(BitString::data_type));
    bitstring[index_in_bitstring] &= ~(BitString::data_type(1) << bit_in_integer);
  }

  bool test_bit(int bit) {
    int index_in_bitstring = bitstring.size() - bit/(8*sizeof(BitString::data_type)) - 1;
    int bit_in_integer = bit%(8*sizeof(BitString::data_type));
    return ((bitstring[index_in_bitstring] & (BitString::data_type(1) << bit_in_integer)) != 0);
  }

  /* Bitwise NOT */

  BitString operator~() const
  {
    BitString rv(len);
    for (size_type i=0; i<bitstring.size(); i++) {
      rv.bitstring[i] = ~ bitstring[i];
    }
    rv.cached_count = -1;
    return BitString(rv);
  }

  /* Bitwise equality */

  bool operator==(const BitString& rhs) const
  {
    if (len != rhs.len) return false;
    for (size_type i=0; i<bitstring.size(); i++) {
      if (bitstring[i] != rhs.bitstring[i])
	return false;
    }
    return true;
  }

  /* less-than comparison so we can use this class as the Key in a std::map
   * (I'm not sure why the default operator< doesn't work, but it doesn't)
   *
   * Treats the BitString as one giant unsigned integer.
   */

  bool operator<(const BitString& rhs) const
  {
    if (len != rhs.len) return (len < rhs.len);
    for (size_type i=0; i<bitstring.size(); i++) {
      if (bitstring[i] != rhs.bitstring[i])
	return (bitstring[i] < rhs.bitstring[i]);
    }
    return false;
  }

  /* Bitwise AND returning a BitString */

  BitString operator&(const BitString& rhs) const
  {
    BitString rv(len);
    for (size_type i=0; i<bitstring.size(); i++) {
      rv.bitstring[i] = bitstring[i] & rhs.bitstring[i];
    }
    rv.cached_count = -1;
    return BitString(rv);
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
    lhs.cached_count = -1;
  }

  BitString operator|=(const BitString& rhs)
  {
    if (len != rhs.len) throw std::runtime_error("incompatible lengths in operator|=");
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] |= rhs.bitstring[i];
    }
    cached_count = -1;
    return *this;
  }

  BitString operator&=(const BitString& rhs)
  {
    if (len != rhs.len) throw std::runtime_error("incompatible lengths in operator&=");
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] &= rhs.bitstring[i];
    }
    cached_count = -1;
    return *this;
  }

  BitString operator^=(const BitString& rhs)
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] ^= rhs.bitstring[i];
    }
    cached_count = -1;
    return *this;
  }

  /* superset test - not strict; calling this method on self will return true */
  bool is_superset_of(const BitString& rhs) const
  {
    if (count() < rhs.count()) return false;
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
    result.cached_count = 0;
    for (size_type i=bitstring.size()-1; i>=0; i--) {
      data_type rightmost_set_bit = bitstring[i] & (-bitstring[i]);
      if (rightmost_set_bit) {
	result.bitstring[i] = rightmost_set_bit;
	result.cached_count = 1;
	return result;
      }
    }
    return result;
  }

  int rightmost_set_bit_index(void) const
  {
    for (int i=bitstring.size()-1; i>=0; i--) {
      for (auto j=0; j<8*sizeof(data_type); j++) {
	if (bitstring[i] & (data_type(1) << j)) {
	  return (bitstring.size() - i - 1)*8*sizeof(data_type) + j;
	}
      }
    }
    return -1;
  }

  /* Simple (but inefficient) method for iterating over all one bits in a BitString */
  std::vector<BitString> all_one_bits(void) const
  {
    std::vector<BitString> result;
    BitString bs = *this;
    BitString rmsb;
    do {
      rmsb = bs.rightmost_set_bit();
      result.emplace_back(rmsb);
      bs ^= rmsb;
    } while (bs);
    return result;
  }

  int count(void) const
  {
    if (cached_count == -1) {
      int value = 0;
      for (int i=0; i<bitstring.size(); i++) {
	value += std::popcount(bitstring[i]);
      }
      cached_count = value;
    }
    return cached_count;
  }

  // Iterator to loop over all one bits in a bitstring

  class bititerator {

  private:
    BitString * bitstring;
    int bit;

  public:

    // These are here to ensure that we can use a bititerator to initialize a std::vector<int> (for testing)
    using difference_type = std::ptrdiff_t;
    using value_type = int;
    using pointer = const int*;
    using reference = const int&;
    using iterator_category = std::input_iterator_tag;

    // Constructor used for begin iterator (bit starts at rightmost one bit)
    explicit bititerator(BitString * ptr) : bitstring(ptr) {
      bit = ptr->rightmost_set_bit_index();
    }

    // Constructor used for end iterator (bit explicitly set to -1)
    explicit bititerator(BitString * ptr, int bit) : bitstring(ptr), bit(bit) {
    }

    // Prefix increment (advance bit to next one bit, or -1 if there is none)
    bititerator& operator++() {
      for (int i=bitstring->bitstring.size()-(bit+1)/(8*sizeof(BitString::data_type))-1; i>=0; i--) {
	for (auto j=(bit+1)%(8*sizeof(BitString::data_type)); j<8*sizeof(BitString::data_type); j++) {
	  if (bitstring->bitstring[i] & (BitString::data_type(1) << j)) {
	    bit = (bitstring->bitstring.size() - i - 1)*8*sizeof(BitString::data_type) + j;
	    return *this;
	  }
	}
      }
      bit = -1;
      return *this;
    }

    // Dereference operator
    //
    // Used by a C++ range-based for loop to construct the loop item
    int operator*() const {
      return bit;
    }

    // Equality operator
    bool operator==(const bititerator& other) const {
      return (bitstring == other.bitstring) && (bit == other.bit);
    }

    // Inequality operator
    bool operator!=(const bititerator& other) const {
      return !(*this == other);
    }

  };

  // Begin iterator to loop over all one bits in a bitstring
  bititerator begin() {
    return bititerator(this);
  }

  // End iterator
  bititerator end() {
    return bititerator(this, -1);
  }
};

std::ostream& operator<<(std::ostream& stream, BitString bs)
{
  for (int i=0; i < bs.bitstring.size(); i++) {
    auto str = std::bitset<8*sizeof(BitString::data_type)>(bs.bitstring[i]).to_string();
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

std::vector<BitString> polys;
std::vector<BitString> inverted_polys;
std::vector<std::vector<BitString>> expanded_polys;
std::vector<std::vector<int>> expanded_polys_bits;

/* All bits for which there is a polynomial containing only that bit.  Such bits must
 * always be present in the output bit strings, and any polynomials that depend on them
 * can be removed from under_consideration, since they will always be satisfied.
 */
BitString single_bit_covers;
int single_bit_polynomials_covered = 0;

/* A basic linked list of BitString's, with a twist: we want to do thread-safe list operations
 * without lock synchronization delays.  Inserting at the head of the list isn't a problem
 * (we'll just use atomics in FinishedBitStrings), but deleting is more difficult.
 * The simplest solution: don't delete, just flag nodes invalid so they get skipped
 * during iteration.
 */

class LinkedBitString : public BitString {
  public:
  LinkedBitString * next = nullptr;
  std::atomic<bool> valid = true;

  LinkedBitString(const BitString& bs) : BitString(bs) {}

  // Nested iterator class (written mostly by GPT-4o)
  class iterator {
  public:
    // Iterator traits
    using iterator_category = std::forward_iterator_tag;
    using value_type = LinkedBitString;
    using difference_type = std::ptrdiff_t;
    using pointer = LinkedBitString*;
    using reference = LinkedBitString&;

    // Constructor
    explicit iterator(pointer ptr = nullptr) : current(ptr) {
      // do this in case the list starts with one or more invalid nodes
      while (current && ! current->valid) {
	current = current->next;
      }
    }

    // Dereference operator
    reference operator*() const {
      return *current;
    }

    pointer operator->() const {
      return current;
    }

    // Prefix increment
    iterator& operator++() {
      if (current) {
        current = current->next; // Move to the next node
	while (current && ! current->valid) {
	  current = current->next;  // skip invalid nodes
	}
      }
      return *this;
    }

    // Postfix increment
    iterator operator++(int) {
      if (current) {
        current = current->next; // Move to the next node
      }
      return *this;
    }

    // Equality operator
    bool operator==(const iterator& other) const {
      return current == other.current;
    }

    // Inequality operator
    bool operator!=(const iterator& other) const {
      return !(*this == other);
    }

  private:
    pointer current; // Pointer to the current node
  };

  // Begin iterator
  iterator begin() {
    return iterator(this);
  }

  // End iterator (nullptr represents the end)
  iterator end() {
    return iterator(nullptr);
  }
};

class FinishedBitStrings
{
public:
  std::vector<std::atomic<LinkedBitString *>> by_count;
  std::atomic<int> valid_bitstrings = 0;
  std::atomic<int> invalid_bitstrings = 0;

  FinishedBitStrings(int max_count) : by_count(max_count+1) { }

  ~FinishedBitStrings()
  {
    for (auto& albs: by_count) {
      auto lbs = albs.load();
      while (lbs) {
	auto lbsnext = lbs->next;
	delete lbs;
	lbs = lbsnext;
      }
    }
  }

  bool contain_a_subset_of(const BitString& bitstring)
  {
    auto count = bitstring.count();

    /* We're checking for subsets of bitstring, so we only need to check
     * finished_bitstrings with fewer bits than bitstring.
     */

    for (auto i=0; i<count; i++) {
      if (by_count[i]) {
	for (auto& fbs: *by_count[i]) {
	  if (bitstring.is_superset_of(fbs))
	    return true;
	}
      }
    }
    return false;
  }

  void add(const BitString &bitstring)
  {
    /* This will make a copy of bitstring, which is what we want, because the current_work.bitstring
     * that was passed in (by reference) will get copy-assigned when we pull the next work item
     * off of backtrack_queue.
     */
    LinkedBitString * new_node = new LinkedBitString(bitstring);
    auto count = bitstring.count();

    /* This code is cribbed from: https://en.cppreference.com/w/cpp/atomic/atomic/compare_exchange
     *
     * We're going to put new_node on the top of the linked list after checking to make sure
     * that it doesn't duplicate anything already on the linked list.  If the head isn't
     * what's stored in new_node->next (either because it's the first time through the loop
     * or because some other thread inserted a node just now) we check for duplicates
     * and try again.  Keep track of the first node on the list that we checked for a duplicate
     * so we don't check twice if we have to run the loop multiple times.
     *
     * If there's anything at all on the list, the compare_exchange_weak will fail the first
     * time through the loop (because new_node->next is NULL and by_count[count] isn't NULL),
     * and we'll check everything on the list for duplicates (because last_checked_node is NULL).
     * Further failures of compare_exchange_weak will only check newly added nodes.
     */

    new_node->next = NULL;
    LinkedBitString * last_checked_node = NULL;

    while (!by_count[count].compare_exchange_weak(new_node->next, new_node,
						  std::memory_order_release,
						  std::memory_order_relaxed)) {
      for (auto check_node = new_node->next; check_node != last_checked_node; check_node = check_node->next) {
	if (*check_node == bitstring) {
	  delete new_node;
	  return;
	}
      }
      last_checked_node = new_node->next;
    }

    valid_bitstrings ++;

#ifdef DEBUG
    /* Let's run a final check for subsets that were added after our last subset check.
     * We can no longer delete the node, but we can invalidate it.
     *
     * We do have to do the entire atomic compare_exchange_strong bit because we can't
     * assume that new_node is still valid!  After all, another thread could be in
     * the next loop below this one and invalidate it.
     */
    for (auto i=0; i<count; i++) {
      if (by_count[i]) {
	for (auto& fbs: *by_count[i]) {
	  if (new_node->valid && bitstring.is_superset_of(fbs)) {
	    bool true_value = true;
	    if (new_node->valid.compare_exchange_strong(true_value, false)) {
	      valid_bitstrings --;
	      invalid_bitstrings ++;
	    }
	    return;
	  }
	}
      }
    }

    /* Now let's run over finished bitstrings with more bits and invalidate those that are supersets of the new one */
    for (count ++; count < by_count.size(); count ++) {
      if (by_count[count]) {
	for (auto& fbs: *by_count[count]) {
	  if (fbs.valid && fbs.is_superset_of(bitstring)) {
	    bool true_value = true;
	    if (fbs.valid.compare_exchange_strong(true_value, false)) {
	      valid_bitstrings --;
	      invalid_bitstrings ++;
	    }
	  }
	}
      }
    }
#endif
  }

  // Nested iterator Class (written mostly by GPT-4o)
  class iterator
  {
  public:
    // Type aliases for standard iterator traits
    using iterator_category = std::forward_iterator_tag;
    using value_type = LinkedBitString;
    using difference_type = std::ptrdiff_t;
    using pointer = LinkedBitString*;
    using reference = LinkedBitString&;

  private:
    FinishedBitStrings* finishedBitStrings;
    size_t outerIndex;                         // Index in `by_count`
    LinkedBitString::iterator innerIterator;   // Iterator for the current LinkedBitString's data
    LinkedBitString::iterator innerEnd;        // End iterator for the current LinkedBitString's data

    // Helper to advance to the next valid LinkedBitString and set up its iterator
    void advanceToNextValid()
    {
      while (outerIndex < finishedBitStrings->by_count.size()) {
	LinkedBitString* current = finishedBitStrings->by_count[outerIndex].load();
	if (current) {
	  innerIterator = current->begin();
	  innerEnd = current->end();
	  if (innerIterator != innerEnd) {
	    return; // Found a valid LinkedBitString with data to iterate over
	  }
	}
	++outerIndex; // Move to the next element in by_count
      }
    }

  public:
    // Constructor
    iterator(FinishedBitStrings* finishedBitStrings, size_t start)
      : finishedBitStrings(finishedBitStrings), outerIndex(start)
    {
      advanceToNextValid(); // Ensure we start at the first valid element
    }

    // Dereference operator
    reference operator*() { return *innerIterator; }

    // Arrow operator
    pointer operator->() { return &(*innerIterator); }

    // Pre-increment operator
    iterator& operator++()
    {
      ++innerIterator;
      if (innerIterator == innerEnd) {
	++outerIndex;
	advanceToNextValid();
      }
      return *this;
    }

    // Post-increment operator
    iterator operator++(int)
    {
      iterator temp = *this;
      ++(*this);
      return temp;
    }

    // Equality operator
    bool operator==(const iterator& other) const
    {
      return outerIndex == other.outerIndex &&
	innerIterator == other.innerIterator &&
	finishedBitStrings == other.finishedBitStrings;
    }

    // Inequality operator
    bool operator!=(const iterator& other) const { return !(*this == other); }
  };

  // Begin and End methods for the iterator
  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, by_count.size()); }
};

/* Covers
 *
 * We partition the input bits into sets that completely cover a group of polynomials,
 * and display the number of bits in the set and the number of polynomials covered.
 *
 * Smaller partitions are easier to handle: 1 is a special case, and 2 through 5
 * can be handled with lookup tables (not implemented, because they don't appear in practice).
 *
 */

class Cover
{
  public:
  BitString cover;                                       /* which bits form the cover */
  std::vector<bool> under_consideration;                 /* which polynomials are candidates */
  std::map<BitString, BitString> single_link_chains;
  int triplets;
  FinishedBitStrings finished_bitstrings;

  /* FinishedBitStrings is a bit difficult to initialize because it includes a vector of atomics,
   * which can't be resized, so we need to know its size at initialization time.  We use a trick
   * recommended here: https://stackoverflow.com/a/61033668/1493790
   * which is to use a delegating constructor to calculate the cover in the initializer list.
   *
   * expand_cover: Given a cover and an initial polynomial "under consideration", expand the
   * cover to include any polynomials under consideration that partially match the cover.
   * Once we're done, all of the polynomials "under consideration" are either completely
   * under the (expanded) cover, or are not under it at all.
   */

  static std::pair<std::vector<bool>, BitString> expand_cover(BitString starting_poly)
  {
    bool expanding_cover;
    std::vector<bool> under_consideration(polys.size(), false);
    BitString cover = starting_poly;
    do {
      expanding_cover = false;
      for (auto i=0; i<polys.size(); i++) {
	if (! under_consideration[i] && ! (single_bit_covers && polys[i]) && (cover && polys[i])) {
	  under_consideration[i] = true;
	  cover |= polys[i];
	  expanding_cover = true;
	}
      }
    } while (expanding_cover);
    return make_pair(under_consideration, cover);
  }

  Cover(BitString starting_poly) : Cover(expand_cover(starting_poly)) { }

  Cover(std::pair<std::vector<bool>, BitString> pair) : cover(pair.second), under_consideration(pair.first), finished_bitstrings(cover.count())
  {
    /* Identify which polynomials are links (2-bit polynomials).
     *
     * If so, there are several possibilities: either it's a "outlier"
     * (the case we want to optimize), or the link is itself the entire cover, or it's a link
     * between two subcovers.  The link has two sides (left and right).  Count up how many polynomials
     * match up with each.  An outlier has one side matching a single polynomial (the link itself)
     * and the other side matches more than one.  A link between two subcovers has both sides
     * matching more than one polynomial.  If both sides match a single polynomial, then the
     * link is the entire cover.
     */

    for (auto i=0; i<polys.size(); i++) {
      if (under_consideration[i] && (polys[i] && cover) && (polys[i].count() == 2)) {
	BitString outlying_point;
	BitString attachment_point;
	for (BitString current_bit: expanded_polys[i]) {
	  int matching_polys = 0;
	  for (auto j=0; j<polys.size(); j++) {
	    if (under_consideration[j] && (polys[j] && current_bit)) matching_polys ++;
	  }
	  if (matching_polys == 1) {
	    outlying_point = current_bit;
	  } else {
	    attachment_point = current_bit;
	  }
	}
	if ((attachment_point.count() > 0) && (outlying_point.count() > 0)) {
	  if (! single_link_chains.contains(attachment_point)) {
	    single_link_chains[attachment_point] = BitString(polys[0].len);
	  }
	  single_link_chains[attachment_point] |= outlying_point;
	}
      }
    }

    /* Count "triplets": Polynomials with two isolated bits and only
     * one touching the rest of the cover.
     */

    triplets = 0;
    for (auto i=0; i<polys.size(); i++) {
      if (under_consideration[i] && (polys[i] && cover) && (polys[i].count() == 3)) {
	int single_poly_bits = 0;
	for (BitString next_bit: expanded_polys[i]) {
	  int matching_polys = 0;
	  for (auto j=0; j<polys.size(); j++) {
	    if (under_consideration[j] && (polys[j] && next_bit)) matching_polys ++;
	  }
	  if (matching_polys == 1) single_poly_bits ++;
	}
	if (single_poly_bits == 2) triplets ++;
      }
    }
  }
};

struct BacktrackPoint
{
  BitString bitstring;
  BitString allowed_bits;
  unsigned int next_index;
  unsigned int next_polynomial;
  Cover * cover;
};

struct BacktrackQueueLess
{
  bool operator()(const BacktrackPoint& l, const BacktrackPoint& r) const { return l.next_polynomial < r.next_polynomial; }
};

LockingPriorityQueue<BacktrackPoint, BacktrackQueueLess> backtrack_queue;

std::list<Cover> all_covers;

/* task() - the main processing routine.  Should be called with at least one
 * BacktrackPoint (probably polynomial 0, index 0) already pushed onto backtrack_queue.
 * Multiple threads can call task() even if only one BacktrackPoint is pushed;
 * it's designed to wait until a BacktrackPoint is available, unless backtrack_queue
 * detects that there is no more work to be done.
 */

void task(void)
{
  /* "current work" is what we're doing now.  "extra work" has to be pushed on
   * the backtrack queue for later processing, either by this thread or by another.
   */
  BacktrackPoint current_work;
  BacktrackPoint extra_work;

  while (true) {
    get_next_work_from_queue:

    /* waitAndPop returns true if we've got work; false if all the work is done */
    if (! backtrack_queue.waitAndPop(current_work)) return;

    (*statistics.last_count_popped)[current_work.next_polynomial] = statistics.candidate_solutions.load();
    (*statistics.last_count_pushed)[current_work.next_polynomial] = 0;

    /* I used to have a test here: finished_bitstrings.contain_a_subset_of(current_work.bitstring),
     * because even though we tested for this before adding the work to the backtrack_queue,
     * we could have added new finished bitstrings while the work was on the queue.  The program
     * spent a lot of time here, and I decided that we don't really need this check because
     * we're going to check contain_a_subset_of if we either add work to the work queue,
     * or add a bitstring to finished_bitstrings.  So the only cost of not doing a
     * subset test at this point is avoiding a single additional pass through this
     * while loop, and that doesn't seem to be worth it.
     *
     * Update: I'm not checking before we add to the backtrack queue currently, because
     * I'm no longer adding the bit to the bitstring before pushing on the backtrack queue.
     * I could add the bit in a temporary just to check if we should add it to the
     * backtrack queue, or maybe it makes more sense to check it here (when we pull it
     * off the backtrack queue), but I haven't put the code back in yet.
     *
     * It's an optimization, because the FinishedBitStrings add() method makes sure
     * that any supersets get invalidated in favor of their subsets.
     */

    /* current_work contains next_polynomial, and the bitstring is obtained from the polys[] array.
     *
     * We could instead put next_polynomial in current_work.  So we pull out next_polynomial
     * and polynomial_bitstring (a truncated version) from current_work.  Or we pull out
     * an index into expanded_polys, which is a std::vector of std::vector<BitString>s.
     *
     * We have next_polynomial and next_index.  We loop next_bit over expanded_polys[next_polynomial]
     * from next_index to its last entry.
     *
     * for (BitString next_bit) loop:
     *   needs to start at zero most of the time, but if we're just coming in from
     *     a new BacktrackPoint we start at current_work.next_index
     *   So, at the end of the loop, where we increment next_polynomial, we set next_index to zero
     *   The first next_bit we find,
     *      we set extra_work's bitstring to current_work.bitstring
     *      we update current_work.bitstring with an assignment, and continue to see if there's a second next_bit
     *   The second next_bit we find,
     *      we save the unmodified bitstring to extra_work (that's why we set it above)
     *      along with the allowed_bits, cover, and next_index (the index of the second next_bit, as we haven't processed it yet)
     *      push extra_work (it's a BacktrackPoint) to the backtrack_queue
     *      break out of the loop and keep going to process the first next_bit
     *   Another additional next_bit's (three or more) will get processed when the BacktrackPoint is processed.
     *      At that point, the second next_bit will now be the first next_bit, so the third next_bit would become the second next_bit
     *      and yet another BacktrackPoint would be created.
     *   If we didn't find any next_bit's, we goto get_next_work_from_queue
     *   Otherwise, keep going with current_work and move on to the next polynomial
     */

    while (current_work.next_polynomial < polys.size()) {
      /* If there's any overlap here, the polynomial is already satisfied, move on */
      // std::cerr << "working on " << current_work.bitstring << "\n";
      if (current_work.cover->under_consideration[current_work.next_polynomial] && ! (current_work.bitstring && polys[current_work.next_polynomial])) {
	bool have_next_work = false;
	/* Extend backtrack queue if there's more than one factor(bit) in the polynomial */
	/* XXX should this loop use a reference for speed? */
	unsigned int initial_index = current_work.next_index;
	for (; current_work.next_index < expanded_polys_bits[current_work.next_polynomial].size(); current_work.next_index ++) {
	  int next_bit = expanded_polys_bits[current_work.next_polynomial][current_work.next_index];
	  /* skip this bit if it's not in the allowed bit set */
	  if (current_work.allowed_bits.test_bit(next_bit)) {
	    // std::cerr << "adding " << next_bit << "\n";
	    if (! have_next_work) {
	      extra_work.bitstring = current_work.bitstring;
	      current_work.bitstring.set_bit(next_bit);
#ifdef SLOW
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      /* XXX this check is optional and time consuming, so we should be more clever about how often we do this */
	      (*statistics.subset_checks_made)[current_work.next_polynomial] ++;
	      if (current_work.cover->finished_bitstrings.contain_a_subset_of(current_work.bitstring)) {
		/* put things back the way they were, and skip this solution */
		(*statistics.subset_checks_succeeded)[current_work.next_polynomial] ++;
		current_work.bitstring = extra_work.bitstring;
		continue;
	      }
#endif
	      /* If the current polynomial's bits are 111, we want to create future work 1xx, 01x, 001,
	       * (not 1xx, x1x, xx1), so we now remove the current bit from the allowed bits.
	       */
	      current_work.allowed_bits.clear_bit(next_bit);

	      /* Check for forced assignments
	       *
	       * All I'm interested in for current purposes are the 01x type assignments, that actually
	       * assign a zero somewhere, since that's the only thing that can trigger forced assignments,
	       * as we have a monotone function, with no negations in the CNF.
	       *
	       * Monotone function: only zeros can trigger forced assignments, and all forced assignments are ones.
	       *
	       * If initial_index is zero, then this is the first bit we're assigning for this polynomial,
	       * and there are no forced assignments, since we're assigning something like "1xx".
	       *
	       * Otherwise (initial_work != 0), this is a later bit in the polynomial, we're assigning
	       * something like "01x" or "001", and the zeros might trigger forced assignments.
	       *
	       * Due to the monotone nature of our CNF, the forced assignments will always be one bits,
	       * so there's no need to loop or recurse here to check for additional forced assignments.
	       *
	       * We want to put the forced ones into current_work, by setting them in bitstring
	       * and clearing them in allowed_bits.
	       *
	       * All polynomials earlier than current_work.next_polynomial have already been satisfied,
	       * so we start our check with current_work.next_polynomial + 1.
	       */

	      /* XXX this is another optional step that we should collect statistics on
	       * If we didn't do this on every pass, it would still pick up all of the
	       * forced assignments next time; i.e, it won't miss any.
	       */

	      if (initial_index != 0) {
		for (int check_poly = current_work.next_polynomial + 1; check_poly < polys.size(); check_poly ++) {
		  /* only check polys in our Cover */
		  if (current_work.cover->under_consideration[check_poly]) {
		    /* if it isn't already satisfied... */
		    if (! (current_work.bitstring && polys[check_poly])) {
		      /* ...and there's only one bit left than can satisfy it... */
		      BitString possible_bits = current_work.allowed_bits & polys[check_poly];
		      if (possible_bits.count() == 1) {
			/* ...then that bit is a forced assignment */
			current_work.bitstring |= possible_bits;
			current_work.allowed_bits &= ~possible_bits;
		      }
		    }
		  }
		}
	      }
	      have_next_work = true;
	    } else {
	      /* extra_work.bitstring was already set above the first time through this loop */
	      extra_work.allowed_bits = current_work.allowed_bits;
	      extra_work.next_index = current_work.next_index;
	      extra_work.next_polynomial = current_work.next_polynomial;
	      extra_work.cover = current_work.cover;
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      /* if (current_work.cover->finished_bitstrings.contain_a_subset_of(extra_work.bitstring)) continue; */
	      (*statistics.last_count_pushed)[extra_work.next_polynomial] = statistics.candidate_solutions.load() + 1;
	      backtrack_queue.push(extra_work);
	      break;
	    }
	  }
	}
	if (! have_next_work) {
	  goto get_next_work_from_queue;
	}
      }
      current_work.next_polynomial ++;
      current_work.next_index = 0;
    }
    /* We've now got a bitstring that's in the DNF.  See if any of its subsets are also in the DNF,
     * by clearing each one bit and checking to see if the bitstring still works.  Does the order
     * of the bits affect the result?  Yes, but we don't expect the final result of the entire
     * program to depend on exactly which subset is selected at this point, since we'll eventually
     * find all solutions.
     */
    for (auto bit: current_work.bitstring) {
      current_work.bitstring.clear_bit(bit);
      int i;
      for (i=0; i<polys.size(); i++) {
	if (current_work.cover->under_consideration[i] && ! (polys[i] && current_work.bitstring)) break;
      }
      if (i != polys.size()) {
	current_work.bitstring.set_bit(bit);
      }
    }
    /* We've now got a prime bitstring in the DNF.  Add it to the finished DNF. */
    /* add() will check first to see if the bitstring is already in finished_bitstrings */
    statistics.candidate_solutions ++;
    current_work.cover->finished_bitstrings.add(current_work.bitstring);
  }
}

/* Remove from consideration all polynomials with only a single bit set, and all others that depend on them.
 *
 * Modifies global variable 'single_bit_covers'
 */

void compute_and_remove_single_bit_covers(void)
{
  single_bit_covers = BitString(polys[0].len);

  std::vector<bool> under_consideration(polys.size(), true);

  for (int i=0; i<polys.size(); i++) {
    if (under_consideration[i] && (polys[i].count() == 1)) {
      single_bit_polynomials_covered ++;
      under_consideration[i] = false;
      single_bit_covers |= polys[i];
    }
  }
  for (int i=0; i<polys.size(); i++) {
    if (under_consideration[i] && (single_bit_covers && polys[i])) {
      single_bit_polynomials_covered ++;
      under_consideration[i] = false;
    }
  }
}

void compute_all_covers(void)
{
  /* For each remaining polynomial, form a cover that is initially just the first polynomial.
   * Then loop over all the polynomials, OR-ing in any polynomials that match the cover.  Keep looping
   * until the cover stabilizes.  Remove all of these polynomials from consideration, and keep doing
   * it until we've completely partitioned the input set.
   */
  BitString union_of_all_covers = single_bit_covers;

  while (true) {
    int polys_covered;
    int i;
    for (i=0; i<polys.size(); i++) {
      if (! (union_of_all_covers && polys[i])) {
	break;
      }
    }
    if (i == polys.size()) {
      break;
    }

    Cover& cover = all_covers.emplace_back(polys[i]);
    union_of_all_covers |= cover.cover;
  }

}

/* Compute and display some statistics about the input data */

void compute_and_display_statistics(void)
{
  if (single_bit_covers.count() == 1) std::cerr << "1 single bit cover covering ";
  else std::cerr << single_bit_covers.count() << " single bit covers covering ";
  if (single_bit_polynomials_covered == 1) std::cerr << "1 polynomial\n";
  else std::cerr << single_bit_polynomials_covered << " polynomials\n";

  for (auto& cover: all_covers) {
    int outliers = 0;
    for (auto const &[attachment_point, bs]: cover.single_link_chains) {
      outliers += bs.count();
    }
    std::cerr << "a " << cover.cover.count() << "-bit cover; ";
    if (outliers == 0) std::cerr << "0 outliers";
    else {
      if (outliers == 1) std::cerr << "1 outlier with ";
      else std::cerr << outliers << " outliers with ";
      if (cover.single_link_chains.size() == 1) std::cerr << "1 attachment point";
      else std::cerr << cover.single_link_chains.size() << " attachment points";
    }

    if (cover.triplets == 1) std::cerr << "; 1 triplet\n";
    else std::cerr << "; " << cover.triplets << " triplets\n";
    /* missing: print number of polynomials covered and max length of chain (if not 1) */
  }

}

// Signal handler
void signalHandler(int signum) {
    if (signum == SIGQUIT) {
        std::cerr << "Backtrack queue: " << backtrack_queue.size() << "\n";
        std::cerr << "Finished bitstrings: ";
	for (auto &cover: all_covers) {
	  std::cerr << cover.finished_bitstrings.valid_bitstrings << " (plus " << cover.finished_bitstrings.invalid_bitstrings << " invalid) ";
	}
	std::cerr << "\n";
	std::cerr << "Candidate solutions: " << statistics.candidate_solutions << "\n";
#if 0
	/* One way of printing the last popped counts: print them all out */
	for (int i = 0; i < polys.size(); i ++) {
	  // std::cerr << (*statistics.last_count_pushed)[i] << " " << (*statistics.last_count_popped)[i] << " ";
	  unsigned int last_pushed = (*statistics.last_count_pushed)[i];
	  unsigned int last_popped = (*statistics.last_count_popped)[i];
	  // if (last_pushed > last_popped) std::cerr << last_pushed << " ";
	  std::cerr << last_pushed << " ";
	  if (i%8 == 7) std::cerr << "\n";
	}
#else
	/* The way I prefer (at the moment): a short indicator string */
	unsigned int last_seen = 1;
	char indicators[] = {'.', ';', '!', '|'};
	int indicator = 0;
	for (int i = 0; i < polys.size(); i ++) {
	  // std::cerr << (*statistics.last_count_pushed)[i] << " " << (*statistics.last_count_popped)[i] << " ";
	  unsigned int last_pushed = (*statistics.last_count_pushed)[i];
	  // if (last_pushed > last_popped) std::cerr << last_pushed << " ";
	  if (last_pushed == last_seen) std::cerr << indicators[indicator];
	  else if (last_pushed > last_seen) {
	    last_seen = last_pushed;
	    if (indicator + 1 < sizeof(indicators)) indicator ++;
	    std::cerr << indicators[indicator];
	  }
	}
#endif
	std::cerr << "\n";
    }
}

void run_some_basic_tests(void)
{
  // std::cout << sizeof(BitString::data_type) << "\n";
  BitString bs7("010000000");
  BitString bs8("100000000");
  BitString bs9("1000000000");
  BitString bs6364("11000000000000000000000000000000000000000000000000000000000000000");
  BitString bs64  ("10000000000000000000000000000000000000000000000000000000000000000");
  BitString bs65 ("100000000000000000000000000000000000000000000000000000000000000000");
  BitString bs3132("110000000000000000000000000000000");
  BitString bs32  ("100000000000000000000000000000000");
  BitString bs33 ("1000000000000000000000000000000000");
  // std::vector<int> bwb(bs6364.begin(), bs6364.end());
  assert(std::vector<int>(bs7.begin(), bs7.end()) == std::vector<int>({7}));
  assert(std::vector<int>(bs8.begin(), bs8.end()) == std::vector<int>({8}));
  assert(std::vector<int>(bs9.begin(), bs9.end()) == std::vector<int>({9}));
  assert(std::vector<int>(bs6364.begin(), bs6364.end()) == std::vector<int>({63,64}));
  assert(std::vector<int>(bs64.begin(), bs64.end()) == std::vector<int>({64}));
  assert(std::vector<int>(bs65.begin(), bs65.end()) == std::vector<int>({65}));
  assert(std::vector<int>(bs3132.begin(), bs3132.end()) == std::vector<int>({31,32}));
  assert(std::vector<int>(bs32.begin(), bs32.end()) == std::vector<int>({32}));
  assert(std::vector<int>(bs33.begin(), bs33.end()) == std::vector<int>({33}));
  assert(bs8.test_bit(8));
  assert(bs9.test_bit(9));
  assert(bs7.test_bit(7));
  assert(! bs9.test_bit(7));
  BitString bs(bs8);
  bs.clear_bit(8);
  assert(! bs.test_bit(8));
  bs.set_bit(7);
  bs.set_bit(0);
  assert(bs.test_bit(0));
  assert(std::vector<int>(bs.begin(), bs.end()) == std::vector<int>({0,7}));
}

int main(int argc, char ** argv)
{
  int nthreads = 1;
  int bitstring_len = 0;
  int opt;

  run_some_basic_tests();

  while ((opt = getopt(argc, argv, "vt:")) != -1) {
    switch (opt) {
    case 'v':
      verbose = true;
      break;
    case 't':
      nthreads = atoi(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-v] [-t nthreads]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  if (optind < argc) {
    fprintf(stderr, "Usage: %s [-v] [-t nthreads]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  std::string bitstring;
  while (std::getline(std::cin, bitstring)) {
    /* BitString's constructor can take a string */
    polys.emplace_back(bitstring);

    BitString bs(bitstring);
    if (bitstring_len == 0) {
      bitstring_len = bs.len;
      // std::cerr << bitstring_len << " byte bitstrings\n";
    } else if (bitstring_len != bs.len) {
      std::cerr << "bitstring lengths aren't consistent\n";
      exit(1);
    }

    inverted_polys.emplace_back(~bs);

    /* Put an empty list at the back of expanded_polys and fill it with the individual bits */
    expanded_polys.emplace_back();
    expanded_polys_bits.emplace_back();
    BitString rmsb;
    do {
      rmsb = bs.rightmost_set_bit();
      expanded_polys.back().emplace_back(rmsb);
      expanded_polys_bits.back().emplace_back(bs.rightmost_set_bit_index());
      bs ^= rmsb;
    } while (bs);
  }

  statistics.subset_checks_made = new std::vector<std::atomic<unsigned int>>(polys.size());
  statistics.subset_checks_succeeded = new std::vector<std::atomic<unsigned int>>(polys.size());
  statistics.last_count_pushed = new std::vector<std::atomic<unsigned int>>(polys.size());
  statistics.last_count_popped = new std::vector<std::atomic<unsigned int>>(polys.size());

  compute_and_remove_single_bit_covers();

  compute_all_covers();
  if (verbose) compute_and_display_statistics();

  if (all_covers.size() > 0) {
    BacktrackPoint initial_work;
    initial_work.next_polynomial = 0;
    initial_work.next_index = 0;
    initial_work.allowed_bits = ~ BitString(bitstring_len);

    for (auto& cover: all_covers) {
      int num_attachment_points = cover.single_link_chains.size();
      /* There's 2^num_attachment_points possible combinations of attachment point "values".
       * For each value, either the attachment point is included and all of the outliers are excluded, or vice versa.
       * Loop through all of them, creating a BacktrackPoint for each and adding it to the work queue.
       */
      for (int attachment_points_value = 0; (attachment_points_value & (1 << num_attachment_points)) == 0; attachment_points_value ++) {
	int working_value = attachment_points_value;
	initial_work.bitstring = BitString(bitstring_len);
	initial_work.allowed_bits = ~ BitString(bitstring_len);
	initial_work.cover = &cover;
	for (auto const &[attachment_point, value] : cover.single_link_chains) {
	  if (working_value & 1) {
	    initial_work.bitstring |= attachment_point;
	    initial_work.allowed_bits &= ~ value;
	  } else {
	    initial_work.bitstring |= value;
	    initial_work.allowed_bits &= ~ attachment_point;
	  }
	  working_value >>= 1;
	}
	backtrack_queue.push(initial_work);
      }
    }

    // Set up the signal handler
    struct sigaction sa;
    sa.sa_handler = signalHandler; // Use custom signal handler
    sa.sa_flags = 0;               // Use default flags
    sigemptyset(&sa.sa_mask);      // No signals blocked during handler
    if (sigaction(SIGQUIT, &sa, nullptr) == -1) {
        perror("sigaction");
        return 1;
    }

    backtrack_queue.set_num_workers(nthreads);
    std::vector<std::thread> threads;

    for (int i=0; i < nthreads; i++) {
      threads.emplace_back(task);
    }

    for (int i=0; i < nthreads; i++) {
      threads[i].join();
    }

    /* We've now got the single_bit_covers, which have to appear in all output product terms,
     * and a list of covers with their finished bitstrings.  We now want to form the product
     * of the finished bitstrings, every possible combination of one finished bitstring from
     * each cover.  To achieve this, we form a vector of iterators over all of the covers
     * and advance the first across its entire vector of bitstrings, then reset the first
     * and advance the second, and so on until we've iterated over all combinations.
     *
     * This code assumes that all of the covers have at least one finished_bitstring.
     * Since a vector of all one bits would satisify all possible input conditions,
     * this seems like a valid assumption (baring a bug).  Since I've hit such a bug,
     * we now check for this condition and throw a runtime error if it happens.
     */

    std::vector<FinishedBitStrings *> finished_bitstrings;
    std::vector<FinishedBitStrings::iterator> iterators;
    for (auto &cover: all_covers) {
      finished_bitstrings.push_back(&cover.finished_bitstrings);
      if (cover.finished_bitstrings.begin() == cover.finished_bitstrings.end()) {
	throw std::runtime_error("Cover with no finished bitstrings");
      }
      iterators.push_back(cover.finished_bitstrings.begin());
    }
    while (true) {
      BitString result = single_bit_covers;
      for (auto &it: iterators) {
	result |= *it;
      }
      std::cout << result << "\n";
      int i;
      for (i=0; i<iterators.size(); i++) {
	iterators[i] ++;
	if (iterators[i] != finished_bitstrings[i]->end()) {
	  break;
	} else {
	  iterators[i] = finished_bitstrings[i]->begin();
	}
      }
      if (i == iterators.size()) break;
    }
  } else {
    std::cout << single_bit_covers << "\n";
  }

  return 0;
}
