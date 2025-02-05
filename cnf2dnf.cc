/*
 * cnf2dnf - conjunctive normal form (CNF) to disjunctive normal form (DNF) converter
 *
 * BUILD: g++ -std=c++2a -march=native -O3 -o cnf2dnf cnf2dnf.cc -lpthread
 *
 * USAGE: ./cnf2dnf NUM_THREADS < input-bit-strings > output-bit-strings
 *
 * GOAL: Given a set of polynomials (that form an ideal), we want to factor them all and
 * form a set of ideals, all formed from irreducible polynomials.  The essence of
 * this operation is converting a logical expression from CNF to DNF.
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
#include <list>
#include <map>
#include <set>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include "LockingQueue.hpp"

// unsigned int bitstring_len = 0;

bool verbose = false;

class BitString
{
public:
  typedef unsigned long long int data_type;
  typedef std::vector<data_type>::size_type size_type;
  unsigned int len;
  std::vector<data_type> bitstring;

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
    if (len != rhs.len) throw std::runtime_error("incompatible lengths in operator|=");
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] |= rhs.bitstring[i];
    }
    return *this;
  }

  BitString operator&=(const BitString& rhs)
  {
    if (len != rhs.len) throw std::runtime_error("incompatible lengths in operator&=");
    for (size_type i=0; i<bitstring.size(); i++) {
      bitstring[i] &= rhs.bitstring[i];
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
    for (size_type i=bitstring.size()-1; i>=0; i--) {
      data_type rightmost_set_bit = bitstring[i] & (-bitstring[i]);
      if (rightmost_set_bit) {
	result.bitstring[i] = rightmost_set_bit;
	return result;
      }
    }
    return result;
  }

  int rightmost_set_bit_index(void) const
  {
    for (int i=bitstring.size()-1; i>=0; i--) {
      for (auto j=0; j<8*sizeof(data_type); j++) {
	if (bitstring[i] & (1<<j)) {
	  return i*sizeof(data_type) + j;
	}
      }
    }
    return -1;
  }

  int count(void) const
  {
    int count = 0;
    for (int i=0; i<bitstring.size(); i++) {
      count += std::popcount(bitstring[i]);
    }
    return count;
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

class CompareByCount {
  public:
  bool operator()(const BitString &lhs, const BitString &rhs) const
  {
    return (lhs.count() < rhs.count());
  }
};

typedef std::multiset<BitString, CompareByCount> FinishedBitStrings2;

class FinishedBitStrings
{
public:
  FinishedBitStrings2 finished_bitstrings;
  std::shared_mutex mutex;

  bool contain_a_superset_of(const BitString& bitstring)
  {
    std::shared_lock<std::shared_mutex> lock(mutex);

    for (auto it = finished_bitstrings.rbegin(); it != finished_bitstrings.rend(); it ++) {
      if (bitstring.count() > it->count())
	return false;
      if (it->is_superset_of(bitstring))
	return true;
    }
    return false;
  }

  bool contain_a_subset_of(const BitString& bitstring)
  {
    std::shared_lock<std::shared_mutex> lock(mutex);

    /* If bitstring has fewer bits than anything in finished_bitstrings, then
     * there's no way that anything in finished_bitstrings is a subset of bitstring.
     */

    for (const BitString& fbs: finished_bitstrings) {
      if (bitstring.count() < fbs.count())
	return false;
      if (bitstring.is_superset_of(fbs))
	return true;
    }
    return false;
  }

  void add(const BitString &bitstring)
  {
    int bitstring_count = bitstring.count();
    std::unique_lock<std::shared_mutex> lock(mutex);

    /* A final extra check for subsets, to avoid problems from the race condition mentioned below.
     * We don't call contain_a_subset_of because it locks the object.
     */

    for (const BitString& fbs: finished_bitstrings) {
      if (bitstring.count() < fbs.count())
	break;
      if (bitstring.is_superset_of(fbs))
	return;
    }

    // std::cerr << "add " << bitstring << "\n";
    /* remove any supersets */
    std::erase_if(finished_bitstrings, [&](const BitString& fbs) { return fbs.is_superset_of(bitstring); });

    finished_bitstrings.insert(bitstring);
  }
};

class Cover
{
  public:
  BitString cover;                                       /* which bits form the cover */
  std::vector<bool> under_consideration;                 /* which polynomials are candidates */
  std::map<BitString, BitString> single_link_chains;
  int triplets;
  FinishedBitStrings finished_bitstrings;
};

struct BacktrackPoint
{
  BitString bitstring;
  BitString allowed_bits;
  unsigned int next_polynomial;
  FinishedBitStrings * finished_bitstrings;
  Cover * cover;
};

std::vector<BitString> polys;
std::vector<BitString> inverted_polys;
std::vector<std::vector<BitString>> expanded_polys;

LockingQueue<BacktrackPoint> backtrack_queue;

std::list<Cover> all_covers;

/* All bits for which there is a polynomial containing only that bit.  Such bits must
 * always be present in the output bit strings, and any polynomials that depend on them
 * can be removed from under_consideration, since they will always be satisfied.
 */
BitString single_bit_covers;
int single_bit_polynomials_covered = 0;

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

    /* I used to have a test here: finished_bitstrings.contain_a_subset_of(current_work.bitstring),
     * because even though we tested for this before adding the work to the backtrack_queue,
     * we could have added new finished bitstrings while the work was on the queue.  The program
     * spent a lot of time here, and I decided that we don't really need this check because
     * we're going to check contain_a_subset_of if we either add work to the work queue,
     * or add a bitstring to finished_bitstrings.
     */

    while (current_work.next_polynomial < polys.size()) {
      /* If there's any overlap here, the polynomial is already satisfied, move on */
      // std::cerr << "working on " << current_work.bitstring << "\n";
      if (current_work.cover->under_consideration[current_work.next_polynomial] && ! (current_work.bitstring && polys[current_work.next_polynomial])) {
	bool have_next_work = false;
	/* Extend backtrack queue if there's more than one factor(bit) in the polynomial */
	/* XXX should this loop use a reference for speed? */
	for (BitString next_bit: expanded_polys[current_work.next_polynomial]) {
	  /* skip this bit if it's not in the allowed bit set */
	  if (next_bit && current_work.allowed_bits) {
	    // std::cerr << "adding " << next_bit << "\n";
	    if (! have_next_work) {
	      /* Like this, but avoids a copy */
	      /* next_work_bitstring = current_work.bitstring | next_bit; */
	      current_work.bitstring.logical_or_assign(next_work_bitstring, next_bit);
	      /* remove all of the current polynomial's bits from the allowed bits for future work */
	      next_work_allowed_bits = current_work.allowed_bits;
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      if (current_work.finished_bitstrings->contain_a_subset_of(next_work_bitstring)) continue;
	      have_next_work = true;
	    } else {
	      /* Like this, but avoids a copy */
	      /* extra_work.bitstring = current_work.bitstring | next_bit; */
	      current_work.bitstring.logical_or_assign(extra_work.bitstring, next_bit);
	      extra_work.allowed_bits = current_work.allowed_bits;
	      extra_work.next_polynomial = current_work.next_polynomial + 1;
	      extra_work.finished_bitstrings = current_work.finished_bitstrings;
	      extra_work.cover = current_work.cover;
	      /* Check first if this is a superset of an existing bit string; skip it if it is */
	      if (current_work.finished_bitstrings->contain_a_subset_of(extra_work.bitstring)) continue;
	      backtrack_queue.push(extra_work);
	    }
	    /* If the current polynomial's bits are 111, we want to create future work 1xx, 01x, 001,
	     * (not 1xx, x1x, xx1), so we now remove the current bit from the allowed bits.
	     */
	    current_work.allowed_bits &= ~next_bit;
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
    current_work.finished_bitstrings->add(current_work.bitstring);
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

/* Partition the input bits into sets that completely cover a group of polynomials,
 * and display the number of bits in the set and the number of polynomials covered.
 *
 * Smaller partitions are easier to handle: 1 is a special case, and 2 through 5
 * can be handled with lookup tables (not implemented, because they don't appear in practice).
 *
 * Given a cover and a set of polynomials "under consideration", expand the
 * cover to include any polynomials under consideration that partially
 * match the cover.  The cover argument is modified.  Once we're done,
 * all of the polynomials "under consideration" are either completely
 * under the (expanded) cover, or are not under it at all.
 * Return the number of polynomials under the expanded cover.
 */
int expand_cover(Cover& cover)
{
    bool expanding_cover;
    int polys_covered = 0;
    do {
      expanding_cover = false;
      for (auto i=0; i<polys.size(); i++) {
	if (! cover.under_consideration[i] && ! (single_bit_covers && polys[i]) && (cover.cover && polys[i])) {
	  cover.under_consideration[i] = true;
	  cover.cover |= polys[i];
	  polys_covered ++;
	  expanding_cover = true;
	}
      }
    } while (expanding_cover);
    return polys_covered;
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

    Cover& cover = all_covers.emplace_back();
    cover.cover = polys[i];
    cover.under_consideration.resize(polys.size(), false);

    polys_covered = expand_cover(cover);
    union_of_all_covers |= cover.cover;

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

    for (i=0; i<polys.size(); i++) {
      if (cover.under_consideration[i] && (polys[i] && cover.cover) && (polys[i].count() == 2)) {
	BitString outlying_point;
	BitString attachment_point;
	for (BitString current_bit: expanded_polys[i]) {
	  int matching_polys = 0;
	  for (auto j=0; j<polys.size(); j++) {
	    if (cover.under_consideration[j] && (polys[j] && current_bit)) matching_polys ++;
	  }
	  if (matching_polys == 1) {
	    outlying_point = current_bit;
	  } else {
	    attachment_point = current_bit;
	  }
	}
	if ((attachment_point.count() > 0) && (outlying_point.count() > 0)) {
	  if (! cover.single_link_chains.contains(attachment_point)) {
	    cover.single_link_chains[attachment_point] = BitString(polys[0].len);
	  }
	  cover.single_link_chains[attachment_point] |= outlying_point;
	}
      }
    }

    /* Count "triplets": Polynomials with two isolated bits and only
     * one touching the rest of the cover.
     */

    cover.triplets = 0;
    for (i=0; i<polys.size(); i++) {
      if (cover.under_consideration[i] && (polys[i] && cover.cover) && (polys[i].count() == 3)) {
	int single_poly_bits = 0;
	for (BitString next_bit: expanded_polys[i]) {
	  int matching_polys = 0;
	  for (auto j=0; j<polys.size(); j++) {
	    if (cover.under_consideration[j] && (polys[j] && next_bit)) matching_polys ++;
	  }
	  if (matching_polys == 1) single_poly_bits ++;
	}
	if (single_poly_bits == 2) cover.triplets ++;
      }
    }
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
      // std::cerr << bitstring_len << " byte bitstrings\n";
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

  compute_and_remove_single_bit_covers();

  compute_all_covers();
  if (verbose) compute_and_display_statistics();

  if (all_covers.size() > 0) {
    BacktrackPoint initial_work;
    initial_work.next_polynomial = 0;
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
	initial_work.finished_bitstrings = & cover.finished_bitstrings;
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
     * We form a std::vector of pointers to the finished_bitstrings because the
     * finished_bitstrings variable in Cover is a std::list, which can't be indexed
     * like a std::vector.  We can't change the std::list in Cover to a std::vector,
     * because std::vector's emplace_back method requires its members to be MoveInsertable,
     * and FinishedBitstrings includes a std::mutex, which can't be moved because
     * threads could be waiting on it.
     *
     * This code assumes that all of the covers have at least one finished_bitstring.
     * Since a vector of all one bits would satisify all possible input conditions,
     * this seems like a valid assumption (baring a bug).
     */

    std::vector<FinishedBitStrings2 *> finished_bitstrings;
    std::vector<FinishedBitStrings2::iterator> iterators;
    for (auto &cover: all_covers) {
      finished_bitstrings.push_back(& cover.finished_bitstrings.finished_bitstrings);
      iterators.push_back(cover.finished_bitstrings.finished_bitstrings.begin());
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
