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
 * Start with a bit string of all zeros (an empty ideal) and
 * and an empty list of finished bit strings, but note that we might remove a bit string from the list of
 * finished bit strings, so they're not completely "finished".  (Or maybe not, as we don't seem to be
 * removing anything from the finished bit strings)
 *
 * We track the number of the polynomial we're working on.  First check to see if there's
 * any overlap between the polynomial's bit string and the working bit string (bitwise AND).
 * If so, this polynomial is satisfied; move on to the next one.
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
 *    - working ideal bit string
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
#include <list>
#include <map>
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

  /* Bitwise equality */

  bool operator==(const BitString& rhs) const
  {
    for (size_type i=0; i<bitstring.size(); i++) {
      if (bitstring[i] != rhs.bitstring[i])
	return false;
    }
    return true;
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

  int count(void) const
  {
    int count = 0;
    for (size_type i=0; i<bitstring.size(); i++) {
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

struct BacktrackPoint
{
  BitString bitstring;
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

class Cover
{
  public:
  BitString cover;
  std::list<int> links;
  std::map<BitString, BitString> single_link_chains;
  bool all_chains_are_single_links;
  int triplets;
};

std::list<Cover> all_covers;

FinishedBitStrings finished_bitstrings;

/* Which polynomials are to be considered in the calculations (initially all of them),
 * but some will be removed for optimization.
 */
std::vector<bool> under_consideration;

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
  BacktrackPoint extra_work;

  while (true) {
    get_next_work_from_queue:

    /* waitAndPop returns true if we've got work; false if all the work is done */
    if (! backtrack_queue.waitAndPop(current_work)) return;

    /* I put this here because while the work was on the queue, we could have added new finished bitstrings */
    /* There's probably a race condition here, as a new finished bitstring could appear at any time! */
    if (finished_bitstrings.contain_a_subset_of(current_work.bitstring)) {
      // std::cerr << "detected superset\n";
      continue;
    }

    while (current_work.next_polynomial < polys.size()) {
      /* If there's any overlap here, the polynomial is already satisfied, move on */
      // std::cerr << "working on " << current_work.bitstring << "\n";
      if (under_consideration[current_work.next_polynomial] && ! (current_work.bitstring && polys[current_work.next_polynomial])) {
	bool have_next_work = false;
	/* Extend backtrack queue if there's more than one factor(bit) in the polynomial */
	/* XXX should this loop use a reference for speed? */
	for (BitString next_bit: expanded_polys[current_work.next_polynomial]) {

	  // std::cerr << "adding " << next_bit << "\n";
	  if (! have_next_work) {
	    /* Like this, but avoids a copy */
	    /* next_work_bitstring = current_work.bitstring | next_bit; */
	    current_work.bitstring.logical_or_assign(next_work_bitstring, next_bit);
	    /* Check first if this is a superset of an existing bit string; skip it if it is */
	    if (finished_bitstrings.contain_a_subset_of(next_work_bitstring)) continue;
	    have_next_work = true;
	  } else {
	    /* Like this, but avoids a copy */
	    /* extra_work.bitstring = current_work.bitstring | next_bit; */
	    current_work.bitstring.logical_or_assign(extra_work.bitstring, next_bit);
	    extra_work.next_polynomial = current_work.next_polynomial + 1;
	    /* Check first if this is a superset of an existing bit string; skip it if it is */
	    if (finished_bitstrings.contain_a_subset_of(extra_work.bitstring)) continue;
	    backtrack_queue.push(extra_work);
	  }
	}
	if (! have_next_work) {
	  goto get_next_work_from_queue;
	}
	current_work.bitstring = next_work_bitstring;
      }
      current_work.next_polynomial ++;
    }
    finished_bitstrings.add(current_work.bitstring);
  }
}

/* Remove from consideration all polynomials with only a single bit set, and all others that depend on them.
 *
 * Modifies global variables 'single_bit_covers' and 'under_consideration'.
 */

void compute_and_remove_single_bit_covers(void)
{
  single_bit_covers = BitString(polys[0].len);

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
 * can be handled with lookup tables (not yet implemented).
 *
 * Given a cover and a set of polynomials "under consideration", expand the
 * cover to include any polynomials under consideration that partially
 * match the cover.  The cover argument is modified.  Once we're done,
 * all of the polynomials "under consideration" are either completely
 * under the (expanded) cover, or are not under it at all.
 * Return the number of polynomials under the expanded cover.
 */
int expand_cover(BitString& cover, std::vector<bool> under_consideration)
{
    bool expanding_cover;
    int polys_covered = 0;
    do {
      expanding_cover = false;
      for (auto i=0; i<polys.size(); i++) {
	if (under_consideration[i] && (cover && polys[i])) {
	  under_consideration[i] = false;
	  cover |= polys[i];
	  polys_covered ++;
	  expanding_cover = true;
	}
      }
    } while (expanding_cover);
    return polys_covered;
}

bool is_link(BitString top_cover, int link, std::vector<bool> under_consideration)
{
  /* under_consideration is a local copy of the global under_consideration */

  if (polys[link].count() != 2) {
    /* Links, by definition, have only two bits set */
    return false;
  }

  /* only consider polynomials under the top cover */
  for (int i=0; i<polys.size(); i++) {
    if (under_consideration[i]) {
      under_consideration[i] = polys[i] && top_cover;
    }
  }
  /* Remove link from consideration */
  under_consideration[link] = false;

  /* Does this cause the top cover to partition in two? */

  BitString cover;
  int i;
  for (i=0; i<polys.size(); i++) {
    if (under_consideration[i]) {
      cover = polys[i];
      break;
    }
  }
  if (i == polys.size()) {
    /* If there's only one polynomial under the cover, then it's a link */
    return true;
  }

  expand_cover(cover, under_consideration);

  /* If the expanded cover equals the top cover, then it's not a link.  Otherwise,
   * the expanded cover is smaller because the top_cover partitioned into two,
   * and it is a link.
   */
  return ! (cover == top_cover);
}

void compute_all_covers(void)
{
  /* For each remaining polynomial, form a cover that is initially just the first polynomial.
   * Then loop over all the polynomials, OR-ing in any polynomials that match the cover.  Keep looping
   * until the cover stabilizes.  Remove all of these polynomials from consideration, and keep doing
   * it until we've completely partitioned the input set.
   */
  BitString union_of_all_covers(polys[0].len);

  while (true) {
    Cover cover;
    int polys_covered;
    int i;
    for (i=0; i<polys.size(); i++) {
      if (under_consideration[i] && ! (union_of_all_covers && polys[i])) {
	cover.cover = polys[i];
	break;
      }
    }
    if (i == polys.size()) {
      break;
    }

    polys_covered = expand_cover(cover.cover, under_consideration);
    union_of_all_covers |= cover.cover;

    /* Identify which polynomials are links (2-bit polynomials that partition the cover when removed) */
    for (i=0; i<polys.size(); i++) {
      if (under_consideration[i] && (polys[i] && cover.cover)) {
	if (is_link(cover.cover, i, under_consideration)) {
	  cover.links.emplace_back(i);
	}
      }
    }

    /* Form the links into chains.  Remove the first link from the unchained_links list, then remove
     * all other links that it can chain with, tracking the left and right ends of the chain.
     * Keep doing this until all links have been formed into chains.
     */
    std::list<int> unchained_links = cover.links;
    cover.all_chains_are_single_links = true;
    while (unchained_links.size() > 0) {
      int left_of_chain = unchained_links.front();
      int right_of_chain = unchained_links.front();
      unchained_links.pop_front();
      int length_of_chain = 1;
      std::list<int>::iterator it;
      bool attached_a_link;
      // std::cerr << "start of chain " << left_of_chain << " " << polys[left_of_chain] << " " << polys[left_of_chain].rightmost_set_bit() << "\n";
      do {
	attached_a_link = false;
	for (it = unchained_links.begin(); it != unchained_links.end(); it++) {
	  if (polys[left_of_chain] && polys[*it].rightmost_set_bit()) {
	    left_of_chain = *it;
	    unchained_links.erase(it);
	    length_of_chain ++;
	    attached_a_link = true;
	    // std::cerr << "attached at left " << left_of_chain << " " << polys[left_of_chain] << "\n";
	    break;
	  }
	  if (polys[right_of_chain].rightmost_set_bit() && polys[*it]) {
	    right_of_chain = *it;
	    unchained_links.erase(it);
	    length_of_chain ++;
	    attached_a_link = true;
	    // std::cerr << "attached at right " << right_of_chain << " " << polys[right_of_chain] << "\n";
	    break;
	  }
	}
      } while (attached_a_link);
      /* If the length of chain is one, there are several possibilities: either it's a "dangling link"
       * (the case we want to optimize), or the link is itself the entire cover, or it's a link
       * between two subcovers.  The link has two sides (left and right).  Count up how many polynomials
       * match up with each.  A dangling link has one side matching a single polynomial (the link itself)
       * and the other sides matches more than one.  A link between two subcovers has both sides
       * matching more than one polynomial.  If both sides match a single polynomial, then the
       * link is the entire cover.
       */
      if (length_of_chain == 1) {
	BitString link_right_bit = polys[left_of_chain].rightmost_set_bit();
	BitString link_left_bit = polys[left_of_chain] & ~link_right_bit;
	int right_bit_polys = 0;
	int left_bit_polys = 0;
	BitString attachment_point;
	BitString outlying_point;
	for (i=0; i<polys.size(); i++) {
	  if (under_consideration[i] && (polys[i] && link_right_bit)) right_bit_polys ++;
	  if (under_consideration[i] && (polys[i] && link_left_bit)) left_bit_polys ++;
	}
	if ((right_bit_polys == 1) && (left_bit_polys > 1)) {
	  attachment_point = link_left_bit;
	  outlying_point = link_right_bit;
	} else if ((right_bit_polys > 1) && (left_bit_polys == 1)) {
	  attachment_point = link_right_bit;
	  outlying_point = link_left_bit;
	} else {
	  /* not a case we handle */
	  cover.all_chains_are_single_links = false;
	}
#if 0
	if (attachment_point.count() > 0) {
	  if (! cover.single_link_chains.contains(attachment_point)) {
	    cover.single_link_chains[attachment_point] = BitString(polys[0].len);
	  }
	  cover.single_link_chains[attachment_point] |= outlying_point;
	}
#endif
      } else {
	cover.all_chains_are_single_links = false;
      }
    }

    for (i=0; i<polys.size(); i++) {
      if (under_consideration[i] && (polys[i] && cover.cover) && (polys[i].count() == 2)) {
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
      if (under_consideration[i] && (polys[i] && cover.cover) && (polys[i].count() == 3)) {
	int single_poly_bits = 0;
	for (BitString next_bit: expanded_polys[i]) {
	  int matching_polys = 0;
	  for (auto j=0; j<polys.size(); j++) {
	    if (under_consideration[j] && (polys[j] && next_bit)) matching_polys ++;
	  }
	  if (matching_polys == 1) single_poly_bits ++;
	}
	if (single_poly_bits == 2) cover.triplets ++;
      }
    }

    all_covers.push_back(cover);
  }

}

/* Compute and display some statistics about the input data */

void compute_and_display_statistics(void)
{
  /* How many covers are there of each size, and how many polynomials are covered */
  std::vector<int> covers(polys[0].len+1, 0);
  std::vector<int> polynomials_covered(polys[0].len+1, 0);
  std::vector<int> count_of_links(polys[0].len+1, 0);
  std::vector<int> count_of_chains(polys[0].len+1, 0);
  std::vector<int> length_of_chains(polys[0].len+1, 0);

  /* Remove all links from consideration, because we're going to repeat the cover calculation later without them. */
  /* for (auto link: links) under_consideration[link] = false; */

  if (single_bit_covers.count() == 1) std::cerr << "1 single bit cover covering ";
  else std::cerr << single_bit_covers.count() << " single bit covers covering ";
  if (single_bit_polynomials_covered == 1) std::cerr << "1 polynomial\n";
  else std::cerr << single_bit_polynomials_covered << " polynomials\n";

  for (auto cover: all_covers) {
    std::cerr << "a " << cover.cover.count() << "-bit cover with " << cover.links.size() << " links";
    if (cover.all_chains_are_single_links) {
      std::cerr << " (all single links) with " << cover.single_link_chains.size() << " attachment points";
    }
    if (cover.triplets == 1) std::cerr << "; 1 triplet\n";
    else std::cerr << "; " << cover.triplets << " triplets \n";
    /* missing: print number of polynomials covered and max length of chain (if not 1) */
  }

#if 0
  /* repeat the cover calculation without the links */

  all_covers = BitString(polys[0].len);
  covers = std::vector<int>(polys[0].len+1, 0);
  polynomials_covered = std::vector<int>(polys[0].len+1, 0);
  while (true) {
    BitString cover;
    int polys_covered;
    int i;
    for (i=0; i<polys.size(); i++) {
      if (under_consideration[i] && ! (all_covers && polys[i])) {
	cover = polys[i];
	break;
      }
    }
    if (i == polys.size()) {
      break;
    }

    polys_covered = expand_cover(cover, under_consideration);
    all_covers |= cover;

    covers[cover.count()] ++;
    polynomials_covered[cover.count()] += polys_covered;
  }
  std::cerr << "After links removed:\n";
  for (int i = 1; i <= polys[0].len; i ++) {
    if (covers[i] == 1) {
      std::cerr << "1 " << i << "-bit cover covering " << polynomials_covered[i] << " polynomials\n";
    } else if (covers[i] > 1) {
      std::cerr << covers[i] << " " << i << "-bit covers covering " << polynomials_covered[i] << " polynomials\n";
    }
  }
#endif
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

  under_consideration.resize(polys.size(), true);
  compute_and_remove_single_bit_covers();

  compute_all_covers();
  compute_and_display_statistics();

  /* If we can optimize (a single cover with all single links), remove the links from consideration */
  bool optimize = (all_covers.size() == 1 && all_covers.front().all_chains_are_single_links);
  optimize = false;
  if (optimize) {
    for (auto link: all_covers.front().links) {
      under_consideration[link] = false;
    }
  }

  BacktrackPoint initial_work;
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
    if (p.len == 0) {
      if (optimize) {
	/* an optimized empty cover with single link chains must be a two bit cover composed of a single link */
	for (auto const &[key, value] : all_covers.front().single_link_chains) {
	  std::cout << key << "\n";
	  std::cout << value << "\n";
	}
      } else {
	std::cout << single_bit_covers << "\n";
      }
    } else {
      p |= single_bit_covers;
      if (optimize) {
	for (auto const &[key, value] : all_covers.front().single_link_chains) {
	  if (! (p && key)) p|= value;
	}
      }
      std::cout << p << "\n";
    }
  }

  return 0;
}
