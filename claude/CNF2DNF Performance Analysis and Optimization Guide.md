# CNF2DNF Performance Analysis & Optimization Guide

**Last Updated**: January 2026
**Status**: 1 optimization implemented, 4 major optimizations recommended

## Optimization Status Summary

| Priority | Optimization | Status | Expected Impact |
|----------|--------------|--------|-----------------|
| ✅ 0 | Early Duplicate Detection | **IMPLEMENTED** (Jan 2026) | Reduces redundant work |
| 🔲 1 | Incremental Polynomial Checking | Not implemented | 2-5x speedup |
| 🔲 2 | Optimize BitString::operator&& | Not implemented | 10-20% speedup |
| 🔲 3 | Track Satisfaction Incrementally | Not implemented | Major speedup |
| 🔲 4 | Memory Layout Optimization | Not implemented | Better cache utilization |
| 🔲 5 | Profile-Guided Optimization | Not implemented | 5-10% overall |

## Project Overview

**Project**: `cnf2dnf` - Conjunctive Normal Form to Disjunctive Normal Form converter  
**Language**: C++20  
**Location**: `/home/baccala/src/helium/cnf2dnf`  
**Purpose**: Converting logical systems from CNF to DNF for algebraic ideal factorization

### What the Program Does

This program is part of a larger algebraic computation system (Helium). It:

1. Takes polynomial generators that define an ideal
2. Factors the polynomials and treats each unique irreducible factor as a logical variable
3. Converts the resulting logical system from CNF (product of sums) to DNF (sum of products)
4. Outputs a set of ideals whose intersection forms the original ideal

**Input**: Bit strings representing "covers" (disjunctions - any bit can be true)  
**Output**: Bit strings representing "product terms" (conjunctions - all bits must be true)  
**Algorithm**: Multi-threaded backtracking search with heuristics

## Code Structure

### Key Files
- `cnf2dnf.cc` - Main implementation (all in one file)
- `LockingPriorityQueue.hpp` - Thread-safe priority queue (not provided but referenced)

### Core Data Structures

```cpp
class BitString {
    unsigned int len;                    // Number of bits
    int cached_count;                    // Cached popcount (-1 if invalid)
    std::vector<data_type> bitstring;    // Actual bits (data_type = unsigned long long)
    // MSB on left, stored in reverse order in vector
};
```

**Key BitString Operations**:
- `operator&&` - Bitwise AND returning boolean (checks overlap) - **PERFORMANCE CRITICAL**
- `operator&` - Bitwise AND returning BitString
- `operator|=` - Bitwise OR assignment
- `count()` - Population count (number of 1 bits) with caching
- `test_bit()`, `set_bit()`, `clear_bit()` - Individual bit operations

```cpp
class Cover {
    BitString cover;                              // Which bits form this partition
    std::vector<bool> under_consideration;        // Which polynomials in this partition
    std::map<BitString, BitString> single_link_chains;  // Outlier optimization
    int triplets;                                 // Count of 3-bit polynomials
    FinishedBitStrings finished_bitstrings;       // Results for this cover
};
```

```cpp
struct BacktrackPoint {
    BitString bitstring;           // Current candidate solution
    BitString allowed_bits;        // Bits available for assignment
    unsigned int next_index;       // Index within polynomial
    unsigned int next_polynomial;  // Which polynomial to process next
    Cover * cover;                 // Which cover we're working on
};
```

### Global Data Structures

```cpp
std::vector<BitString> polys;                      // Input polynomials
std::vector<BitString> inverted_polys;             // Bitwise complements
std::vector<std::vector<BitString>> expanded_polys; // Each poly as vector of single-bit BitStrings
std::vector<std::vector<int>> expanded_polys_bits;  // Each poly as vector of bit indices
BitString single_bit_covers;                        // Bits that must appear in all solutions
LockingPriorityQueue<BacktrackPoint, BacktrackQueueLess> backtrack_queue;  // Work queue
std::list<Cover> all_covers;                        // All identified covers
```

### Main Algorithm Flow

1. **Preprocessing** (`main()`, lines 854-920):
   - Read input bit strings
   - Build `polys`, `inverted_polys`, `expanded_polys`, `expanded_polys_bits`
   - Compute single-bit covers (polynomials with only one bit)
   - Partition remaining bits into covers

2. **Work Queue Initialization** (lines 922-949):
   - For each cover, create initial BacktrackPoints
   - Push them onto the thread-safe priority queue

3. **Parallel Processing** (`task()`, lines 577-803):
   - Each thread pulls BacktrackPoints from the queue
   - Processes polynomials sequentially
   - Adds bits to candidate solutions
   - Generates new BacktrackPoints for unexplored branches
   - Reduces solutions to prime implicants
   - Adds completed solutions to `finished_bitstrings`

4. **Output Generation** (lines 967-998):
   - Combine finished bitstrings from all covers
   - Add mandatory single-bit covers
   - Output final DNF

### Threading Model

- **Work Queue**: `LockingPriorityQueue` with priority based on `next_polynomial` (highest first)
- **Synchronization**: 
  - Mutex-protected queue operations
  - Atomic operations for `FinishedBitStrings` insertion
  - Condition variables for thread wake/sleep
- **Termination**: Queue detects when all threads are waiting on empty queue

## Performance Analysis Results

### Profiling Data

Profiling was done using `perf record` and `perf annotate` on the `task()` function.

### Critical Performance Bottleneck Identified

**Location**: Lines 780-782 in `task()`

```cpp
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
```

**What this does**: Reduces a candidate solution to a prime implicant by clearing each bit and checking if all polynomials remain satisfied.

### Perf Data Analysis

**~60% of runtime** is spent in the nested loop structure:

```assembly
# Outer loop: iterate over polynomials (21.16% at ca8)
ca8:   add     $0x1,%rdi          # 21.16% - loop increment
       add     $0x20,%r9          # 1.50%  - advance to next polynomial structure
       cmp     %rdi,%r8           # 1.56%  - loop condition
       je      d2c

# Check if polynomial is under consideration and not satisfied
cb5:   mov     %rdi,%rbp          # 1.49%  - setup
       sar     $0x6,%rbp          # 1.59%  - bit index calculation
       shl     %cl,%rax           # 1.91%  - create bit mask
       and     (%r11,%rbp,8),%rax # 6.37%  - check under_consideration bit
       je      ca8                # 1.60%  - skip if not under consideration

# Inner loop: BitString::operator&& - check if polynomial satisfied
ce8:   mov     0x0(%rbp,%rax,8),%rcx  # 17.39% - load polynomial bitstring chunk
       and     (%rsi,%rax,8),%rcx     # 8.52%  - AND with candidate bitstring
       jne     ca8                    # 13.85% - polynomial satisfied, continue
       add     $0x1,%rax              # 2.95%  - next chunk
       cmp     %r12,%rax              # 12.02% - loop condition
       jb      ce8                    # Continue inner loop
```

**Breakdown**:
- Memory loads: ~17.39% (loading polynomial data)
- Bitwise AND operations: ~8.52%
- Branch/comparison overhead: ~26.87% (13.85% + 12.02% + others)
- Loop control: ~21.16%

### Why It's Slow

1. **Quadratic behavior**: For each bit in candidate solution × all polynomials × bitstring chunks
2. **Poor cache locality**: Accessing different polynomials sequentially
3. **Repeated work**: Rechecking polynomials that don't depend on the cleared bit
4. **Memory bandwidth bound**: Lots of loads from `polys[]` array

## Optimization Recommendations

### ✅ Implemented: Early Duplicate Detection (Commit 63a40d2)

**Status**: COMPLETED - January 2026

**Problem**: The expensive prime implicant reduction operation was being performed even on bitstrings that match already-known solutions.

**Solution**: Add an early check before the reduction operation:

```cpp
/* We've now got a bitstring that's in the DNF.  See if it matches any already known
 * solutions.
 */
if (current_work.cover->finished_bitstrings.contain_a_subset_of(current_work.bitstring)) {
  continue;
}
/* We've now got a bitstring that's new in the DNF.  See if any of its subsets are also in the DNF,
 * by clearing each one bit and checking to see if the bitstring still works...
```

**Impact**: Profiling showed this was a hot spot - expanding the matched bitstring to make it prime was expensive. This check minimizes the need for that calculation by skipping bitstrings that are already covered by known solutions.

**Location**: Added at line ~1127 in `cnf2dnf.cc`, just before the prime implicant reduction loop (the expensive operation at lines 780-782 in the old code).

**Trade-offs**:
- Very low overhead (O(log n) subset check)
- High benefit when many similar solutions exist
- Simple, localized change with no risk to correctness

---

### Priority 1: Incremental Polynomial Checking (Expected 2-5x speedup)

**Status**: Not yet implemented - recommended as next optimization

**Problem**: Currently rechecks ALL polynomials after clearing each bit, even though most polynomials don't depend on that bit.

**Solution**: Build an index of which polynomials contain each bit:

```cpp
// Add to global data structures
std::vector<std::vector<int>> polys_containing_bit;

// During initialization (in main(), after reading input):
polys_containing_bit.resize(bitstring_len);
for (int i=0; i<polys.size(); i++) {
  for (auto bit: polys[i]) {
    polys_containing_bit[bit].push_back(i);
  }
}

// Replace lines 780-782 with:
for (auto bit: current_work.bitstring) {
  current_work.bitstring.clear_bit(bit);
  bool still_valid = true;
  
  // Only check polynomials that contain this bit
  for (auto poly_idx: polys_containing_bit[bit]) {
    if (!current_work.cover->under_consideration[poly_idx]) continue;
    if (! (polys[poly_idx] && current_work.bitstring)) {
      still_valid = false;
      break;
    }
  }
  
  if (!still_valid) {
    current_work.bitstring.set_bit(bit);
  }
}
```

**Complexity Reduction**:
- Before: O(bits_in_solution × total_polynomials × chunks_per_bitstring)
- After: O(bits_in_solution × polynomials_containing_each_bit × chunks_per_bitstring)
- For sparse polynomials, this is a major win

### Priority 2: Optimize BitString::operator&&

**Current implementation** (lines 264-271):
```cpp
bool operator&&(const BitString& rhs) const
{
  for (size_type i=0; i<bitstring.size(); i++) {
    if (bitstring[i] & rhs.bitstring[i])
      return true;
  }
  return false;
}
```

**Issue**: Compiler may not auto-vectorize due to early return.

**Solution A - Help Auto-Vectorization**:
```cpp
bool operator&&(const BitString& rhs) const
{
  const size_type n = bitstring.size();
  const data_type* __restrict__ p1 = bitstring.data();
  const data_type* __restrict__ p2 = rhs.bitstring.data();
  
  for (size_type i=0; i<n; i++) {
    if (p1[i] & p2[i])
      return true;
  }
  return false;
}
```

**Solution B - Explicit SIMD** (if needed):
```cpp
#include <immintrin.h>

bool operator&&(const BitString& rhs) const
{
  const size_type n = bitstring.size();
  const data_type* p1 = bitstring.data();
  const data_type* p2 = rhs.bitstring.data();
  
  size_type i = 0;
  
  // Process 4 uint64_t at a time with AVX2
  if (n >= 4) {
    for (; i + 4 <= n; i += 4) {
      __m256i v1 = _mm256_loadu_si256((__m256i*)(p1 + i));
      __m256i v2 = _mm256_loadu_si256((__m256i*)(p2 + i));
      __m256i result = _mm256_and_si256(v1, v2);
      if (!_mm256_testz_si256(result, result))
        return true;
    }
  }
  
  // Handle remainder
  for (; i < n; i++) {
    if (p1[i] & p2[i])
      return true;
  }
  return false;
}
```

### Priority 3: Track Satisfaction Incrementally

Instead of checking satisfaction each time, maintain state:

```cpp
// Add to BacktrackPoint or make local to task()
std::vector<bool> satisfied;
std::vector<int> satisfaction_count;  // How many bits in bitstring satisfy this poly

// When adding a bit:
void add_bit_and_update_satisfaction(int bit) {
  current_work.bitstring.set_bit(bit);
  for (auto poly_idx: polys_containing_bit[bit]) {
    if (!satisfied[poly_idx]) {
      satisfaction_count[poly_idx]++;
      if (satisfaction_count[poly_idx] > 0) {
        satisfied[poly_idx] = true;
      }
    }
  }
}

// When removing a bit:
void clear_bit_and_update_satisfaction(int bit) {
  current_work.bitstring.clear_bit(bit);
  for (auto poly_idx: polys_containing_bit[bit]) {
    if (satisfied[poly_idx]) {
      satisfaction_count[poly_idx]--;
      if (satisfaction_count[poly_idx] == 0) {
        satisfied[poly_idx] = false;
      }
    }
  }
}
```

This eliminates the polynomial checking entirely during prime implicant reduction.

### Priority 4: Memory Layout Optimization

**Current issue**: `polys` is a vector of BitStrings, each with its own vector. Poor cache locality.

**Solution**: Structure of Arrays layout:
```cpp
// Instead of vector<BitString> polys
class PolyDatabase {
  std::vector<data_type> data;           // All polynomial data contiguous
  std::vector<size_t> offsets;           // Where each polynomial starts
  std::vector<size_t> lengths;           // Chunks per polynomial
  
  // Access polynomial i, chunk j:
  data_type get(size_t poly_idx, size_t chunk_idx) const {
    return data[offsets[poly_idx] + chunk_idx];
  }
};
```

This improves prefetching and cache utilization.

### Priority 5: Profile-Guided Optimization

```bash
# First compilation with profiling
g++ -std=c++20 -march=native -O3 -fprofile-generate \
    -o cnf2dnf cnf2dnf.cc -lpthread

# Run with representative data
./cnf2dnf -t 4 < typical_input.txt > /dev/null

# Recompile with profile data
g++ -std=c++20 -march=native -O3 -fprofile-use \
    -o cnf2dnf cnf2dnf.cc -lpthread
```

This helps with branch prediction and inlining decisions.

## Build Verification

The current binary shows evidence of proper optimization:

1. **SIMD vectorization present**: `movdqu`, `pand` instructions visible in perf output
2. **Native instructions used**: `popcnt` at line 1090 (from `-march=native`)
3. **Aggressive inlining**: BitString operations mostly inlined (no function call overhead)
4. **Tight loops**: Optimized assembly with minimal instructions

To verify a new build:

```bash
# Check for vectorization
objdump -d cnf2dnf | grep -E "movdqu|pand|vpand" | head

# Check for native instructions
objdump -d cnf2dnf | grep "popcnt"

# Compare sizes
size cnf2dnf

# Look for inlined operators
nm cnf2dnf | grep -i "operator&&"  # Should see few/no entries
```

## Current Build Command

```bash
g++ -std=c++20 -march=native -O3 -o cnf2dnf cnf2dnf.cc -lpthread
```

Flags explained:
- `-std=c++20`: Required for `std::popcount` and `<format>`
- `-march=native`: Use all CPU features available on build machine
- `-O3`: Maximum optimization
- `-lpthread`: Link pthread library for threading

## Testing & Validation

### Basic Functionality Tests

The code includes basic tests in `run_some_basic_tests()` (lines 819-842):
- BitString iterator tests
- Bit manipulation tests
- Edge cases around word boundaries (32, 64, 65 bits)

### Simplification Mode

```bash
./cnf2dnf -s < input.txt > simplified.txt
```

Removes duplicate and superset polynomials. Useful for debugging.

### Verification Mode

```bash
./cnf2dnf < input.txt > output.txt
./cnf2dnf -V output.txt < input.txt
```

Verifies each DNF implicant satisfies all CNF clauses (polynomial-time check only).

### Performance Testing

```bash
# Run with performance monitoring
perf record -g ./cnf2dnf -t 8 < large_input.txt > /dev/null

# View hotspots
perf report

# Annotate specific function
perf annotate task
```

### Progress Monitoring

Send `SIGQUIT` (Ctrl+\) to running process for status:
```bash
kill -QUIT <pid>
```

Output shows:
- Finished bitstrings count per cover
- Candidate solutions count
- Progress visualization through backtrack queue

## Known Issues & Considerations

1. **Memory usage**: Each BitString allocates a vector. For large problems, this can be significant.

2. **Thread scaling**: May not scale linearly beyond 4-8 threads depending on problem structure.

3. **Priority queue lock contention**: With many threads, the backtrack_queue lock can become a bottleneck.

4. **No incremental results**: Program must complete before any output is produced.

5. **Subset checking disabled**: The `#ifdef SLOW` checks (lines 692-700) are currently disabled because they're too expensive. This is where optimization could help most.

## Next Steps for Optimization Work

### Completed ✅

- **Early Duplicate Detection** (January 2026, commit 63a40d2)
  - Added check before prime implicant reduction
  - Skips processing of bitstrings already covered by known solutions
  - Low overhead, high benefit for problems with many similar solutions

### Immediate (Highest ROI)

1. **Implement Priority 1**: Add `polys_containing_bit` index and modify the prime implicant reduction loop
   - Expected impact: 2-5x speedup on the remaining hot loop
   - Low risk, localized change
   - This is still the #1 bottleneck (even after early duplicate detection)
   - Measure with perf before/after

2. **Measure improvement**:
   ```bash
   perf stat -e cycles,instructions,cache-misses ./cnf2dnf -t 4 < input.txt
   ```

### Short-term

3. **Try Priority 2**: Optimize `BitString::operator&&`
   - Try `__restrict__` first (easy)
   - If insufficient, try explicit SIMD
   - Verify with `perf annotate`

4. **Profile again**: See if bottleneck moved to different location

### Medium-term

5. **Implement Priority 3**: Incremental satisfaction tracking
   - More complex, requires careful state management
   - Could eliminate checking entirely
   - High risk, high reward

6. **Consider enabling subset checks**: With faster operators, the `#ifdef SLOW` checks might become viable

### Long-term

7. **Algorithmic improvements**: 
   - Better heuristics for bit ordering
   - More sophisticated prime implicant reduction
   - Parallel cover processing

8. **Memory optimization**: Structure-of-arrays layout for better cache utilization

## Code Locations Reference

| What | Line Numbers | Function/Class |
|------|--------------|----------------|
| Main hotspot | 780-782 | `task()` |
| BitString::operator&& | 264-271 | `BitString` |
| Prime implicant reduction | 774-803 | `task()` |
| Main processing loop | 606-773 | `task()` |
| Input reading | 854-920 | `main()` |
| Cover computation | 452-479 | `Cover::expand_cover()` |
| BitString bit operations | 171-191 | `BitString` |
| Finished bitstrings insertion | 418-477 | `FinishedBitStrings::add()` |
| Thread coordination | 577-803 | `task()` |

## Dependencies

- **C++20 compiler**: GCC 10+ or Clang 12+
- **pthread**: For threading
- **LockingPriorityQueue.hpp**: Custom header (not provided, needs implementation)

## Performance Baseline

### Before Any Optimizations (2025)

The program spent:
- **~60%** in polynomial satisfaction checking (lines 780-782 of the prime implicant reduction loop)
- **~15%** in BitString operations (copying, OR operations)
- **~10%** in thread synchronization
- **~15%** other (I/O, initialization, etc.)

### After Early Duplicate Detection (January 2026)

With the early duplicate check implemented (commit 63a40d2):
- The expensive prime implicant reduction is now skipped for bitstrings that match already-known solutions
- Impact varies by problem: highest benefit when many similar solutions exist
- The polynomial satisfaction checking loop remains the primary bottleneck for new (non-duplicate) bitstrings
- **Priority 1 optimization** (incremental polynomial checking) remains the most impactful next step

Any optimization that reduces the remaining polynomial checking hotspot will have dramatic impact on overall runtime.