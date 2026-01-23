# Plan: Capture CNF2DNF Input Data During Processing

## Overview

Instrument `minprimes.sage` to automatically save cnf2dnf input data to disk during normal processing. This captures real-world test cases from production runs without interrupting the workflow. The cnf2dnf program runs normally and writes results to SQL as usual.

## Requirements

- Save input bitstrings to disk BEFORE each cnf2dnf call
- Execute cnf2dnf normally (don't skip it)
- Save ALL cnf2dnf inputs during a processing stage (including recursive calls)
- Use unique filenames with identifier + sequence number
- Save metadata about each call
- Normal processing continues (results written to SQL)
- Minimal performance impact

## Implementation Approach

### Modify `cnf2dnf_external()` Function

**Location:** `/home/baccala/src/helium/minprimes.sage`, lines 195-227

**Key changes:**
1. Add parameters `save_input=False` and `output_dir='./cnf2dnf_test_data/'`
2. Before calling subprocess, save input bitsets to disk
3. Generate unique filename using identifier + sequence number
4. Continue with normal cnf2dnf execution

### File Naming Strategy

```
cnf2dnf_test_data/
  input_{identifier}_{sequence}.txt      # CNF input bitstrings
  input_{identifier}_{sequence}_meta.txt # Metadata
```

- `identifier` = `staging.identifier` from SQL table
- `sequence` = incrementing counter (0, 1, 2, ...) for each cnf2dnf call during this system's processing

**Example:** System 42 calls cnf2dnf 3 times:
```
cnf2dnf_test_data/input_42_0.txt      # First call
cnf2dnf_test_data/input_42_1.txt      # Second call (recursive)
cnf2dnf_test_data/input_42_2.txt      # Third call (recursive)
```

## Implementation Steps

### 1. Add Configuration Flags

**Location:** Near top of `minprimes.sage` (~line 92, after imports)

```python
# Configuration for CNF2DNF input capture
SAVE_CNF2DNF_INPUTS = True
CNF2DNF_INPUT_DIR = './cnf2dnf_test_data/'
```

### 2. Track Call Sequence in Statistics

**Location:** In `inner_processing_stage()` initialization (~line 813)

```python
# Initialize cnf2dnf call counter
stats['cnf2dnf_call_count'] = 0
```

### 3. Modify cnf2dnf_external() Function

**Location:** Lines 195-227

**Add to function signature:**
```python
def cnf2dnf_external(cnf_bitsets, simplify=False, parallel=False,
                     stats=None, verbose=False, save_input=False, output_dir='./cnf2dnf_test_data/'):
```

**Add before subprocess call (~line 200):**
```python
# Save input if requested
if save_input and stats and 'identifier' in stats:
    identifier = stats['identifier']

    # Track call sequence
    if 'cnf2dnf_call_count' not in stats:
        stats['cnf2dnf_call_count'] = 0
    sequence = stats['cnf2dnf_call_count']
    stats['cnf2dnf_call_count'] += 1

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save input file
    input_filename = os.path.join(output_dir, f'input_{identifier}_{sequence}.txt')
    with open(input_filename, 'w') as f:
        for bs in cnf_bitsets:
            f.write(str(bs))
            f.write('\n')

    # Save metadata
    meta_filename = os.path.join(output_dir, f'input_{identifier}_{sequence}_meta.txt')
    with open(meta_filename, 'w') as f:
        f.write(f"Staging Identifier: {identifier}\n")
        f.write(f"Sequence Number: {sequence}\n")
        f.write(f"Number of CNF clauses: {len(cnf_bitsets)}\n")
        f.write(f"Number of Boolean variables: {cnf_bitsets[0].capacity() if cnf_bitsets else 0}\n")
        f.write(f"Simplify mode: {simplify}\n")
        f.write(f"Parallel mode: {parallel}\n")
        f.write(f"Input file: {input_filename}\n")
        f.write(f"\nTo run manually:\n")
        num_threads = num_processes if parallel else 1
        f.write(f"  ./cnf2dnf -t {num_threads} < {input_filename}\n")

    if verbose:
        print(f"Saved CNF input: {input_filename}")

# Continue with normal cnf2dnf execution...
cmd = ['./cnf2dnf', '-t', str(num_processes if parallel else 1)]
# ... rest of existing code ...
```

### 4. Enable at Call Site

**Location:** In `inner_processing_stage()` ~line 856

**Modify the cnf2dnf_external call:**

Before:
```python
dnf = cnf2dnf_external(cnf_bitsets, simplify=False, parallel=parallel_execution,
                       stats=stats, verbose=verbose)
```

After:
```python
dnf = cnf2dnf_external(cnf_bitsets, simplify=False, parallel=parallel_execution,
                       stats=stats, verbose=verbose,
                       save_input=SAVE_CNF2DNF_INPUTS, output_dir=CNF2DNF_INPUT_DIR)
```

## Critical Files to Modify

### `/home/baccala/src/helium/minprimes.sage`

**Add new function after line 882** (after `unpack_eqns` definition)

**Dependencies:**
- Database connection: `conn` and `conn2` (already global)
- Existing functions: `unpack_eqns()`, `simplifyIdeal()`, `factor_eqns()`
- Global caches: `persistent_data`, `persistent_data_inverse`

## File Format Specifications

### CNF Input File (`cnf2dnf_input_{identifier}.txt`)

```
110101
010011
111000
...
```

- One bitstring per line
- Each line represents a CNF clause (disjunction of factors)
- Bit position corresponds to a factor in `all_factors`
- No headers, no comments, no trailing newlines after last clause
- Compatible with: `./cnf2dnf < file.txt`

### Metadata File (`cnf2dnf_input_{identifier}_meta.txt`)

```
Staging Identifier: 42
Stage (depth): 2
Origin: 17
Status: finished
Number of CNF clauses: 156
Number of Boolean variables: 234
Number of polynomials in system: 45
Number of simplifications applied: 12

Input file: ./cnf2dnf_test_data/cnf2dnf_input_42.txt

To run cnf2dnf manually:
  ./cnf2dnf -t 1 < ./cnf2dnf_test_data/cnf2dnf_input_42.txt

Timing statistics available in staging_stats table
```

## Usage Examples

### Enable Input Capture

```python
load('minprimes.sage')

# Input capture is enabled by default via SAVE_CNF2DNF_INPUTS = True
# Process normally - inputs will be saved automatically
SQL_stage2_parallel()
```

### Disable Input Capture

```python
# Edit minprimes.sage or set at runtime
SAVE_CNF2DNF_INPUTS = False
SQL_stage2_parallel()  # No inputs will be saved
```

### Identify Interesting Test Cases

```python
# Find long-running cnf2dnf cases
cursor.execute("""
    SELECT identifier, cnf2dnf_time
    FROM staging_stats
    WHERE cnf2dnf_time > 60
    ORDER BY cnf2dnf_time DESC
""")

for row in cursor.fetchall():
    print(f"System {row[0]}: {row[1]}s")
```

### Test Captured Inputs

```bash
# Find all inputs for system 42
ls cnf2dnf_test_data/input_42_*.txt

# Test with single thread
time ./cnf2dnf -t 1 < cnf2dnf_test_data/input_42_0.txt > output.txt

# Test with multiple threads
time ./cnf2dnf -t 8 < cnf2dnf_test_data/input_42_0.txt > output.txt

# Verify correctness
./cnf2dnf -V output.txt < cnf2dnf_test_data/input_42_0.txt
```

## Error Handling

- **File I/O errors**: Wrapped in try/except to prevent processing failure
- **Missing identifier**: Gracefully skip saving (cnf2dnf still runs)
- **Directory creation failure**: Log warning but continue processing
- **Disk space issues**: Will fail gracefully with OS error

## Verification Plan

### 1. Test Basic Functionality

```python
# Enable input capture
SAVE_CNF2DNF_INPUTS = True

# Process a single system
cursor.execute("SELECT identifier FROM staging WHERE current_status = 'queued' LIMIT 1")
test_id = cursor.fetchone()[0]

# Process it
SQL_stage2(specific_identifier=test_id, verbose=True)

# Check that files were created
import os
import glob
files = glob.glob(f'cnf2dnf_test_data/input_{test_id}_*.txt')
print(f"Created {len(files)} input files for system {test_id}")
assert len(files) >= 1  # At least one cnf2dnf call

# Check metadata files exist
meta_files = glob.glob(f'cnf2dnf_test_data/input_{test_id}_*_meta.txt')
assert len(meta_files) == len(files)  # One meta per input
```

### 2. Verify File Format

```python
# Check first input file
input_file = f'cnf2dnf_test_data/input_{test_id}_0.txt'
with open(input_file) as f:
    lines = f.readlines()
    assert all(set(line.strip()).issubset({'0', '1'}) for line in lines)
    print(f"Input file has {len(lines)} CNF clauses")

# Check metadata
meta_file = f'cnf2dnf_test_data/input_{test_id}_0_meta.txt'
with open(meta_file) as f:
    content = f.read()
    assert f"Staging Identifier: {test_id}" in content
    assert "Sequence Number: 0" in content
    print("Metadata file looks good")
```

### 3. Test cnf2dnf Compatibility

```bash
# Should run without errors
./cnf2dnf -t 1 < cnf2dnf_test_data/input_*_0.txt > /tmp/test_output.txt

# Verify output is valid bitstrings
head /tmp/test_output.txt
```

### 4. Verify Normal Processing

```python
# Processing should complete normally
cursor.execute("SELECT current_status FROM staging WHERE identifier = %s", (test_id,))
status = cursor.fetchone()[0]
assert status == 'finished', f"Expected 'finished', got '{status}'"

# Results should be in database
cursor.execute("SELECT COUNT(*) FROM prime_ideals_tracking WHERE origin = %s", (test_id,))
count = cursor.fetchone()[0]
print(f"System {test_id} generated {count} prime ideals")
```

### 5. Performance Impact Test

```python
import time

# Test with input saving enabled
SAVE_CNF2DNF_INPUTS = True
cursor.execute("SELECT identifier FROM staging WHERE current_status = 'queued' LIMIT 5")
test_ids = [row[0] for row in cursor.fetchall()]

start = time.time()
for id in test_ids[:3]:
    SQL_stage2(specific_identifier=id, verbose=False)
time_with_saving = time.time() - start

# Test with input saving disabled
SAVE_CNF2DNF_INPUTS = False
start = time.time()
for id in test_ids[3:]:
    SQL_stage2(specific_identifier=id, verbose=False)
time_without_saving = time.time() - start

overhead = (time_with_saving - time_without_saving) / time_without_saving * 100
print(f"Input saving overhead: {overhead:.1f}%")
# Should be < 5% for typical systems
```

## Integration Notes

**Changes to existing code:**
- Modify `cnf2dnf_external()` signature (add optional parameters)
- Add input saving logic before subprocess call
- Modify one call site in `inner_processing_stage()`
- Add two global configuration variables

**Backward compatibility:**
- Default `save_input=False` means no behavior change
- Can be disabled by setting `SAVE_CNF2DNF_INPUTS = False`
- No changes to database schema or SQL operations
- No impact on parallel processing or results

**Safety:**
- File I/O errors won't crash processing (can wrap in try/except)
- Only writes to filesystem, never modifies SQL
- Minimal performance overhead (< 5%)

## Future Enhancements

1. **Conditional saving**: Only save if CNF size > threshold
2. **Timing annotation**: Add actual execution time to metadata after cnf2dnf completes
3. **Cleanup utility**: Script to keep only long-running cases
4. **Compression**: Automatically gzip large files
5. **Database tracking**: Optional SQL table to log saved files
