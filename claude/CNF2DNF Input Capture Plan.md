# Plan: Capture CNF2DNF Input Data During Processing

## Overview

Modify `minprimes.sage` to automatically save cnf2dnf input data to disk during normal processing. This captures real-world test cases from production runs without interrupting the workflow.

## Goal

When processing a system (via `SQL_stage2` or `SQL_stage2_parallel`), save the input to every `cnf2dnf` call so we can later use these files for testing and optimization of the cnf2dnf program.

## Key Insight

During `inner_processing_stage()`, cnf2dnf is called recursively multiple times:
1. Initial call on the full system (potentially long-running)
2. Recursive calls on simplified subsystems (usually faster)

**We can't know ahead of time which will be slow**, so we save input from ALL cnf2dnf calls. The last/longest one will be the most interesting for optimization testing.

## Implementation Approach

### Modify `cnf2dnf_external()` function

**Location:** `/home/baccala/src/helium/minprimes.sage`, lines 195-227

**Changes:**
1. Add optional parameter `save_input=True` and `identifier=None`
2. Before calling cnf2dnf, save input bitsets to disk
3. Generate unique filename for each call using identifier + sequence number
4. Continue with normal cnf2dnf execution

### File Naming Strategy

```
cnf2dnf_test_data/
  input_{identifier}_{sequence}.txt      # The CNF input bitstrings
  input_{identifier}_{sequence}_meta.txt # Metadata about this call
```

Where:
- `identifier` = `staging.identifier` from the SQL table
- `sequence` = incrementing counter for each cnf2dnf call during this system's processing (0, 1, 2, ...)

**Example files:**
```
cnf2dnf_test_data/input_42_0.txt      # First cnf2dnf call for system 42
cnf2dnf_test_data/input_42_0_meta.txt
cnf2dnf_test_data/input_42_1.txt      # Second cnf2dnf call (recursive)
cnf2dnf_test_data/input_42_1_meta.txt
cnf2dnf_test_data/input_42_2.txt      # Third call (recursive)
...
```

## Detailed Implementation

### 1. Add Call Counter to Statistics

Track how many times cnf2dnf is called per system:

```python
# In inner_processing_stage(), initialize counter
stats['cnf2dnf_call_count'] = 0
```

### 2. Modify cnf2dnf_external() Signature

```python
def cnf2dnf_external(cnf_bitsets, simplify=False, parallel=False,
                     stats=None, verbose=False, save_input=False, output_dir='./cnf2dnf_test_data/'):
```

### 3. Add Input Saving Logic

**Insert at line ~200 (before calling subprocess):**

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
        f.write(f"  ./cnf2dnf -t {'4' if parallel else '1'} < {input_filename}\n")

    if verbose:
        print(f"Saved CNF input to: {input_filename}")

# Continue with normal cnf2dnf execution...
cmd = ['./cnf2dnf', '-t', str(num_processes if parallel else 1)]
...
```

### 4. Enable Input Saving in processing_stage()

Modify the call to `cnf2dnf_external()` around line 856:

**Before:**
```python
dnf = cnf2dnf_external(cnf_bitsets, simplify=False, parallel=parallel_execution,
                       stats=stats, verbose=verbose)
```

**After:**
```python
dnf = cnf2dnf_external(cnf_bitsets, simplify=False, parallel=parallel_execution,
                       stats=stats, verbose=verbose, save_input=True)
```

### 5. Add Configuration Flag (Optional)

Add a global flag to enable/disable input saving without code changes:

```python
# Near the top of minprimes.sage (after imports, ~line 92)
SAVE_CNF2DNF_INPUTS = True  # Set to False to disable input saving
CNF2DNF_INPUT_DIR = './cnf2dnf_test_data/'
```

Then in `cnf2dnf_external()`:
```python
if SAVE_CNF2DNF_INPUTS and stats and 'identifier' in stats:
    # ... save input logic ...
```

## File Format Specifications

### Input File (`input_{id}_{seq}.txt`)

```
110101
010011
111000
...
```

- Plain text, one bitstring per line
- Identical to what cnf2dnf receives on stdin
- Can be piped directly: `./cnf2dnf < input_42_0.txt`

### Metadata File (`input_{id}_{seq}_meta.txt`)

```
Staging Identifier: 42
Sequence Number: 0
Number of CNF clauses: 156
Number of Boolean variables: 234
Simplify mode: False
Parallel mode: True
Input file: ./cnf2dnf_test_data/input_42_0.txt

To run manually:
  ./cnf2dnf -t 4 < ./cnf2dnf_test_data/input_42_0.txt
```

## Usage

### Enable Input Capture

```python
load('minprimes.sage')

# Set the flag (if using configuration approach)
SAVE_CNF2DNF_INPUTS = True

# Run normal processing - inputs will be saved automatically
SQL_stage2_parallel()
```

### Process Runs Normally

- Systems are fetched from staging table
- Processing happens as usual
- CNF inputs are saved before each cnf2dnf call
- CNF2DNF runs normally, results go to SQL
- Multiple files created per system (one per recursive call)

### Identify Long-Running Cases

After processing completes:

```bash
# Check timing in staging_stats table
SELECT identifier, cnf2dnf_time
FROM staging_stats
WHERE cnf2dnf_time > 60  -- Systems where cnf2dnf took > 60 seconds
ORDER BY cnf2dnf_time DESC;
```

For each interesting identifier, you'll have saved inputs:
```bash
ls cnf2dnf_test_data/input_42_*
# Shows: input_42_0.txt, input_42_1.txt, input_42_2.txt, ...
# The last sequence number is usually the longest-running one
```

## Testing the Captured Inputs

```bash
# Test with single thread
time ./cnf2dnf -t 1 < cnf2dnf_test_data/input_42_0.txt > output.txt

# Test with multiple threads
time ./cnf2dnf -t 8 < cnf2dnf_test_data/input_42_0.txt > output.txt

# Test with optimized version
time ./cnf2dnf_optimized -t 8 < cnf2dnf_test_data/input_42_0.txt > output.txt

# Verify correctness
./cnf2dnf -V output.txt < cnf2dnf_test_data/input_42_0.txt
```

## Performance Impact

**Minimal overhead:**
- Writing text files is fast compared to cnf2dnf execution time
- Only affects systems being processed (not the database)
- Can be disabled by setting `SAVE_CNF2DNF_INPUTS = False`

**Disk space:**
- Most input files are small (< 1 KB)
- Long-running cases may have larger inputs (few MB)
- Can clean up after extracting interesting cases

## Verification Plan

### 1. Test Basic Functionality

```python
# Process a single system
cursor.execute("SELECT identifier FROM staging WHERE current_status = 'queued' LIMIT 1")
test_id = cursor.fetchone()[0]

# Enable saving
SAVE_CNF2DNF_INPUTS = True

# Process it
SQL_stage2(specific_identifier=test_id, verbose=True)

# Check that files were created
import os
files = [f for f in os.listdir('cnf2dnf_test_data') if f.startswith(f'input_{test_id}_')]
print(f"Created {len(files)} files for system {test_id}")
assert len(files) >= 2  # At least one input + one meta file
```

### 2. Verify File Format

```python
# Check input file format
with open(f'cnf2dnf_test_data/input_{test_id}_0.txt') as f:
    lines = f.readlines()
    assert all(set(line.strip()).issubset({'0', '1'}) for line in lines)
    print(f"Input file has {len(lines)} clauses")

# Check metadata file
with open(f'cnf2dnf_test_data/input_{test_id}_0_meta.txt') as f:
    content = f.read()
    assert f"Staging Identifier: {test_id}" in content
    assert "Sequence Number: 0" in content
```

### 3. Test Captured Input with cnf2dnf

```bash
# Should run without errors
./cnf2dnf -t 1 < cnf2dnf_test_data/input_*_0.txt > /tmp/test_output.txt

# Output should be valid bitstrings
head /tmp/test_output.txt
```

### 4. Verify Normal Processing Still Works

```python
# Processing should complete normally
# Results should still be written to SQL
cursor.execute("SELECT current_status FROM staging WHERE identifier = %s", (test_id,))
status = cursor.fetchone()[0]
assert status == 'finished'
```

## Critical Files to Modify

1. **`/home/baccala/src/helium/minprimes.sage`**
   - Modify `cnf2dnf_external()` (~line 195-227)
   - Modify call site in `inner_processing_stage()` (~line 856)
   - Optionally add global configuration flags (~line 92)

## Integration Notes

- **Non-breaking change**: Default `save_input=False` means no behavior change unless explicitly enabled
- **Backward compatible**: Existing code continues to work
- **Independent**: Doesn't affect database operations or parallel processing
- **Configurable**: Can be enabled/disabled with a global flag
- **Safe**: Only writes to filesystem, never modifies SQL

## Future Enhancements

1. **Conditional saving**: Only save if CNF is larger than threshold
2. **Timing annotation**: Add actual cnf2dnf execution time to metadata after completion
3. **Cleanup utility**: Script to remove small/fast cases, keep only interesting ones
4. **Compression**: Automatically gzip large input files
5. **Database logging**: Track which files were saved in a separate SQL table
