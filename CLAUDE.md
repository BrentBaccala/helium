# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Helium is a computational algebra system for solving the non-relativistic Schrödinger equation for hydrogen and helium atoms. It searches for exact solutions using symbolic differential algebra techniques. The code attempts to find exact solutions by plugging trial solutions with free coefficients into differential equations, expanding them, and solving the resulting system of polynomial equations using Groebner bases and ideal decomposition.

The system combines SageMath for symbolic mathematics, C++ for high-performance CNF-to-DNF conversion, and PostgreSQL for distributed processing of polynomial system decomposition.

The core algorithm:
1. Define a trial solution (ansatz) with free coefficients
2. Plug it into the Schrödinger equation and collect coefficients
3. Find minimal associated primes of the resulting polynomial ideal
4. Each prime ideal represents an irreducible solution variety

## Build Commands

Compile the C++ CNF-to-DNF converter:
```bash
g++ -std=c++20 -march=native -O3 -o cnf2dnf cnf2dnf.cc -lpthread
```

## Interactive Usage

**Basic hydrogen calculation:**
```sage
sage
load('helium.sage')
prep_hydrogen(5)
init()
ideal(eqns_RQQ).minimal_associated_primes()
```

**Database-driven parallel computation (for large systems):**
```sage
load('minprimes.sage')
delete_database()    # only needed to remove a previous run
create_database()
SQL_stage1(eqns)
SQL_stage2_parallel()
consolidate_ideals(load_prime_ideals())
```

Use `SQL_stage2()` for single-threaded processing.

**Parallel worker output:**
When running `SQL_stage2_parallel()`, each worker process writes its output (including stats) to a separate log file named `sql_stage2_worker_<pid>.log`. After parallel processing completes, merge the logs:
```sage
merge_worker_logs()  # Creates sql_stage2_combined.log and removes individual worker logs
merge_worker_logs('my_output.log', cleanup=False)  # Custom output file, keep worker logs
```

**Helium solution (spherically symmetric):**
```sage
load('helium.sage')
prep_helium(-7)  # Negative number for spherically symmetric (faster)
init()
ideal(eqns_RQQ).minimal_associated_primes()
```

**Running journal paper computations:**
```sage
load('joca.sage')     # or joca2.sage
# Follow inline documentation
```

## Key Architecture

### Main Components

**helium.sage** - The primary computational engine containing:
- Trial solution construction (`trial_polynomial` function)
- Multiple "ansatz" types (numbered 1-16 with variants) defining different forms of trial solutions
- Preparation functions `prep_hydrogen()` and `prep_helium()` to set up specific problems
- Polynomial ring creation and management (uses both Singular and FLINT backends)
- Symbolic Ring (SR) → convertField → reduceRing → equation system conversion pipeline
- `init()` function that performs the computationally intensive setup

**minprimes.sage** - Advanced ideal decomposition:
- Factorization followed by CNF to DNF conversion to split ideal systems
- PostgreSQL database integration for distributed/parallel computation
- The `cnf2dnf` C++ program for efficient logical formula conversion
- Functions to manage staging, processing, and consolidation of decomposed ideals

**dynamic_pool.py** - Process management infrastructure:
- Dynamic worker pool that responds to Unix signals (SIGUSR1/SIGUSR2)
- Signal-based worker addition/removal for long-running computations
- Fork-based multiprocessing compatible with SageMath

**cnf2dnf.cc** - High-performance C++20 CNF/DNF converter:
- Converts conjunctive normal form (factored polynomials) to disjunctive normal form (ideal intersections)
- Multi-threaded with work-stealing queue architecture
- Custom bitset operations and heuristics for polynomial ideal factorization

### Key Concepts

**Ansatz System**: Each numbered ansatz (1, 5, 5.1, 16, 16.31, etc.) defines a different form for the trial solution. Higher numbers generally represent more complex forms. Negative ansatz numbers use spherically symmetric coordinates for faster computation.

**Variables Classification**:
- `coordinates`: Independent spatial variables (x,y,z or r1,r2,r12)
- `roots`: Square root expressions (r1, r2, r12)
- `coeff_vars`: Free coefficients to solve for (includes E for energy)
- `ODE_vars`: Variables representing ODE solutions (Zeta, DZeta, etc.)
- `alg_exts`: Algebraic extension variables (gamma)

**Homogenization**: Forces specific polynomial coefficients to be 1 to prevent zero polynomials. Can be iterated through all possibilities with `homogenize=-1`.

### Mathematical Pipeline

- `helium.sage` - Main module: defines ansatzen, applies Hamiltonian, creates polynomial systems
  - `prep_hydrogen(n)` / `prep_helium(n)` - Set up problem with ansatz number n
  - `init()` - Initialize polynomial rings and reduce equations
  - Negative ansatz numbers use spherically symmetric Hamiltonian (faster FLINT implementation)

- `minprimes.sage` - Distributed prime ideal computation
  - Stage 1: Variable substitution, factorization, CNF encoding
  - Stage 2: CNF-to-DNF conversion, recursive simplification
  - Final: GTZ algorithm via `minimal_associated_primes()`

### CNF-to-DNF Conversion

The `cnf2dnf.cc` program converts factored polynomial systems from conjunctive normal form (AND of ORs) to disjunctive normal form (OR of ANDs). This NP-complete problem uses heuristics including single-bit extraction, partitioning into independent covers, and multi-threaded parallel processing.

### Process Management

`dynamic_pool.py` provides a signal-based dynamic worker pool:
- `SIGUSR1` - Add a worker process
- `SIGUSR2` - Gracefully stop one worker
- `SIGTERM/SIGINT` - Shutdown all workers

### Database Schema (PostgreSQL)

**Main Tables:**
- `prime_ideals` - Irreducible varieties (pickled polynomial lists)
- `prime_ideals_tracking` - Tracks which staging systems generated which prime ideals
- `staging` - Systems being processed (queued, running, finished, interrupted, failed)
- `staging_stats` - Timing and resource statistics
- `globals` - Pickled rings and polynomials (deduplicated persistent data)

**Configuration:**
- Database name: `helium-16.6` (or pattern `helium-{ansatz}`)
- Default user: `baccala`
- Default hosts: localhost (for 'samsung' or 'edge' machines), or 192.168.2.201 (edge server)

### Connecting to the PostgreSQL database

There are two primary computation environments in use, each with its own PostgreSQL installation.

- 'samsung' is my development laptop.  I usually use it for a fast calculation called "hydrogen-5" for testing.
- 'edge' is a Cisco C200 with 12 cores and 96 GB of RAM.  It's running a long difficult calculation called "helium-16.6".

Both databases are called `helium-16.6`, even though the one on 'samsung' is actually used for "hydrogen-5".

### Custom SageMath build

A custom SageMath build is needed for two reasons:

- faster simplifyIdeal function
- a specialized pickle format

If you get this message:

```
AttributeError: Can't get attribute 'unpickle_MPolynomial_libsingular_bytestring' on <module 'sage.rings.polynomial.multi_polynomial_libsingular' from '/home/baccala/miniforge3/lib/python3.12/site-packages/sage/rings/polynomial/multi_polynomial_libsingular.cpython-312-x86_64-linux-gnu.so'>
```

It's probably because you're not running a custom version of SageMath with a custom unpickling thing that I've written

## Dependencies

- SageMath 9.x+
- PostgreSQL 10+ with psycopg2
- GCC 11+ (C++20 support)
- Optional: clint (progress bars) - `%pip install clint`
- Optional: WolframScript (verification)
- Optional: DifferentialAlgebra (for joca scripts) - `%pip install DifferentialAlgebra`
- Custom Sage branch recommended: https://github.com/BrentBaccala/sage (simplifyIdeal branch)

## Important Notes

- Never use `git add -i`, `git rebase -i`, or other interactive git commands (not supported in this environment)
- The code requires specific Sage versions and may need the simplifyIdeal branch from https://github.com/BrentBaccala/sage for optimal performance
- Database configuration: Use `SQL_use_edge()` to explicitly connect to the edge server (192.168.2.201) when not running on samsung or edge directly

## Ring Structure

**Ring Pipeline**:
- Symbolic Ring (SR) for initial symbolic manipulations
- convertRing/convertField for expression conversion (FLINT or Singular backend)
- idealRing for reduction ideal construction
- reduceRing for final polynomial reduction
- RQQ (QQ coefficients) and R32003 (finite field) for equation systems

Variables are created in the Symbolic Ring (x1/y1/z1/x2/y2/z2), with r1/r2/r12 as square roots. The workflow maps expressions through convertField to reduceRing, reducing modulo reductionIdeal (e.g., r^2 - x^2 - y^2 - z^2).

## Testing

**CNF2DNF Tests:**
```sage
load('minprimes.sage')
test_cnf2dnf()        # single-threaded
test_cnf2dnf(parallel=True)  # multi-threaded
```

**CNF2DNF Round-trip Test:**
```bash
./cnf2dnf-check <input-file>
```

## Performance Considerations

- The `init()` function is the first computationally expensive operation
- FLINT rings are faster than Singular for some operations but don't support Groebner bases
- Ring selection affects performance significantly; code tries to use FLINT when possible
- The `.polynomial(ring=convertField)` method is much faster than constructing via Symbolic Ring
- Negative ansatz numbers avoid roots/Groebner basis calculations and run much faster
- `num_processes` defaults to `cpu_count() / 2` to account for hyperthreading

## Key File Types

- `.sage` files: SageMath Python scripts (run with `sage filename.sage` or `load('filename.sage')`)
- `.sage.py` files: Pre-parsed SageMath scripts
- `.ssi` files: Sage serialized ideals (ASCII text with very long lines)
- `.dict` files: Pickled Python dictionaries
- `.cc`/`.hpp` files: C++20 source code

## Signal Handling

When using `DynamicProcessManager` or parallel processing:
- Send SIGUSR1 to add worker processes
- Send SIGUSR2 to remove worker processes gracefully
- SIGTERM/SIGINT trigger shutdown
