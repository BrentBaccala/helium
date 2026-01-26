# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Helium is a computational algebra system for solving the non-relativistic Schrödinger equation for hydrogen and helium atoms. It combines SageMath for symbolic mathematics, C++ for high-performance CNF-to-DNF conversion, and PostgreSQL for distributed processing of polynomial system decomposition.

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

## Key Architecture

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

- `systems` - Polynomial systems with encodings
- `staging` - Systems being processed
- `staging_stats` - Timing and resource statistics

Database name follows pattern `helium-{ansatz}` (e.g., `helium-16.6`).

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
- Optional: clint (progress bars), WolframScript (verification)
- Custom Sage branch recommended: https://github.com/BrentBaccala/sage (simplifyIdeal branch)

## Ring Structure

Variables are created in the Symbolic Ring (x1/y1/z1/x2/y2/z2), with r1/r2/r12 as square roots. The workflow maps expressions through convertField to reduceRing, reducing modulo reductionIdeal (e.g., r^2 - x^2 - y^2 - z^2).
