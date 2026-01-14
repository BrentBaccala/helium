# -*- mode: python -*-
#
# Python code to compute minimal associated prime ideals
# of a system of equations using factorization and monotone
# dualization, followed by the GTZ algorithm.
#
# Sage's minimal_associated_primes() can run out of memory on large
# polynomial systems.  In the case, we can run a series of
# simplification steps that produce a family of simpler systems that
# are stored in a PostgreSQL database, processed in parallel, then
# recombined for the final result.
#
# The basic simplication step is to look for simple variable
# substitutions (i.e, if "a=0" is in the system, set a equal to zero
# everywhere in the system), then factor the polynomials, then run
# monotone dualization, which is a cnf-to-dnf (conjunctive normal form
# to disjunctive normal form) conversion to obtain a set of factored
# ideals.  We'll recurse and repeat the basic simplification step
# again on every output ideal.  Once we've obtained an ideal that
# can't be simplified further using this procedure, run the GTZ
# algorithm of minimal_associated_primes() to complete the
# calculation.  The final result is the union of all of the ideals so
# obtained.
#
# INTERACTIVE USAGE:
#
# load('minprimes.sage')
# delete_database()    # only needed to remove a previous run
# create_database()
# SQL_stage1(eqns)
# SQL_stage2_parallel()
# consolidate_ideals(load_prime_ideals())
#
# SQL_stage2_parallel() has a single threaded variant SQL_stage2().
# Additionally, these routines can be run on a cluster,
# using the SQL database as a centralized data store.
#
# SQL_stage1(eqns) runs single-threaded, accepts the initial
# system of equations, runs a single iteration of the
# variable substitution / factor / cnf2dnf algorithm to
# populate multiple systems in SQL database 'staging'.
#
# For speed, it's best to use a custom version of Sage that has an
# optimized simplifyIdeal function, available as the simplifyIdeal
# branch of https://github.com/BrentBaccala/sage.
#
# This algorithm should produce identical results to the standard algorithm,
# assuming that the system is simple enough that the standard algorithm works.
# In particular, the output of the following two commands should be identical:
#
# sorted([j.groebner_basis() for j in ideal(eqns).minimal_associated_primes()])
# sorted([j.groebner_basis() for j in consolidate_ideals(load_prime_ideals())])
#
# by Brent Baccala
#
# first version - August 2019
# latest version - January 2025
#
# no rights reserved; you may freely copy, modify, or distribute this
# program

import itertools
from itertools import *

import subprocess

import threading

import pickle
import io

# for running cnf2dnf tests
import re

import os
import psutil
import datetime

import time

from sage.data_structures.bitset import FrozenBitset

import concurrent.futures

import traceback

# We expect to have twice as many threads as cores due to Intel Hyperthreading,
# so completely occupy the cores while leaving the other threads available.
num_processes = os.cpu_count() / 2

try:
    from clint.textui.progress import Bar
    # like Bar, but don't display ProgressBar at all if the whole event takes less than a second
    class ProgressBar(Bar):
        def __init__(self, label, expected_size):
            self.timer = threading.Timer(float(1.0), lambda: self.show(-1))
            self.timer.start()
            super().__init__(label=label, expected_size=expected_size, hide=True)
        def show(self, i):
            if i == -1:
                self.hide = False
                super().show(self.last_progress)
            else:
                super().show(i)
        def done(self):
            if self.timer:
                self.timer.cancel()
            else:
                super().done()
except ModuleNotFoundError:
    print("clint package not available; no progress bars will be displayed")
    class ProgressBar:
        def __init__(self, label, expected_size):
            pass
        def show(self, i):
            pass
        def done(self):
            pass

postgres_connection_parameters = {
    'user':     'baccala',
}

# If we're on my laptop (samsung, for testing) or the postgres server (edge), use localhost
# Otherwise, connect over IP to the postgres server (192.168.2.201, edge)
def SQL_use_edge():
    postgres_connection_parameters['host'] = '192.168.2.201'

if not os.uname()[1] in ('samsung', 'edge'):
    SQL_use_edge()

try:
    import psycopg2
except ModuleNotFoundError:
    print("psycopg2 package not available; no SQL database support")

def postgres_connect():
    try:
        try:
            global conn, conn2
            conn = psycopg2.connect(**postgres_connection_parameters)
            conn2 = psycopg2.connect(**postgres_connection_parameters)
            print(f"connected to SQL database {conn.info.dbname} at {conn.info.host}")
        except psycopg2.OperationalError as ex:
            print('SQL OperationalError during connection attempt; no SQL database support')
    except NameError as ex:
        if ex.name == 'psycopg2':
            # if we couldn't load psycopg2, we already printed a warning
            pass
        else:
            raise

postgres_connection_parameters['database'] = 'helium-16.6'
postgres_connect()

# Solving a system of polynomials
#
# For small systems, we can just call minimal_associated_primes, but for larger systems, we wish to repeatedly
# apply two simplification steps before calling that function:
#
# 1. Factor all of the polynomials in the system of equations and build a set of systems,
# all with irreducible polynomials, that generates the same variety.
#
# 2. Find simple linear substitutions that allow a variable to be eliminated.

def factor_eqns(eqns, verbose=False):
    retval = []
    if verbose:
        pb = ProgressBar(label='factor equations', expected_size=len(eqns))
    for i, eqn in enumerate(eqns):
        retval.append(tuple(f for f,m in eqn.factor()))
        if verbose:
            pb.show(i)
    if verbose:
        pb.done()
        print()
    return tuple(retval)

def factor_eqn(eqn):
    return eqn.factor()

# After we factor all of the polynomials, we have a system of equations which all have to be zero,
# and each equation is a product of factors, only one of which needs to be zero to make that
# entire equation zero.  From a standpoint of logic theory, treating each factor as a logic
# variable, this is a AND-of-ORs, or a product-of-sums, a conjunctive normal form (CNF).
# We wish to convert this to disjunctive normal form (DNF), an OR-of-ANDs, sum-of-products,
# that will give us multiple systems of irreducible factors.
#
# I've tried several ways to do this, starting with native Python code and later using
# the program "espresso".  Most recently, I've been using a custom C++ program.

def cnf2dnf_external(cnf_bitsets, simplify=False, parallel=False, stats=None, verbose=False):

    cmd = ['./cnf2dnf', '-t', str(num_processes if parallel else 1)]
    if simplify:
        cmd.append('-s')
    time1 = time.time()
    if verbose:
        print(time.ctime(), f"system {stats['identifier']} : cnf2dnf starting")
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr)
    for bs in cnf_bitsets:
        proc.stdin.write(str(bs).encode())
        proc.stdin.write(b'\n')
    proc.stdin.close()
    # There is some concern that if we proc.wait() here, proc.stdout could fill with data and both processes will deadlock.
    # So we read proc.stdout first, then wait for the subprocess to terminate.
    retval = list(FrozenBitset(bs.decode().strip()) for bs in proc.stdout)
    proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, './cnf2dnf')
    time2 = time.time()
    if stats:
        stats['cnf2dnf_time'] += time2-time1
    # Consensus on stack overflow https://stackoverflow.com/questions/4154571
    # is that since the first thing sorted() does is to convert the argument
    # to a list, there's no reason not to convert the generator to a list
    # (seven lines above), then sort it, then convert it to a tuple.
    #
    # One reason we want it as a tuple (the main reason we want it sorted, too),
    # is so that cnf2dnf_checking can easily check it.
    #
    # Interesting little discovery: A FrozenBitset can't be sorted
    # directly.  You have to convert them to strings to get keys that
    # can be compared.
    retval.sort(key = lambda x:str(x))
    return tuple(retval)

# In the literature, this logic operation is called "dualization of a monotone
# Boolean function", which means CNF-to-DNF conversion of a monotone Boolean
# function.  In the monotone case, the output is unique and unordered (so we sort
# to compare it) and running the CNF-to-DNF algorithm again on the output (dualization)
# reproduces the input (up to a simplification step - remove all duplicates and supersets
# from the original input).
#
# The C++ cnf2dnf program has exhibited enough bugs that I now perform dualization
# verification on every cnf2dnf calculation.

def cnf2dnf_checking(cnf_bitsets, parallel=False, stats=None, verbose=False):
    # it might be a generator, so convert it to a list since we're going to loop over it twice
    cnf_bitsets = list(cnf_bitsets)
    retval = cnf2dnf_external(cnf_bitsets, parallel=parallel, stats=stats, verbose=verbose)
    time1 = time.time()
    if verbose:
        print(time.ctime(), f"system {stats['identifier']} : cnf2dnf verifying")
    simplified = cnf2dnf_external(cnf_bitsets, simplify=True, parallel=parallel, stats=None, verbose=False)
    verification = cnf2dnf_external(retval, parallel=parallel, stats=None, verbose=False)
    time2 = time.time()
    if stats:
        stats['cnf2dnf_verification_time'] += time2-time1
    # simplified (the simplified input) should equal verification (the output run back through the algorithm again)
    # since they're both tuples, comparing them is really this simple
    if (simplified != verification):
        raise RuntimeError("cnf2dnf verification failed")
    return retval

# Which version of cnf2dnf should we use?
#
# cnf2dnf          - pure Python version (slow, maybe buggy) - buried back in Git history
# cnf2dnf_external - my own C++ program
# cnf2dnf_checking - front end to cnf2dnf_external that reverses and verifies the calculation (can run much slower)

cnf2dnf = cnf2dnf_external

# CNF to DNF tests
#
# Each test is a pair of lists of strings (input and output)
#
# Converting CNF to DNF in this format is almost exactly inverse to itself,
# so we run these tests in both direction - each list of strings is used
# as both the input and the output.
#
# Some of these lists include duplicates (11 and 11) and subsets (11 and 10)
# for testing purposes in the input, but duplicates and subsets should never
# appear in the output, so a little bit of processing is done on the output
# strings before we check them.

cnf2dnf_tests = [
    (['1'], ['1']),
    (['1', '1'], ['1']),
    (['01', '10'], ['11']),
    (['11'], ['01', '10']),
    (['11', '01'], ['01']),
    (['11', '01', '10'], ['11']),
    (['111'], ['100', '010', '001']),
    # check 5f28bf
    (['100000000', '000000001'], ['100000001']),
    # check 2c332f
    (['10000000', '00000001'], ['10000001']),
    # tests for links and chains
    (['110000000', '011000000', '001100000', '000010000', '000001100', '000000110', '000000011'],
     ['101010101', '010110101', '101011010', '101010110', '011010101', '010111010', '010110110', '011011010', '011010110']),
    (['110000000', '011000000', '001100000', '000010000', '000001100', '000000110', '000000011', '000010100'],
     ['101010101', '010110101', '101011010', '101010110', '011010101', '010111010', '010110110', '011011010', '011010110']),
    (['110000000', '011000000', '001100000', '000001100', '000010000', '000000110', '000000011', '000010100', '000001000'],
     ['101011101', '010111101', '101011010', '011011101', '010111010', '011011010']),
    # test for cover optimization
    (['1110', '0011'],
     ['1001', '0101', '0010']),
    (['01110', '00011', '11000'],
     ['01001', '10101', '10010', '01010']),
]

cnf2dnf_tests2 = [
    ('001000010', '001111000', '100000001', '101000100')
]


def make_unique_and_discard_subsets(iter):
    retval = set(iter)
    discards = []
    # we can't discard while iterating over the set
    for s in retval:
        regex = re.compile(s.replace('0','.'))
        for s2 in retval:
            if s != s2 and regex.match(s2):
                discards.append(s2)
    for s2 in set(discards):
        retval.remove(s2)
    return retval

def test_cnf2dnf(parallel=False):
    for test in cnf2dnf_tests:
        for input,output in ((0,1), (1,0)):
            inputs = [FrozenBitset(s) for s in test[input]]
            outputs = [FrozenBitset(s) for s in make_unique_and_discard_subsets(test[output])]
            actual_outputs = cnf2dnf(inputs, parallel)
            # why do we need to convert to str here for the sort to work?
            if sorted(actual_outputs, key=lambda x:str(x)) != sorted(outputs, key=lambda x:str(x)):
                print("Input:", inputs)
                print("Expected:", sorted(outputs, key=lambda x:str(x)))
                print("Actual:", sorted(actual_outputs, key=lambda x:str(x)))
                raise RuntimeError('cnf2dnf test failed')
    for test in cnf2dnf_tests2:
        # for these tests, we just run the cnf2dnf_checking tests (which are actually pretty good)
        inputs = [FrozenBitset(s) for s in test]
        outputs = cnf2dnf(inputs, parallel)
    print('All cnf2dnf tests passed')

# Once we've built all of the systems, then we do this:
#
# consolidate_ideals(reduce(lambda a,b: a.union(b), [set(ideal(s).minimal_associated_primes()) for s in systems]))
#
# i.e, compute the minimum associated prime ideals of each system, union the results all together,
# and discard any ideals which are a subset of another ideal
#
# For hydrogen-5, this method produces the same set of five ideals that were generated by computing
# the radical and primary decomposition directly on the original system.  My hope is that this
# method will be usable for helium-16.6, where the direct approach runs out of memory, though
# it is slower for hydrogen-5.

# consolidate_ideals(list_of_ideals)
#
# For ideals, A < B if A is a subset of B
#
# We seek to represent our original ideal as an intersection of prime ideals.
# If any ideal is a strict subset of another in the set, we
# can discard the larger of the two without affecting the intersection.

def consolidate_ideals(list_of_ideals):
    consolidated_ideals = []
    for next_ideal in list_of_ideals:
        if any(I < next_ideal for I in consolidated_ideals):
            continue
        consolidated_ideals = [I for I in consolidated_ideals if not next_ideal < I]
        # simplify the visual presentation of the ideal by using a reduced Groebner basis
        consolidated_ideals.append(ideal(next_ideal.groebner_basis()))
    return consolidated_ideals

# The "simplifyIdeal" procedure in Singular's primdec.lib (primary decomposition library) checks
# for equations with simple variable substitutions, but doesn't get all linear relations.
# It's used as a preprocessing step before starting into something like the GTZ algorithm
# to compute a primary decomposition.  Let's do that step here, but also check for the
# more complicated linear relations.
#
# It finds things like v+p()=0, where p() doesn't involve v, but I also want to get
# q()v+p()=0, which can be split into two systems, one where q and p are both zero,
# and the other where v=-p/q.

try:
    # This version of simplifyIdeal (the fastest) requires the input to be a list, not a tuple,
    # and returns a modified version of that list.
    from sage.rings.polynomial.multi_polynomial_ideal_libsingular import simplifyIdeal_libsingular as simplifyIdeal
except ImportError:
    try:
        from sage.libs.singular.function_factory import ff
        singularSimplifyIdeal = ff.simplifyIdeal__lib.simplifyIdealBWB
        def simplifyIdeal(I):
            return singularSimplifyIdeal(ideal(I))
        print("Optimized simplifyIdeal not available; failing back on Singular version")
    except NameError:
        print("Neither optimized nor Singular simplifyIdealBWB available; falling back on slow Python version")
        def simplifyIdeal(I):
            # I should be a list or a tuple, not an ideal
            # returns a pair: a list of equations and a list of substitutions
            # The substitutions are equations with a simple linear term that were eliminated from the first list of equations
            simplifications = []
            for v in I[0].parent().gens():
                for p in I:
                    if p == 0:
                        pass
                    elif p/p.lc() == v:
                        #print(v, "=", 0)
                        I = tuple(map(lambda p: p.subs({v: 0}), I))
                        simplifications.append(v)
                    elif p.degree(v) == 1:
                        q,r = p.quo_rem(v)
                        if r == 0:
                            # We should pick up this case with another run through cnf2dnf
                            # print("reducible polynomial detected")
                            pass
                        elif q.is_constant() and r.number_of_terms() == 1:
                            # polynomial is qv+r; qv+r=0; replace v with -r/q
                            #print(v, "=", -r/q)
                            #start_time = time.time()
                            I = tuple(map(lambda p: p.subs({v: -r/q}), I))
                            #execution_time = time.time() - start_time
                            #print(f'subs done in {execution_time} seconds')
                            simplifications.append(q*v+r)
                            break
            return I,tuple(simplifications)

# We use simple simplifications (factoring polynomials and substituting for linear variables) to split a big
# system of polynomial equations into subsystems, each of which are then pickled and stored into a SQL
# database.  Then we'll come back and simplify each subsystem using Singular's GTZ algorithm.
#
# Each system is separated into two subsets - complex polynomials (no linear terms) and simple polynomials (a linear term in each).
# This is done to simplify the GTZ calculations, which only need to be done on the complex set, as
# once you know the solutions to the complex set, the solutions to the linear set can be easily calculated.

sql_schema='''
CREATE TYPE status AS ENUM ('queued', 'running', 'finished', 'interrupted', 'failed');

-- irreducible varieties

CREATE TABLE prime_ideals (
      identifier INTEGER GENERATED ALWAYS AS IDENTITY,
      ideal BYTEA,                -- a pickle of a list of polynomials
      degree INTEGER,             -- the maximum degree of the polynomials in the ideal
      num INTEGER,                -- how many 'staging' systems generated this prime ideal
      simplified BOOLEAN          -- has this ideal been a subset of another one?
);

CREATE UNIQUE INDEX ON prime_ideals(md5(ideal));

-- tracking from staging to prime_ideals

CREATE TABLE prime_ideals_tracking (
      origin INTEGER,
      destination INTEGER
);

-- systems that have only be partially simplified

CREATE TABLE staging (
      identifier INTEGER GENERATED ALWAYS AS IDENTITY,
      stage INTEGER,              -- the original system is stage 0, each system derived from it is one integer higher
      origin INTEGER,             -- for stage 1, 0; for later stages, which identifier in the previous stage this system came from
      system BYTEA,               -- a pickle of a tuple pair of tuples of polynomials; the first complex, the second simple
      current_status status,
      node VARCHAR,
      pid INTEGER
);

CREATE INDEX ON staging(identifier);
CREATE INDEX ON staging(identifier) where current_status = 'queued' or current_status = 'interrupted';

-- This table contains pickled rings and pickled polynomials, to keep down the size of the pickled systems.
-- We expect lots of duplicate polynomials, so we only want to store each one once.

CREATE TABLE globals (
      identifier INTEGER GENERATED ALWAYS AS IDENTITY,
      pickle BYTEA
);

-- Making "pickle" unique causes the resulting index to exceed a PostgreSQL limit, so we make the md5 hash unique instead.

CREATE UNIQUE INDEX ON globals(md5(pickle));

CREATE INDEX ON globals(identifier);

-- staging_stats will be ALTERed by the Statistics class if it doesn't have a _time field that corresponds to a statistic

CREATE TABLE staging_stats (
      identifier INTEGER,
      node VARCHAR,
      pid INTEGER,
      unpickle_time INTERVAL,
      factor_time INTERVAL,
      cnf2dnf_time INTERVAL,
      save_global_time INTERVAL,
      simplifyIdeal_time INTERVAL,
      insert_into_staging_time INTERVAL,
      total_time INTERVAL,
      memory_utilization BIGINT
);

CREATE INDEX ON staging_stats(identifier);
'''

def delete_database():
    persistent_data.clear()
    persistent_data_inverse.clear()
    with conn.cursor() as cursor:
        cursor.execute("DROP OWNED BY current_user")
    conn.commit()

def create_database():
    with conn.cursor() as cursor:
        cursor.execute(sql_schema)
    conn.commit()

# To keep the size of our pickled objects down, we use the "persistent data" feature of the pickle
# protocol, which allows objects to be tagged with an identifier string that is saved into
# the pickle instead of the pickled object itself.
#
# We pickle the polynomials (and the ring that the polynomials come from) and save them in
# a SQL table called 'globals', then use an integer in the table as the identifier string.
#
# We also have a UNIQUE constraint on the pickled data in 'globals' so that identical
# polynomials are only stored once.
#
# Inserting polynomials into 'globals' during a long-running transaction creates deadlocks
# on the 'globals' index if there are other processes doing the same thing.
#
# To avoid this, we do all of our 'globals' operations on a separate SQL connection (conn2)
# and commit after every polynomial is inserted.  The long-running transaction is on the
# original connection (conn).

persistent_data = {}
persistent_data_inverse = {}

# GlobalWithTag is a wrapper around an object (a polynomial) with an associated tag.
# The tag is going to be a numeric string that matches 'identifier' in the SQL `globals` table.

class GlobalWithTag:
    def __init__(self, obj, tag):
        self.obj = obj
        self.tag = tag

    def __getattr__(self, name):
        # Delegate to the wrapped object
        return getattr(self.obj, name)

# Saves the object (a polynomial) to the SQL `globals` table and returns it wrapped in a GlobalWithTag

def save_global(obj, stats=None):
    if obj in persistent_data_inverse:
        return GlobalWithTag(obj, persistent_data_inverse[obj])
    p = persistent_pickle(obj)
    # See this stackoverflow post: https://stackoverflow.com/questions/34708509/how-to-use-returning-with-on-conflict-in-postgresql
    # for issues with the simple "ON CONFLICT DO NOTHING RETURNING identifier"
    time1 = time.time()
    with conn2:
        with conn2.cursor() as cursor:
            cursor.execute("INSERT INTO globals (pickle) VALUES (%s) ON CONFLICT DO NOTHING", (p,))
            cursor.execute("SELECT identifier FROM globals WHERE pickle = %s", (p,))
            id = cursor.fetchone()[0]
            persistent_data[str(id)] = obj
            persistent_data_inverse[obj] = str(id)
    time2 = time.time()
    if stats:
        stats['save_global_sql_time'] += time2-time1
    return GlobalWithTag(obj, str(id))

# This routine isn't called anywhere (it's just for debugging) because the globals table is large
# and we don't want to load the whole thing into memory.  But this is the idea.

def load_globals():
    with conn2:
        with conn2.cursor() as cursor:
            cursor.execute("SELECT identifier, pickle FROM globals")
            for id, p in cursor:
                obj = GlobalWithTag(pickle.loads(p), str(id))
                persistent_data[str(id)] = obj
                persistent_data_inverse[obj] = str(id)

# Note that no attempt is made to save objects that aren't already in the persistent data tables.
# If you want something saved into the globals table, you have to call save_global() explicitly.

def persistent_id(obj):
    # It's hard to tell if an arbitrary Python object is hashable
    #    See https://stackoverflow.com/a/3460725/1493790
    # The stackoverflow suggestion (try/except) is quite slow
    if isinstance(obj, GlobalWithTag):
        return obj.tag
    if isinstance(obj, sage.rings.ring.Ring) or isinstance(obj, sage.rings.polynomial.multi_polynomial.MPolynomial):
        return persistent_data_inverse.get(obj, None)
    return None

# This is pretty much how the dumps code works, except that it tries first to use an optimized version from _pickle

def persistent_pickle(val):
    src = io.BytesIO()
    p = pickle.Pickler(src)
    p.persistent_id = persistent_id
    p.dump(val)
    return src.getvalue()

def persistent_load(id):
    # We don't start a new transaction here "with conn2" because the unpickling step
    # can trigger a recursive call to persistent_load while the transaction is still open,
    # and that would fail.
    if id not in persistent_data:
        with conn2.cursor() as cursor:
            cursor.execute("SELECT pickle FROM globals WHERE identifier = %s", (int(id),))
            if cursor.rowcount == 0:
                raise pickle.UnpicklingError("Invalid persistent id")
            else:
                obj = unpickle_internal(cursor.fetchone()[0])
                persistent_data[id] = obj
                persistent_data_inverse[obj] = id
    return persistent_data[id]

def unpickle_internal(p):
    dst = io.BytesIO(p)
    up = pickle.Unpickler(dst)
    up.persistent_load = persistent_load
    retval = up.load()
    return retval

def unpickle(p, verbose=False):
    dst = io.BytesIO(p)
    persistent_ids = []
    def persistent_load_counter(id):
        if id not in persistent_data:
            persistent_ids.append(int(id))
    up = pickle.Unpickler(dst)
    up.persistent_load = persistent_load_counter
    up.load()
    if verbose:
        print(time.ctime(), len(persistent_ids), 'persistent_ids')
    first_one = True
    if persistent_ids:
        block_size = 10000
        for blocknum in range((len(persistent_ids)+block_size-1)//block_size):
            block = persistent_ids[blocknum*block_size:(blocknum+1)*block_size]
            with conn2.cursor() as cursor:
                cursor.execute("SELECT identifier, pickle FROM globals WHERE identifier IN %s", (tuple(block),))
                while pickl := cursor.fetchone():
                    if first_one:
                        if verbose:
                            print(time.ctime(), 'pickle fetch from SQL done; unpickling')
                        first_one = False
                    id = pickl[0]
                    obj = unpickle_internal(pickl[1])
                    persistent_data[str(id)] = obj
                    persistent_data_inverse[obj] = str(id)
    retval = unpickle_internal(p)
    conn2.commit()
    return retval

# Statistics is like a dict, but with the ability to call += on non-existent keys,
# so we can do things like stats['cnf2dnf_time'] += time6-time5
#
# Can be initialized by passing it a dictionary.
#
# Also has an insert_into_SQL method that adds time columns to the stats table if they don't already exist,
# and converts the integers stored in '_time' keys into the INTERVAL stored in SQL.

class Statistics(dict):
    def __init__(self, d=None):
        if d:
            assert type(d) == dict
            self.update(d)
    def __getitem__(self, key):
        return self.setdefault(key, 0)
    def insert_into_SQL(self, cursor):
        for k in self:
            if k.endswith('_time'):
                self[k] = datetime.timedelta(seconds = int(self[k]))
                cursor.execute(f"ALTER TABLE staging_stats ADD COLUMN IF NOT EXISTS {k} INTERVAL")
        cursor.execute(f"INSERT INTO staging_stats ({','.join(self.keys())}) VALUES %s", (tuple(self.values()),))

def dropZeros(eqns):
    return tuple(e for e in eqns if e != 0)

def normalize(eqns):
    return tuple(e/e.lc() for e in eqns)

# I'm experimenting with turning these on or off for efficiency
#
# We can dump polynomials to the SQL table 'globals' and only store persistent ids in the 'systems' table.
# This saves disk space because duplicate polynomials are only stored once.
#
# We can save tracking information that tracks which systems came from which stage2 systems.
#
# Once a stage hits stage_processing_time, we let it finish its current calculation, then dump it to SQL

dump_polynomials_to_globals_table = True
depth_first_processing = True
tracking = True
stage_processing_time = 60

def eliminateZeros(I):
    # I should be a list or a tuple of polynomials, not an ideal
    # returns a list of equations after substituting zero for any variables that appear alone in the system
    simplifications = []
    for v in I[0].parent().gens():
        for p in I:
            if p == 0:
                pass
            elif p/p.lc() == v:
                I = tuple(map(lambda p: p.subs({v: 0}), I))
                simplifications.append(v)
    return simplifications + [p for p in I if p != 0]

def polish_system(system, simplifications, origin, stats=None):
    time1 = time.time()
    # should we add simplifications to the system before running GTZ on it?
    minimal_primes = ideal(system + simplifications).minimal_associated_primes()
    time2 = time.time()
    if stats:
        stats['GTZ_time'] += time2-time1
    with conn.cursor() as cursor:
        for I in minimal_primes:
            ss = I.gens()
            for p in ss:
                save_global(p)
            degree = max(p.degree() for p in ss)
            p = persistent_pickle(sorted(ss))
            # This version creates deadlocks after I introduced the code that runs
            # a bunch of processing stages until stage_processing_time is reached
            #
            #cursor.execute("""INSERT INTO prime_ideals (ideal, degree, num) VALUES (%s, %s, 1)
            #                  ON CONFLICT (md5(ideal)) DO UPDATE SET num = prime_ideals.num + 1
            #                  RETURNING identifier""",
            #               (p, int(degree)))
            cursor.execute("""INSERT INTO prime_ideals (ideal, degree, num) VALUES (%s, %s, 1)
                              ON CONFLICT (md5(ideal)) DO NOTHING""",
                           (p, int(degree)))
            cursor.execute("SELECT identifier FROM prime_ideals WHERE ideal = %s", (p,))
            id = cursor.fetchone()[0]
            cursor.execute("""INSERT INTO prime_ideals_tracking (origin, destination) VALUES (%s, %s)""",
                           (origin, id))

def initial_processing_stage(system, initial_simplifications, origin, verbose=False, stats=None):
    if origin == 0:
        stage = 1
    else:
        stage = 2
    if verbose: print(time.ctime(), 'system', origin, ': simplifyIdeal')
    time1 = time.time()
    eqns,s = simplifyIdeal(list(system))
    time2 = time.time()
    if stats:
        stats['simplifyIdeal_time'] += time2 - time1
    if verbose: print(time.ctime(), 'system', origin, ':', len(s), 'simplifications')
    if origin != 0 and len(s) == 0:
        # If origin != 0, then the last thing we did to this system was to factor and cnf2dnf it.
        # If there are now no simplifications, then it's ready to be "polished" with GTZ
        # If origin == 0, then fall through and factor it.  Worst thing that can happen is that
        #    it gets reinserted as system #1 and we come back here with origin == 1
        if verbose: print(time.ctime(), 'system', origin, ': polishing')
        polish_system(system, initial_simplifications, origin, stats=stats)
        time2a = time.time()
        stats['polishing_time'] += time2a - time2
        if verbose: print(time.ctime(), 'system', origin, ': done')
        return (None, None)
    # See comment below for why we like to keep things sorted
    # We need simplifications to be a tuple because we're going to pickle it
    global simplifications
    simplifications = tuple(sorted(initial_simplifications + normalize(s)))
    # We can't save simplifications itself, since persistent_id() only works on rings and polynomials, not tuples
    for s in simplifications:
        save_global(s)
    time3 = time.time()
    if stats:
        stats['save_global_time'] += time3-time2
    eqns = normalize(dropZeros(eqns))
    if len(eqns) == 0:
        if verbose: print(time.ctime(), 'system', origin, ': polishing')
        polish_system(eqns, simplifications, origin, stats=stats)
        time3a = time.time()
        stats['polishing_time'] += time3a - time3
        if verbose: print(time.ctime(), 'system', origin, ': done')
        return (None, None)
    if any(eqn == 1 for eqn in eqns):
        # the system is inconsistent and needs no further processing
        return (None, None)
    # global for debugging purposes
    global eqns_factors
    if verbose: print(time.ctime(), 'system', origin, ': factoring')
    eqns_factors = factor_eqns(eqns)
    time4 = time.time()
    if stats:
        stats['factor_time'] += time4-time3
    # all_factors is global so that the subprocesses in the next two ProcessPools can access it
    global all_factors
    # By sorting all_factors, we ensure that the systems inserted into stage2 are sorted,
    # because iterating over a FrozenBitset generates integers
    # in ascending order, which are then used as indices to all_factors.
    # We like sorted systems because they help us detect duplicate systems and reduce duplicate work.
    all_factors = set(f for l in eqns_factors for f in l)
    all_factors = sorted(all_factors)

    # bitsets is global so the subprocesses in the ProcessPool below can access it.
    # We create cnf_bitsets now in case we want to dump it for debugging purposes.
    # If we aren't debugging, then we won't use cnf_bitsets until we call cnf2dnf further below.
    global bitsets
    cnf_bitsets = [FrozenBitset(tuple(all_factors.index(f) for f in l), capacity=len(all_factors)) for l in eqns_factors]

    # cnf2dnf records its own timing statistics and prints its own status messages
    dnf_bitsets = cnf2dnf(cnf_bitsets, stats=stats, verbose=verbose)

    return (simplifications, tuple(tuple(all_factors[j] for j in bs) for bs in dnf_bitsets))

def dump_to_SQL(eqns, simplifications, origin, stats=None, verbose=False):
    # This is a list that matches all_factors, but the polynomials are tagged. The tags
    # (short strings that map to the polynomials) will be used to form the pickled objects
    # that are saved to SQL.
    eqns_tagged = []
    simplifications_tagged = []

    time4 = time.time()
    for eqn in eqns:
        tagged_eqn = save_global(eqn, stats=stats)
        eqns_tagged.append(tagged_eqn)
    for eqn in simplifications:
        tagged_eqn = save_global(eqn, stats=stats)
        simplifications_tagged.append(tagged_eqn)
    time5 = time.time()
    if stats:
        stats['save_global_time'] += time5-time4

    time6 = time.time()
    if verbose: print(time.ctime(), 'system', origin, ': insert into SQL')
    with conn.cursor() as cursor:
        system = persistent_pickle((tuple(eqns_tagged), tuple(simplifications_tagged)))
        cursor.execute("INSERT INTO staging (system, origin, current_status) VALUES (%s, %s, 'queued')",
                       (system, int(origin)))
    time7 = time.time()
    if stats:
        stats['insert_into_staging_time'] += time7-time6

# To avoid the situation where we're running short processing stages and spending most of our time
# dumping the results to SQL, we recurse on our results if we've spent less than stage_processing_time
# seconds on this processing stage.  Once a stage hits stage_processing_time seconds, then we wait for
# the current stage to finish (i.e, we don't kill it), but then everything else gets dumped to SQL.
# Of course, the remaining stages might run very quickly, but we don't know that, unless
# we're willing to start them running and then kill them, which we currently don't do.

def processing_stage(system, initial_simplifications, origin, start_time=None, verbose=False, stats=None):
    if not start_time:
        start_time = time.time()
    simplifications, subsystems = initial_processing_stage(system, initial_simplifications, origin, verbose=verbose, stats=stats)
    elapsed_time = time.time() - start_time
    if subsystems:
        for subsystem in subsystems:
            # origin 0 is special because it won't polish in initial_processing_stage; instead, it'll loop forever
            # If we've been processing this stage for more than stage_processing_time seconds, dump to SQL
            if origin == 0 or elapsed_time > stage_processing_time:
                dump_to_SQL(subsystem, simplifications, origin, stats=stats, verbose=verbose)
            else:
                processing_stage(subsystem, simplifications, origin, start_time=start_time, verbose=verbose, stats=stats)

def SQL_stage1(eqns):
    # To keep the size of the pickles down, we save the ring as a global since it's referred to constantly.
    save_global(eqns[0].parent())
    stats = Statistics({'identifier' : 0})
    processing_stage(eqns, tuple(), int(0), verbose=True, stats=stats)
    conn.commit()
    print(stats)

def SQL_stage2(requested_identifier=None, verbose=False):
    with conn.cursor() as cursor:
        while True:
            # This post explains the subquery and the use of "FOR UPDATE SKIP LOCKED"
            # https://dba.stackexchange.com/a/69497
            if requested_identifier:
                cursor.execute("""UPDATE staging
                                  SET current_status = 'running', pid = %s, node = %s
                                  WHERE identifier = %s AND ( current_status = 'queued' OR current_status = 'interrupted' )
                                  RETURNING system, identifier""", (os.getpid(), os.uname()[1], int(requested_identifier)) )
                conn.commit()
            else:
                if depth_first_processing:
                    order='DESC'
                else:
                    order='ASC'
                cursor.execute(f"""UPDATE staging
                                   SET current_status = 'running', pid = %s, node = %s
                                   WHERE identifier = (
                                       SELECT identifier
                                       FROM staging
                                       WHERE ( current_status = 'queued' OR current_status = 'interrupted' )
                                       ORDER BY identifier {order}
                                       LIMIT 1
                                       FOR UPDATE SKIP LOCKED
                                       )
                                   RETURNING system, identifier""", (os.getpid(), os.uname()[1]) )
                conn.commit()
            if cursor.rowcount == 0:
                break
            pickled_system, identifier = cursor.fetchone()
            try:
                start_time = time.time()
                if verbose: print(time.ctime(), 'system', identifier, ': unpickling')
                system, simplifications = unpickle(pickled_system, verbose=verbose)
                unpickle_time = time.time() - start_time

                # The keys in stats (a fancy dictionary) must match column names in the 'staging_stats' SQL table,
                # except that keys ending in '_time' are processed differently because they correspond to SQL INTERVALs;
                # the other keys correspond to SQL INTEGERs, BIGINTs, or VARCHARs.  For the '_time' variables,
                # we store seconds in the dictionary, then convert them to datetime.timedelta's right
                # before we pass them to SQL, and any missing '_time' columns are added to the SQL on the fly.

                stats = Statistics({'identifier' : identifier,
                                    'pid' : os.getpid(),
                                    'node' : os.uname()[1],
                                    'unpickle_time' : unpickle_time})

                processing_stage(system, simplifications, identifier, stats=stats)

                cursor.execute("""UPDATE staging
                                  SET current_status = 'finished'
                                  WHERE identifier = %s""", (identifier, ))

                stats['total_time'] = time.time() - start_time
                stats['memory_utilization'] = psutil.Process(os.getpid()).memory_info().rss

                print(stats)

                stats.insert_into_SQL(cursor)

                conn.commit()

            except KeyboardInterrupt:
                conn.rollback()
                cursor.execute("""UPDATE staging
                                  SET current_status = 'interrupted'
                                  WHERE identifier = %s""", (identifier,))
                cursor.execute("""DELETE FROM staging WHERE origin = %s""", (identifier,))
                conn.commit()
                raise
            except:
                conn.rollback()
                cursor.execute("""UPDATE staging
                                  SET current_status = 'failed'
                                  WHERE identifier = %s""", (identifier,))
                cursor.execute("""DELETE FROM staging WHERE origin = %s""", (identifier,))
                conn.commit()
                raise

            # keep our memory down by clearing our cached polynomials
            persistent_data.clear()
            persistent_data_inverse.clear()

# SQL_stage2_parallel uses ProcessPool for parallelization and needs an initializer and a done callback.

def done_callback(future):
    if future.exception():
        # I'm not sure what to do here to get a full traceback printed
        traceback.print_exception(None, future.exception(), None)

def process_pool_initializer():
    # futures can't share the original SQL connection, so create a new one
    # We DON'T want to close the existing connection, since it's still being used by the parent process
    global conn, conn2
    conn = psycopg2.connect(**postgres_connection_parameters)
    conn2 = psycopg2.connect(**postgres_connection_parameters)

def SQL_stage2_parallel(max_workers = num_processes):
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, initializer=process_pool_initializer) as executor:
            futures = [executor.submit(SQL_stage2) for _ in range(max_workers)]
            for future in futures:
                future.add_done_callback(done_callback)
            num_completed = 0
            while num_completed < max_workers:
                concurrent.futures.wait(futures, timeout=1)
                num_completed = tuple(future.done() for future in futures).count(True)
        for future in futures:
            future.result()
    except:
        # The only exception I've actually seen here is a BrokenProcessPool when my attempt to raise CalledProcessError
        # in one of the futures failed due to a TypeError because I didn't call the BrokenProcessPool constructor correctly.
        conn.rollback()
        with conn.cursor() as cursor:
            cursor.execute("""UPDATE staging
                              SET current_status = 'failed'
                              WHERE current_status = 'running' AND node = %s""", (os.uname()[1],))
        conn.commit()
        raise

def load_prime_ideals():
    retval = []
    with conn.cursor() as cursor:
        cursor.execute("SELECT ideal FROM prime_ideals WHERE simplified IS NOT TRUE")
        for sys in cursor:
            retval.append(ideal(unpickle(sys[0])))
    return retval

def simplify_ideals():
    ideals = load_prime_ideals()
    for I in ideals:
        if any(I < any_ideal for any_ideal in ideals):
            # flag I as simplified
            p = persistent_pickle(I.gens())
            with conn.cursor() as cursor:
                cursor.execute("UPDATE prime_ideals SET simplified = TRUE WHERE ideal = %s", (p,))
                print(f"{cursor.rowcount} rows updated")

def simplify_ideals_2():
    ideals = []
    with conn.cursor() as cursor:
        cursor.execute("SELECT identifier, ideal FROM prime_ideals")
        for sys in cursor:
            ideals.append((sys[0], sys[1], ideal(unpickle(sys[1]))))
    for identifier,pickle,I in ideals:
        if any(I < any_ideal for _,_,any_ideal in ideals):
            # flag I as simplified
            p = persistent_pickle(sorted(I.gens()))
            #print(p)
            #print(bytes(pickle))
            with conn.cursor() as cursor:
                #cursor.execute("UPDATE prime_ideals SET simplified = TRUE WHERE ideal = %s", (p,))
                cursor.execute("UPDATE prime_ideals SET simplified = TRUE WHERE identifier = %s", (identifier,))
                print(f"{cursor.rowcount} rows updated")
    conn.commit()

# Various utility functions for debugging the SQL database from the command line

def list_staging():
    with conn.cursor() as cursor:
        cursor.execute("SELECT identifier FROM staging")
        return tuple(sys[0] for sys in cursor)

def load_stage1(identifier, verbose=False):
    with conn:
        with conn.cursor() as cursor:
            cursor.execute("SELECT system FROM staging WHERE identifier = %s", (int(identifier),))
            return unpickle(cursor.fetchone()[0], verbose=verbose)

def SQL_stage2_reset():
    # this will delete everything except the results of running the first pass
    with conn:
        with conn.cursor() as cursor:
            cursor.execute("DELETE FROM staging WHERE origin != 0")
            cursor.execute("DELETE FROM staging_stats")
            cursor.execute("UPDATE staging SET current_status = 'queued', pid = NULL, node = NULL")
