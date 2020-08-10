# -*- mode: python -*-
#
# Python code to search for solutions of Hydrogen and Helium
# non-relativistic time-invariant Schrodinger equation
#
# ALGORITHM:
#
# We start with a trial solution with free coefficients (coeff_vars)
# that we hope will yield an exact solution if they are set to the
# correct real numbers.  We don't try to reduce the complexity of this
# system to the point where we have single isolated solutions, so we
# expect our solutions to be higher-dimensional algebraic varieties.
#
# We plug the trial solution into the differential equation that we're
# trying to solve, expand it out, and collect coefficients of like
# terms.  This gives us a set of polynomial equations in the
# coefficients that define an algebraic variety in the space generated
# by the coefficient variables.  We pick a point at random and use a
# numerical root-finding algorithm in the hopes of finding an
# approximate point on that variety.  If successful, we look for
# algebraic relationships between the point's coordinates, hoping to
# find the algebraic equations that define the irreducible component
# that the point lies on.
#
# A big computational bottleneck is the space complexity of expanding
# out the differential equation.  To solve this problem, we use
# multiple processes to expand out the polynomial incrementally.
#
# by Brent Baccala
#
# first version - August 2019
# latest version - Auguest 2020
#
# no rights reserved; you may freely copy, modify, or distribute this
# program
#
# TODO list:
# - collector should parse directly into matrix form
# - compile a better regular expression to parse terms
# - more easily turn multiprocessing on and off
# - allow worker processes on different hosts
# - remove unused code like LU decomposition and the regular expression parser
# - check collected coefficient polynomials to see if they factor
# - automate finding polynomial relations
# - save checkpoints of optimization iterations
# - optimize scipy sparse matrix by vector multiplication
# - optimize creation of multi-vectors
# - allow workers to be added or removed on the fly

# scipy.minimize tries to minimize a scalar function: the sum of the
# squares of the coefficient polynomials, divided by the value of the
# zero variety (sum of a's squared).  It works pretty well for
# hydrogen, but doesn't seem to converge well for helium.

use_scipy_minimize = False

# scipy.root minimizes the least squares of a vector function

use_scipy_root = True

# divide by an exponential instead of a norm

use_divExp = True

# If neither of the use_scipy options are True, we use my (first)
# custom root finding algorithm, based on the Numerical Recipes text.
# We can either do an exact line search by fitting to a rational
# function, use the scipy built-in line search, or default to the NR
# text algorithm.  We also have the option of using my distributed LU
# decomposition code to solve the least squares problem, or use
# scipy's built-in routine.

use_exact_linesearch = False
use_scipy_line_search = True
use_scipy_lstsq = False
use_scipy_lu = True

# If True, use TCP/IP connections to interconnect Python sub-processes,
# otherwise use UNIX sockets.  TCP/IP has much slower connection setups,
# but allows multiple hosts to be used.

use_tcpip_multiprocessing = False

import platform
import glob
import psutil

current_process = psutil.Process(os.getpid())

from itertools import *
import scipy.optimize

from sage.symbolic.operators import add_vararg, mul_vararg

# from python docs
def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)

# from python docs
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))



var('x1,y1,z1,x2,y2,z2')

r1 = sqrt(x1^2+y1^2+z1^2)
r2 = sqrt(x2^2+y2^2+z2^2)
r12 = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)

SRr_s = (SR.var('r1'), SR.var('r2'), SR.var('r12'))

def trial_polynomial(base, coordinates, radii, degree):
    # cterms are coefficient terms
    # rterms are radius terms
    cterms = flatten([combinations_with_replacement(coordinates, d) for d in range(degree+1)])
    # use this 'terms' for real
    terms = list(map(mul, (product(map(mul, cterms), map(mul, powerset(radii))))))
    # use this 'terms' for testing
    # terms = list(map(mul,cterms)) + list(radii)
    coefficients = tuple(var(base+str(c)) for c in range(len(terms)))
    poly = sum([var(base+str(c))*v for c,v in enumerate(terms)])
    return (coefficients, poly)

# Energy constant in Schroedinger's equations
var('E')

def Del(Psi,vars):
    return sum([diff(Psi,v,2) for v in vars])

def finish_prep(ansatz):
    global eq, H, coeff_vars, coordinates, radii
    global zero_variety, zero_variety_masks

    (Avars, A) = trial_polynomial('a', coordinates, radii, 1)
    (Bvars, B) = trial_polynomial('b', coordinates, radii, 1)
    (Cvars, C) = trial_polynomial('c', coordinates, radii, 1)
    (Dvars, D) = trial_polynomial('d', coordinates, radii, 1)
    (Fvars, F) = trial_polynomial('f', coordinates, radii, 1)
    (Gvars, G) = trial_polynomial('g', coordinates, radii, 1)

    coeff_vars = (E,) + Avars + Bvars + Cvars + Dvars + Fvars + Gvars

    Phi = function('Phi')(*coordinates)
    Xi = function('Xi')(*coordinates)
    Chi = function('Chi')(*coordinates)
    DChi = function('DChi')(*coordinates)

    if ansatz == 1:
        Psi = A*Phi
        zero_variety = sum(map(square, Avars))
    elif ansatz == 2:
        Psi = A*Xi
        zero_variety = sum(map(square, Avars))
    elif ansatz == 3:
        Psi = A*Chi
        zero_variety = sum(map(square, Avars))
    elif ansatz == 4:
        Psi = Chi
        zero_variety = sum(map(square, flatten((Cvars, Dvars, Fvars, Gvars)))) * sum(map(square, flatten((Bvars[1:], Cvars))))
    else:
        raise 'Bad ansatz'

    eq = H(Psi) - E*Psi

    # This is where we set the differential relationships that define
    # how Phi, Xi, and Chi differentiate.

    # Phi is an exponential; Phi = e^B, so diff(Phi,B) = Phi and diff(Phi,v) = diff(B,v)*Phi
    # Xi is a logarithm; Xi = ln C, so diff(Xi,C) = 1/C and diff(Xi,v) = diff(C,v)/C
    # Chi is a second-order ODE: C d^2 Chi/dB^2 - D dChi/dB - F Chi - G = 0

    dict1 = {diff(Phi,v): diff(B,v)*Phi for v in coordinates}
    dict1.update({diff(Xi,v): diff(C,v)/C for v in coordinates})
    dict1.update({diff(Chi,v): diff(B,v)*DChi for v in coordinates})

    dict2 = {diff(Phi,v,2): diff(dict1[diff(Phi,v)],v) for v in coordinates}
    dict2.update({diff(Xi,v,2): diff(dict1[diff(Xi,v)],v) for v in coordinates})
    dict2.update({diff(Chi,v,2): diff(B,v,2)*DChi + diff(B,v)^2*(D/C*DChi+F/C*Chi+G/C) for v in coordinates})

    # replace Phi(x1,y1,z1) with Phi to reduce ginac's memory utilization
    eq = eq.subs(dict2).subs(dict1).subs({Phi: SR.var('Phi'), Xi: SR.var('Xi'), Chi: SR.var('Chi'), DChi: SR.var('DChi')})

    # reduce coeff_vars to those which actually appear in the equation
    coeff_vars = tuple(sorted(set(eq.free_variables()).intersection(coeff_vars), key=lambda x:str(x)))

    # We seek to avoid the zero variety; it's the trivial solution to the DE.
    #
    # However, in the simple case where the zero variety is just a
    # product of factors, each the sum of squares of some coefficient
    # variables, then we can compute the zero variety's function by
    # masking off those variables and computing their norm.
    #
    # In fact, that's the only kind of zero variety the code currently
    # supports (thus the assert).

    def variety_to_masklist(v):

        if v == 1:
            return ()
        elif v.operator() is add_vararg:
            vars = tuple(set(v.free_variables()).intersection(coeff_vars))
            return (np.array(tuple(c in vars for c in coeff_vars)), )
        elif v.operator() is operator.pow:
            return variety_to_masklist(v.operands()[0]) * v.operands()[1]
        elif v.operator() is mul_vararg:
            return tuple(flatten(map(variety_to_masklist, v.operands())))
        else:
            raise "Unknown operator"

    zero_variety_masks = variety_to_masklist(zero_variety)

    coeff_vec = np.array(coeff_vars)
    assert zero_variety == mul(sum(map(square, tuple(coeff_vec * mask))) for mask in zero_variety_masks)

def prep_hydrogen(ansatz=1):
    global H, coordinates, radii

    coordinates = (x1,y1,z1)
    radii = (r1,)

    def H(Psi):
        return - 1/2 * Del(Psi,[x1,y1,z1]) - (1/r1)*Psi

    finish_prep(ansatz=ansatz)

def prep_helium(ansatz=4):
    global H, coordinates, radii

    coordinates = (x1,y1,z1, x2,y2,z2)
    radii = (r1,r2,r12)

    def H(Psi):
        return - 1/2 * Del(Psi,[x1,y1,z1]) - 1/2 * Del(Psi,[x2,y2,z2]) - (2/r1)*Psi - (2/r2)*Psi + (1/r12)*Psi

    finish_prep(ansatz=ansatz)


# Now we want to replace all of the sqrt(...) factors with 'r',
# and we use a clever Python trick to build a dictionary
# that maps expressions to variable names.

def varName(var):
    for name,value in globals().items():
        if id(var) == id(value):
            return name
    return None

def mk_maps(radii):
    return {v.operands()[0] : SR.var(varName(v)) for v in radii}

# convert all (x^2+y^2+z^2)^(n/2) expressions to r^n
# What if we have multiple x^2+y^2+z^2 expressions in a single power?
def roots_to_rs(expr):
    if isinstance(expr, Expression) and expr.operator():
       if expr.operator() == operator.pow and bool(expr.operands()[0] in maps):
           return maps[expr.operands()[0]]^(expr.operands()[1] * 2)
       else:
           return expr.operator()(*map(roots_to_rs, expr.operands()))
    else:
       return expr

def create_eq_a():
    # first, build the dictionary that maps expressions like (x1^2+y1^2+z1^2) to variables like r1
    # make 'maps' global to simplify the map function inside roots_to_rs()
    global maps, eq_a
    maps = mk_maps(radii)
    # next, convert all of the roots in the equation to use the r-variables
    eq_a = roots_to_rs(eq)

# Create a polynomial ring to hold our expressions.
#
# Sage does this using Singular, which stores polynomials internally in standard
# form (i.e, fully expanded).

def create_polynomial_ring():
    radii = SR.var('r1,r2,r12')
    ODE_vars = ('Chi', 'DChi')

    global R,F
    R = PolynomialRing(ZZ, names=tuple(flatten((radii, coeff_vars, coordinates, ODE_vars))))
    F = Frac(R)

def recursive_convert(eq, F):
    # Converting a Sage expression from the Symbolic Ring to a polynomial
    # ring or field is more efficient when we do it using this function.
    def recursion(eq):
        if eq.operator() == add_vararg:
            return sum(map(recursion, eq.operands()))
        elif eq.operator() == mul_vararg:
            return mul(map(recursion, eq.operands()))
        else:
            return F(eq)
    return recursion(eq)

def create_lcm_denominator():
    # Find the least common denominator of all of the terms, then
    # clear the denominators and expand out all of the powers.
    # This is faster than expand(eq_a.numerator()).  One test
    # on helium ran in 61 sec, vs 337 sec for expand/numerator.
    global lcm_denominator
    lcm_denominator = lcm(map(denominator, eq_a.operands()))

def create_polynomial_eq():
    # Next, we want to multiple by lcm_denominator and expand out the resulting polynomial,
    # but the expansion step runs into memory exhaustion issues.  Therefore, we want
    # to distribute the multiplication across eq_a's top-level sum and expand out
    # each term individually.
    global polynomial_eq
    polynomial_eq = expand(eq_a * lcm_denominator)
    # polynomial_eq is now a polynomial, but it's got higher powers of r's in it

def reduce_polynomial_eq():
    # Next, convert powers of r's to x,y,z's
    global reduced_polynomial_eq
    sdict = {SR.var(v)^d : (globals()[v]^d, SR.var(v)*globals()[v]^(d-1))[d%2] for d in range(2,8) for v in ('r1','r2','r12')}
    reduced_polynomial_eq = polynomial_eq.subs(sdict)

def extract_ops():
    global ops
    ops = reduced_polynomial_eq.operands()

def expand_poly(p, v):
    r"""
    Use the custom _split_poly method of Singular polynomials to expand
    a polynomial based on variable number `v` (0-based numbering),
    substituting for even powers of the variable using the global
    variable of the same name, which is expected to be a square root
    in the Symbolic Ring, and can therefore be converted into the
    polynomial ring R if it is raised to an even power.

    DESTROYS INPUT POLYNOMIAL since _split_poly does
    """
    splits = p._split_poly(v+1)
    r = R.gens()[v]
    result = 0
    for power,terms in enumerate(splits):
        podd = power % 2
        peven = power - podd
        multiple = R(globals()[str(r)]^peven) * r^podd
        result += multiple * terms
        print("Finished expanding power", power, "in", r)
    return result

def dump_fast(obj, fn):
    r"""
    Dump an object to a file using 'fast' pickling option to avoid
    excess memory consumption for large polynomials.  Results in a
    larger file since the Integer objects that are the polynomial
    coefficients are not stored as references to earlier Integers if
    they are duplicated.
    """
    f = open(fn, 'wb')
    pickler = pickle.Pickler(f)
    pickler.fast = True
    pickler.dump(obj)
    f.close()

def split_and_dump(p, v, base_name, target_size):
    r"""
    Use the custom _split_poly and _split_poly_evenly methods of
    Singular polynomials to split a polynomial based on variable
    number `v` (0-based numbering), dumping to a set of pickle files
    starting with `base_name`, based on an estimate of expanding `v`
    to a maximum of `target_size` terms in each file,

    Used instead of expand_poly() when expand_poly() would exhaust
    available memory.  Output format is a tuple, the first element of
    which is the monomial we've factored out of the polynomial, which
    is the second term.

    DESTROYS INPUT POLYNOMIAL (to prevent input exhaustion)
    """

    r = R.gens()[v]
    expr = globals()[str(r)]

    for power, terms in enumerate(p._split_poly(v+1)):
        num_terms = sum(1 for _ in terms)
        num_partitions = ceil(sum(1 for _ in R(expr^(2*floor(power/2)))) * num_terms / target_size)
        if num_partitions > 1:
            partitions = terms._split_poly_evenly(num_partitions)
            for partnum, part in enumerate(partitions):
                dump_fast((r^power, part), '{}-{}-{}.pickle'.format(base_name, power, partnum))
        else:
            dump_fast((r^power, terms), '{}-{}.pickle'.format(base_name, power))

def expand_file_to_file(infn, outfn):
    f = open(infn, 'rb')
    # If the input pickle wasn't dumped with the 'fast' option, then
    # it'd make sense to allocate a dedicated Unpickler here, so that
    # we could delete it when we're done and drop the memory
    # references to the memoized sub-objects.  Since the code is built
    # to have dumped using 'fast', I don't do that here.
    (monomial, poly) = pickle.load(f)
    f.close()

    exps = tuple(monomial.dict().keys())
    assert len(exps) == 1
    nzps = exps[0].nonzero_positions()
    assert len(nzps) <= 1

    if len(nzps) == 0:
        dump_fast(poly, outfn)
        print('Current RSS: {:6.1f} GB'.format(float(current_process.memory_info().rss/(2^30))))
    else:
        v = nzps[0]
        power = exps[0][v]
        r = R.gens()[v]
        assert r^power == monomial

        podd = power % 2
        peven = power - podd
        multiple = R(globals()[str(r)]^peven) * r^podd
        result = multiple * poly

        dump_fast(result, outfn)

        print('Current RSS: {:6.1f} GB'.format(float(current_process.memory_info().rss/(2^30))))

def expand_fileset(in_prefix, out_prefix):
    for infn in sorted(glob.glob(in_prefix + '*.pickle')):
        (begin, end) = infn.split(in_prefix, 1)
        assert begin == ""
        outfn = out_prefix + end
        print("Expanding", infn, "to", outfn)
        expand_file_to_file(infn, outfn)

import time

def timefunc(func, *args):
    start_time = time.perf_counter()
    func(*args)
    end_time = time.perf_counter()
    print('{:30} {:10.2f} sec'.format(func.__name__, end_time - start_time))

def multi_init():
    timefunc(create_eq_a)
    timefunc(create_lcm_denominator)
    timefunc(create_polynomial_eq)
    timefunc(reduce_polynomial_eq)
    timefunc(extract_ops)

# Now expand out powers, and collect like x,y,z's terms together to
# get a system of polynomials
#
# This is a slow step, so I've tried several different ways to do it.
#
# We perform the expansion step using multiple processes.
#
# Python multithreading isn't useful for parallel processing because
# of Python's global interpreter lock.  Multiple processes are used
# instead, but this introduces considerable overhead from
# serialization.  To avoid unnecessary serialization, we handle
# polynomials as strings as much as possible.
#
# The desire to avoid serialization also leads to an optimization that
# seems a bit wierd at first.  We don't spawn any of the subprocesses
# until we've calculated `ops`, the terms of the polynomial that we
# want to expand.  This avoids having to serialize `ops`, since each
# worker already has a copy.  This also has the negative side-effect
# of preventing workers from running on other machines, since each
# worker has to be forked from the main process.
#
# We have a manager process (separate from the main process) that acts
# as a traffic cop, and worker processes that pass data directly to
# other workers, avoiding turning either the main process or the
# manager process into a bottleneck.  You might think that the main
# process could be used as the manager, but that isn't well supported
# by Python's multithreading library, and there's very little overhead
# associated with running the manager in its own task.
#
# The proxy to the manager process's ManagerClass is stored in a
# global variable `mc`, and from it methods like `get_workers` can be
# used to obtain proxies to the worker processes.
#
# Right now, it works like this:
#
# load('helium.sage')
# prep_hydrogen()
# multi_init()
# multi_expand()

# number of operands to process in each expander worker
blocksize = 1000

# number of collector processes
num_collectors = 1

# number of simultaneous expander processes
num_expanders = 1

# I've found that starting workers as sub-sub-processes from the
# manager subprocess is error prone because you get copies of the
# manager class in the worker process, and this causes problems with
# reference tracking and the like.  It's better to start all of the
# processes "clean" from the master process, so we handshake between
# the master process and the manager to feed new workers to the
# manager as the old ones terminate (to avoid their memory consumption
# growing without bound).  It's also best to index the `managers`
# dictionary using the proxy object's token, rather than the proxy
# object itself, since otherwise the main process is holding proxy
# objects which get duplicated every time it forks.  This triggers
# reference counter increments on the worker processes, and if they're
# busy, the increments can take some time to run, creating a very
# noticeable delay when starting a new worker.

if not 'managers' in vars():
    managers = {}

if use_tcpip_multiprocessing:
    hostname = platform.node() + ".fios-router.home"
    collector_address = (hostname, 0)
    worker_address = (hostname, 0)
    manager_address = (hostname, 50000)
else:
    collector_address = None
    worker_address = None
    manager_address = None

def start_collector():
    collector_manager = BaseManager(address=collector_address)
    collector_manager.start()
    cc = collector_manager.CollectorClass()
    mc.register_collector(cc)
    managers[cc._token] = collector_manager
    return cc

def start_worker():
    worker_manager = BaseManager(address=worker_address)
    worker_manager.start()
    wc = worker_manager.ExpanderClass()
    mc.register_worker(wc)
    managers[wc._token] = worker_manager

def reap_worker():
    wc = mc.get_finished_worker()
    wc.join_threads()
    # Tricky, tricky... we want to destroy wc, then its manager,
    # in that order.  "del managers[wc]" destroys in the wrong
    # order (since wc is still alive after the `del`)
    worker_manager = managers[wc._token]
    del managers[wc._token]
    del wc
    del worker_manager

def multi_expand():

    start_manager_process()
    mc.set_range(len(ops), blocksize)

    for i in range(num_collectors):
        start_collector()

    for i in range(num_expanders):
        start_worker()

    while True:
        reap_worker()
        if mc.thousands_len() > 0:
            start_worker()
        else:
            break

    for i in range(num_expanders-1):
        reap_worker()

    global ccs
    ccs = mc.get_collectors()

    for cc in ccs:
        cc.join_threads()
        cc.convert_to_matrix()
    for cc in ccs:
        cc.join_threads()

import glob

def multi_load(wildcard):
    for fn in glob.glob(wildcard):
        start_collector().load_matrix(fn)
    global ccs
    ccs = mc.get_collectors()

def remove_duplicates():
    r""""
    Remove duplicate equations from collector processes.

    Detect duplicates by feeding each collector a random vector and
    looking for identical results.  It's possible (but hopefully not
    likely) that different polynomials could return identical results;
    nothing is done to detect this case.  Try to remove duplicates
    from the collectors with the largest number of equations first,
    to even out the time required by each collector to run.
    """
    iv = [random.random() for i in coeff_vars]
    global results
    results = [cc.eval_fns(iv) for cc in ccs]
    fetched_results = list(map(lambda x: x.get(), results))
    sizes = [cc.nrows() for cc in ccs]
    value_to_count = dict(zip(*np.unique(np.hstack(fetched_results), return_counts=True)))
    dups = set((u for u,c in value_to_count.items() if c>1))
    generators = [((j,e) for j,e in enumerate(fetched_results[i]) if e in dups) for i in range(len(ccs))]
    indices = [[] for i in range(len(ccs))]
    remaining_iterations = sum(sizes) - len(value_to_count)
    while len(dups) > 0:
        if remaining_iterations % 1000 == 0:
            print(remaining_iterations, 'iterations remaining')
        maxsize = max(sizes)
        try:
            i = (i for i,c in enumerate(sizes) if c == maxsize and generators[i] is not None).__next__()
        except StopIteration:
            break
        try:
            (index, value) = generators[i].__next__()
        except StopIteration:
            generators[i] = None
            continue
        value_to_count[value] -= 1
        sizes[i] -= 1
        if value_to_count[value] == 1:
            dups.remove(value)
        indices[i].append(index)
        remaining_iterations -= 1
    for i in range(len(ccs)):
        ccs[i].delete_rows(indices[i])

import multiprocessing
multiprocessing.current_process().authkey = b"genius"

from multiprocessing import shared_memory

# SharedMemoryManager is a subclass of BaseManager
#
# maybe it should be used here instead of BaseManager (don't know)
#
# Its __init__ method includes the following comment,
# which justifies the next block of code:
#
# bpo-36867: Ensure the resource_tracker is running before
# launching the manager process, so that concurrent
# shared_memory manipulation both in the manager and in the
# current process does not create two resource_tracker
# processes.

from multiprocessing import resource_tracker
resource_tracker.ensure_running()


import logging
logger = multiprocessing.get_logger()
logger.setLevel(logging.INFO)
if len(logger.handlers) == 0:
    multiprocessing.log_to_stderr()

# I want some of longer-running methods in worker processes to return
# immediately after arranging to start processing in background.  This
# decorator causes a method to run asynchronously in a background
# thread, as well as keeping a record of such threads.  The class
# implements a `join_threads` method that joins (i.e, waits for) all
# of these threads, and should be called before terminating the worker
# process.

import threading

def async_method(method):
    def inner(self, *args):
        th = threading.Thread(target = method, args=(self,) + args)
        th.start()
        if hasattr(self, 'threads'):
            self.threads.append(th)
        else:
            self.threads = [th]
    return inner

# The @async_result decorator is used like @async_method, except that
# a proxy object is generated as the immediately returned result.  The
# actual result is obtained by calling the `get` method on the proxy,
# and the time required to run the method is obtained from proxy's
# `time` method.

class AsyncResult:
    def __init__(self, outerself, method, *args):
        self._cond = threading.Condition()
        self._ready = False

        def target(outerself, method, *args):
            start_time = time.time()
            self._result = method(outerself, *args)
            self._time = time.time() - start_time
            self._cond.acquire()
            self._ready = True
            self._cond.notify_all()
            self._cond.release()

        th = threading.Thread(target = target, args=(outerself,method) + args)
        th.start()
        if hasattr(outerself, 'threads'):
            outerself.threads.append(th)
        else:
            outerself.threads = [th]
        time.sleep(0)

    def get(self):
        self._cond.acquire()
        if not self._ready:
            self._cond.wait()
        self._cond.release()
        return self._result

    def time(self):
        self._cond.acquire()
        if not self._ready:
            self._cond.wait()
        self._cond.release()
        return self._time

def async_result(method):
    def inner(self, *args):
        classname = 'AsyncResult'
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        # If we find a local multiprocessing server, use it to create
        # the AsyncResult object and return a proxy to it.  Otherwise,
        # return an AsyncResult object directly.  This allows a
        # decorated method to be called both locally and remotely.
        if server:
            (ident, exposed) = server.create(None, classname, *((self, method) + args))
            token = multiprocessing.managers.Token(typeid=classname, address=server.address, id=ident)
            proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
            server.decref(None, ident)
            return proxy
        else:
            return AsyncResult(self, method, *args)

    return inner

# Standard Python proxy objects (from the multiprocessing library) and
# their tokens can't be used as keys in dictionaries because multiple
# proxies can point to the same object.  This code adjusts the hashing
# and equality methods to ensure that proxies and tokens for the same
# object test equal and can then be used as dictionary keys.
#
# Probably needs to be reported as a bug in the multiprocessing library.
#
# TODO list for Python multiprocessing library:
#   - handle CNTL-C without killing subprocesses
#   - add rconsole to manager/server and use Sage interact for rconsole
#   - hashing proxies and tokens
#   - autoself() and more general auto-objects
#   - async return objects

def BaseProxy_hash(self):
    return hash((self._token.address, self._token.id))
def BaseProxy_eq(self, other):
    return hash(self) == hash(other)

from multiprocessing.managers import BaseProxy
BaseProxy.__hash__ = BaseProxy_hash
BaseProxy.__eq__ = BaseProxy_eq

def Token_hash(self):
    return hash((self.address, self.id))
def Token_eq(self, other):
    return hash(self) == hash(other)

from multiprocessing.managers import Token
Token.__hash__ = Token_hash
Token.__eq__ = Token_eq

# Shared objects in Python's multiprocessing library are usually
# wrapped in proxy objects, but `self` is never a proxy, so how do we
# pass `self` to another shared object?  `autoself` creates a proxy
# object that refers to `self` and can be passed to another shared
# object.
#
# It's a bit inefficient since a multiprocessing server only maintains
# an id-to-obj mapping and now we need to perform an obj-to-id lookup.
# It also assumes that our serializer is 'pickle' (the default), since
# the server doesn't save the name of serializer, though we could
# reverse lookup the type of server.Listener in the listener_client
# table to figure out what the serializer name was.
#
# https://stackoverflow.com/a/55303121/1493790 suggests that this
# might not be needed in newer versions of Python.

# from rfoo.utils import rconsole

class JoinThreads:
    r"""
    Inherit from this class to obtain the `join_threads` method, which
    waits for all @async_method and @async_result methods to complete,
    except for the thread that called `join_threads` (if
    `join_threads` was called from an @async_method or @async_result
    method)
    """

    def join_threads(self):
        if hasattr(self, 'threads'):
            for th in self.threads:
                if th is not threading.currentThread():
                    logger.debug('join %s', th)
                    th.join()

class Autoself(JoinThreads):
    def autoself(self):
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        classname = self.__class__.__name__
        if server:
            for key,value in server.id_to_obj.items():
                if value[0] == self:
                    token = multiprocessing.managers.Token(typeid=classname, address=server.address, id=key)
                    proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
                    return proxy
        else:
            return self

    def autocreate(self, classname, *args, **kwds):
        r"""
        If called from a multiprocessing context, return a proxy
        object that wraps a new object.  The proxy is suitable
        for being returned over a remote procedure call.

        If not called from a multiprocessing context, just
        create the new object and return it.
        """
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        if server:
            (ident, exposed) = server.create(None, classname, *args, **kwds)
            token = multiprocessing.managers.Token(typeid=classname, address=server.address, id=ident)
            proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
            server.decref(None, ident)
            return proxy
        else:
            return globals()[classname](*args, **kwds)

    # convenience functions for development
    def getpid(self):
        return os.getpid()
    def load(self, filename):
        load(filename)
#    def rconsole(self, port=54321):
#        rconsole.spawn_server(port=port)

import queue

class ManagerClass(Autoself):
    # a list of indices waiting to be expanded
    thousands = []
    # a set of running workers
    running_workers = set()
    # a queue of workers that have finished and are waiting for the
    # main process to reap them
    workers_finished = queue.Queue()
    # a list of collectors
    collectors = []
    # a dictionary mapping workers to the number of data transfers
    # remaining until this worker is ready to start a 'combine'
    worker_data_count = {}

    # This can't be put in an __init__ method because autoself()
    # doesn't work until after self has been created and registered
    # with the multiprocessing server.
    def register_self(self):
        global mc
        mc = self.autoself()
    def __del__(self):
        for i in range(len(self.collectors)):
            del self.collectors[-1]

    def set_range(self, limit, blocksize):
        self.thousands.extend(reversed(range(0, limit, blocksize)))
    def register_collector(self, cc):
        cc.register_manager(self.autoself())

        logger.debug('register collector %s', cc._token)
        self.collectors.append(cc)
    def get_collectors(self):
        return self.collectors
    def register_worker(self, wc):
        self.running_workers.add(wc)
        wc.register_manager(self.autoself())
        wc.register_collectors(self.collectors)

        logger.debug('register_worker %s', wc._token)

        if len(self.thousands) > 0:
            thousand = self.thousands.pop()
            logger.debug('%s expand (%d,%d)', self, thousand, thousand+blocksize)
            wc.start_expand(thousand, thousand+blocksize)
        else:
            self.shutdown_worker(wc)

    def get_finished_worker(self):
        return self.workers_finished.get()
    @async_method
    def notify_expand_done(self, wc):
        logger.debug('notify_expand_done %s', wc._token)
        self.running_workers.remove(wc)
        self.workers_finished.put(wc)
    def thousands_len(self):
        return len(self.thousands)

import time

class ExpanderClass(Autoself):
    r"""
    The equation we're trying to solve had now been converted to a big
    polynomial with thousands of terms, each of which needs to be
    expanded out to a sum of monomials.  Asking Sage to do this
    directly results in memory exhaustion, plus we want to do this
    operation in parallel anyway, so this worker process will expand
    out only a subset of the polynomial's terms, split each resulting
    monomial into a key (the variables we're grouping by) and a value
    (the coefficients), and finally transfer the resulting key-value
    dictionaries to a set of collector processes (each key is assigned
    to a single collector process).

    Since everything has to be serialized for IPC transfer to the
    collector processes, it makes sense to handle the monomials
    as strings, as counter-intuitive as that may seem.
    """

    def register_manager(self, mc):
        self.mc = mc
    def register_collectors(self, collectors):
        self.collectors = collectors
        self.dicts = list(map(dict, [[]] * len(collectors)))
    @async_method
    def shutdown(self):
        for i in range(len(self.collectors)):
            del self.collectors[-1]
    @async_method
    def start_expand(self, start, stop):
        # sleep for a fraction of a second here to let this method
        # return back to the manager before we begin work
        time.sleep(float(0.1))
        expr = expand(sum(islice(ops, start, stop)))
        # Each monomial is of the form 2*c0^2*x^3.  We will split on
        # '*' to get the factors.  The coefficient factors are the
        # ones that start with either a number or a coefficient
        # variable.  The rest are key variables.
        coeff_strs = tuple(map(str, coeff_vars)) + tuple("0123456789")
        for monomial in expr.operands():
            s = str(monomial)
            if s[0] == '-':
                sign='-'
                s=s[1:]
            else:
                sign='+'
            # sort to make sure that we can't get separate keys for
            # something like x*y vs y*x.  Probably unnecessary; I
            # suspect that the factors are already consistently
            # sorted.  The collector code currently assumes that any
            # number will be the leading factor, but sorted puts
            # numbers before letters, so we're OK.
            vs = sorted(s.split('*'))
            # the coefficient variables go in the values
            value = '*'.join(filter(lambda x: any([x.startswith(c) for c in coeff_strs]), vs))
            # all the other variables (x,y,z,r,Phi,Xi) form the key
            key = '*'.join(filterfalse(lambda x: any([x.startswith(c) for c in coeff_strs]), vs))
            # hash the key to figure which collector its assigned to
            dictnum = hash(key) % len(self.collectors)
            # XXX the collector code will later pick this stuff apart
            # again using more regex/string operations.  Maybe it
            # would be efficient to combine the sign, number, and
            # variables into a tuple and store these dictionary values
            # as a list of tuples here.
            self.dicts[dictnum][key] = self.dicts[dictnum].get(key, '') + sign + value
        # expansion finished; send the expanded monomials to the collector processes
        for i in range(len(self.collectors)):
            self.collectors[i].combine_data(self.dicts[i])
        self.mc.notify_expand_done(self.autoself())

import json
import pickle

import re

import numpy as np

import scipy.sparse

# from https://stackoverflow.com/questions/46126840
import scipy.sparse as sp
def sp_unique(sp_matrix, axis=0, new_format=None):
    ''' Returns a sparse matrix with the unique rows (axis=0)
    or columns (axis=1) of an input sparse matrix sp_matrix'''
    if axis == 1:
        sp_matrix = sp_matrix.T

    old_format = sp_matrix.getformat()
    dt = sp_matrix.dtype
    ncols = sp_matrix.shape[1]

    if old_format != 'lil':
        sp_matrix = sp_matrix.tolil()

    _, ind = np.unique(sp_matrix.data + sp_matrix.rows, return_index=True)
    rows = sp_matrix.rows[ind]
    data = sp_matrix.data[ind]
    nrows_uniq = data.shape[0]

    sp_matrix = sp.lil_matrix((nrows_uniq, ncols), dtype=dt)  #  or sp_matrix.resize(nrows_uniq, ncols)
    sp_matrix.data = data
    sp_matrix.rows = rows

    if new_format is None:
        new_format = old_format

    ret = sp_matrix.asformat(new_format)
    if axis == 1:
        ret = ret.T
    return ret

# from https://stackoverflow.com/a/26504995
#

def delete_rows_csr(mat, indices):
    """
    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
    """
    if not isinstance(mat, sp.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]

class LUMatrix(JoinThreads):
    r"""
    Proxy class that wraps a matrix and provides methods to operate on
    it to do an LU decomposition.  An instance of this class is
    expected to be one of many, each managing a block of rows in a
    larger matrix.  The matrix is modified in place.  At the end of
    the calculation the remaining matrices in the LUMatrix instances
    hold the lower blocks of the L matrix in the decomposition.

    Optionally, a second matrix can be supplied, which is expected to
    have the same number of rows as the matrix being decomposed.  Row
    pivoting operations are mimicked on the second matrix, and its
    rows are returned along with candidate rows from the main matrix.
    The expected use is that the second matrix is the right-hand side
    of a matrix equation.
    """

    def __init__(self, matrix, f=None):
        self.M = matrix
        self.f = f
        (self.rows, cols) = matrix.shape
        # XXX might be a bit slow
        self.scaling = np.amax(abs(matrix), axis=1)
        self.U = np.zeros((cols, cols))

    def get(self):
        self.join_threads()
        return self.M

    def get_shm_name(self):
        r"""
        Create a shared memory block containing this matrix and return
        its name.  Calling this function multiple times returns the
        same name.  We detach from the block before returning and
        after copying our matrix into the block.  It's the caller's
        responsibility to unlink the block.

        Could be improved.  The matrix could be created so that it's
        already in a shared memory block ready to go.  Not sure about
        the idea of letting the caller unlink the block even though
        this routine might keep returning the name of an unlinked
        block.
        """
        self.join_threads()
        if not hasattr(self, 'shm_name'):
            shm = shared_memory.SharedMemory(create=True, size=self.M.nbytes)
            shm_M = np.ndarray(self.M.shape, dtype=self.M.dtype, buffer=shm.buf)
            shm_M[:] = self.M[:]
            self.shm_name = shm.name
            shm.close()
        return self.shm_name

    def get_shm_info(self):
        return platform.node(), self.get_shm_name(), self.M.shape, self.M.dtype

    def shape(self):
        self.join_threads()
        return self.M.shape

    @async_result
    def LU_step(self, col, oldval, Ucol):
        r"""
        Execute one parallel step in the LU decomposition algorithm.

        First, if col > 0, multiply `oldval` through column `col-1`
        (the last substep in the previous step).

        Then, update our copy of the U matrix by adding the next column,
        that was included in our arguments.

        Next, search our portion of the matrix to find our candidate
        for the next row to pivot.

        Finally, return the index of that row along with the value
        used to select it, and the value that it would produce on the
        diagonal, ordering with the selection value first so that a
        simple `max` operation in the main process will allow it to
        select the best overall candidate.
        """
        self.join_threads()
        if col > 0:
            self.multiply_column(col-1, 1/oldval)
        self.U[:,col] = Ucol
        b = np.array([self.M[row,col] - self.M[row].dot(Ucol) for row in range(self.rows)])
        (val, row, f_val) = self.select_next_row(b)
        return (val, row, b[row], f_val)

    def select_next_row(self, b_vector):
        r"""
        Default routine to select our candidate for the next row to
        pivot.  Picks the row that will produce the largest absolute
        value on the U diagonal, scaled by the maximum value in each
        row.

        INPUT:

        - ``b_vector`` -- a precomputed vector of the diagonal
        elements that would appear on the U matrix for any given row

        """
        abs_b = abs(b_vector)
        # pick row to maximize diagonal, scaled by maximum value in each row
        row = (abs_b / self.scaling).argmax()
        # pick row with largest value of f
        #row = self.f.argmax()
        # pick row with largest value of f/gradient
        #row = (self.f / abs_b).argmax()
        if self.f is not None:
            return (abs_b[row], row, self.f[row])
        else:
            return (abs_b[row], row, None)

    def fetch_and_remove_row(self, row):
        res = self.M[row]
        self.M = np.delete(self.M, row, axis=0)
        self.scaling = np.delete(self.scaling, row)
        if self.f is not None:
            self.f = np.delete(self.f, row)
        self.rows -= 1
        return res

    def multiply_column(self, col, val):
        for row in range(self.rows):
            self.M[row,col] = val * (self.M[row,col] - self.M[row].dot(self.U[:,col]))

class JacobianMatrix(LUMatrix):
    r"""
    Proxy class that wraps a Jacobian matrix.

    The matrix is stored locally on a worker process and the main
    process operates on it using method calls.  The `get` method,
    while provided, should be avoided if the matrix is large.
    Instead, methods are provided that allow a LU decomposition
    algorithm to be implemented in the main process, transferring only
    those rows that need to be permuted to the top.
    """

    @async_method
    def _start_calculation(self, vec, shm_name):
        if not shm_name:
            mdv = np.stack([self.collector.generate_multi_D_vector(vec, var) for var in coeff_vars], axis=1)
        else:
            shm = shared_memory.SharedMemory(shm_name)
            mdv_shape = (self.collector.ncols(), len(coeff_vars))
            mdv = np.ndarray(mdv_shape, dtype=vec.dtype, buffer=shm.buf)

        M = self.collector.dot(mdv)

        if shm_name:
            shm.close()

        f = self.collector.eval_fns(vec).get()
        LUMatrix.__init__(self, M, f)

    def __init__(self, collector, vec, shm_name):
        self.collector = collector
        self._start_calculation(vec, shm_name)

class JacobianDivAMatrix(LUMatrix):
    r"""
    Proxy class that wraps a Jacobian matrix.

    This Jacobian is for the squares of the functions divided by the sum of the squares of the A vectors.

    The idea is to drive the solution away from the hyperplane where
    the A coefficients are all zero.
    """

    @async_method
    def _start_calculation(self, vec):
        # n is number of functions; c is number of variables
        # the values - vec.shape = (c)

        # the functions - N.shape = (n)
        N = self.collector.eval_fns(vec).get()
        # their derivatives - dN.shape = (n,c)
        dN = np.stack([self.collector.dot(self.collector.generate_multi_D_vector(vec, var)) for var in coeff_vars], axis=1)

        # the A values - Av.shape = (c)
        Av = vec * zero_variety_mask
        # sum A^2 - a scalar
        Adenom = sum(Av*Av)

        # output - M.shape = (n,c)
        M = 2*(N.reshape(-1,1))*dN/Adenom - 2*np.outer(N*N,Av)/(Adenom^2)

        LUMatrix.__init__(self, M, N*N/Adenom)

    def __init__(self, collector, vec):
        self.collector = collector
        self._start_calculation(vec)

class CollectorClass(Autoself):
    result = {}
    rows = {}
    indices = {}
    max_degree = 0
    dok = scipy.sparse.dok_matrix((0, 0), np.int64)
    M = None
    multi_D_vec_times = 0
    multi_D_vec_calls = 0
    multi_vec_times = 0
    multi_vec_calls = 0
    dot_times = 0
    dot_calls = 0

    def register_manager(self, mc):
        self.mc = mc

    @async_method
    def combine_data(self, SRd):
        for key,value in SRd.items():
            self.result[key] = self.result.get(key, '') + value

    def dump(self):
        fn = '/tmp/' + str(os.getpid()) + '.json'
        fp = open(fn, 'w')
        json.dump(self.result, fp)
        fp.close()
        logger.info('result dumped to %s', fn)
    @async_method
    def load_json(self, fn):
        fp = open(fn)
        self.result = json.load(fp)
        fp.close()
        logger.info('result loaded from %s', fn)
    @async_method
    def dump_matrix(self):
        fn = '/tmp/' + str(os.getpid()) + '.pickle'
        fp = open(fn, 'wb')
        pickle.dump(self.M, fp)
        fp.close()
        logger.info('matrix dumped to %s', fn)
    @async_method
    def load_matrix(self, fn):
        fp = open(fn, 'rb')
        if self.M is None:
            self.M = pickle.load(fp)
        else:
            self.M = sp_unique(sp.vstack((self.M, pickle.load(fp))))
        fp.close()
        logger.info('matrix loaded from %s', fn)
        logger.info(repr(self.M))

        # compute max_degree based on the number of columns in the loaded matrix,
        # which should match up with the number of elements in a multi_vector
        self.max_degree = 0
        columns = 0
        while columns < self.M.shape[1]:
            self.max_degree += 1
            columns += binomial(self.max_degree+len(coeff_vars)-1,self.max_degree)
        logger.info("max_degree set to %d", self.max_degree)
        if columns != self.M.shape[1]:
            logger.info("Warning: max_degree inconsistent with matrix dimensions: %d != %d", columns, self.M.shape[1])

    def times(self):
        return (self.dot_calls, self.dot_times, self.multi_vec_calls, self.multi_vec_times, self.multi_D_vec_calls, self.multi_D_vec_times)
    def nrows(self):
        return self.M.shape[0]
    def ncols(self):
        return self.M.shape[1]
    def dtype(self):
        return self.M.dtype
    def len_result(self):
        return len(self.result)
    def len_values(self):
        return sum([len(v) for v in a.result.values()])
    def len_keys(self):
        return sum([len(v) for v in a.result.keys()])
    def nterms(self):
        return sum([v.count('+')+v.count('-') for v in self.result.values()])

    # Now we want to evaluate possibly millions of polynomials, as
    # well as their first derivatives.  Using Sage symbolic
    # expressions is slow and consumes much memory, so we use scipy
    # sparse matrices instead.  Each row in the matrix is a polynomial
    # and each column corresponds to a monomial.  The matrix entry is
    # that monomial's coefficient.  We can then evaluate all the
    # polynomials at once by forming a vector of all the monomials (up
    # to a given maximum degree) and multiplying it by the matrix.

    # Given a vector, generate a multi-vector of the vector itself,
    # then all the biproducts of the vector elements, then all the
    # triproducts of the vector elements, etc.
    #
    # generate_multi_vector([a,b,c])
    #   -> [a,b,c,a^2,a*b,a*c,b^2,b*c,c^2,a^3,a^2*b,a^2*c,a*b^2,a*b*c,a*c^2,a*b^2,b^2*c,b*c^2,c^3]

    # XXX don't really need this in a method

    def generate_multi_vector(self, v):
        start_time = time.time()
        npv = np.array(v)
        stack = [npv]
        for d in range(1, self.max_degree):
            stack.append([np.hstack(stack[-1][i:]) * npv[i] for i in range(len(v))])
        res = np.hstack(tuple(flatten(stack)))
        self.multi_vec_times += time.time() - start_time
        self.multi_vec_calls += 1
        return res

    # same idea, but for a first derivative

    def generate_multi_D_vector(self, v, var):
        start_time = time.time()
        ind = coeff_vars.index(var)

        npv = np.array(v)

        firsts = np.zeros(npv.size)
        firsts[ind] = int(1)
        D_stack = [firsts]
        stack = [npv]

        for d in range(1, self.max_degree):
            D_stack.append([np.hstack(stack[-1][i:]) * firsts[i] + np.hstack(D_stack[-1][i:]) * npv[i]
                            for i in range(len(v))])
            stack.append([np.hstack(stack[-1][i:]) * npv[i] for i in range(len(v))])

        res = np.hstack(tuple(flatten(D_stack)))

        self.multi_D_vec_times += time.time() - start_time
        self.multi_D_vec_calls += 1
        return res

    # unused in operation; only for testing

    def verify_D_vector(self):
        return all([all([bool(diff(e,v)==d) for e,d in zip(self.generate_multi_vector(coeff_vars),
                                                           self.generate_multi_D_vector(coeff_vars, v))]) for v in coeff_vars])

    @async_method
    def convert_to_matrix(self):
        self.i = 0
        # I haven't calculated in advance the maximum degree of the monomials,
        # so start at 1 and reshape the matrix every time we hit a monomial
        # that's too big.
        self.max_degree = 1
        self.indices = {str(pair[1]) : pair[0] for pair in enumerate(self.generate_multi_vector(coeff_vars))}
        veclen = len(self.indices)
        self.dok = scipy.sparse.dok_matrix((len(self.result), veclen), np.int64)
        terms_re = re.compile(r'([+-][^+-]+)')
        number_re = re.compile(r'^[0-9]*$')
        for value in self.result.values():
            terms = re.split(terms_re, value)
            for term in filter(bool, terms):
                if term[0] == '-':
                    sign = -1
                    term = term[1:]
                else:
                    sign = 1
                    if term[0] == '+':
                        term = term[1:]
                try:
                    (coeff, monomial) = term.split('*', 1)
                except ValueError:
                    coeff = '1'
                    monomial = term
                if not re.match(number_re, coeff):
                    monomial = coeff + '*' + monomial
                    coeff = '1'
                coeff = sign * int(coeff)

                while True:
                    try:
                        index = self.indices[monomial]
                        break
                    except KeyError:
                        self.max_degree += 1
                        self.indices = {str(pair[1]) : pair[0] for pair in enumerate(self.generate_multi_vector(coeff_vars))}
                        veclen = len(self.indices)
                        self.dok.resize((len(self.result), veclen))

                self.dok[self.i, index] += coeff
            self.i += 1

        #self.M = np.unique(self.M, axis=0)
        #self.M = self.dok.tocsr()
        self.M = sp_unique(self.dok, axis=0, new_format='csr')
        logger.info('convert_to_matrix done')

    def load_from_pickle(self, fn, section=0, total_sections=1):
        r"""
        Loads a polynomial from a pickle file that was dumped using
        the new Singular polynomial pickling code and stores it in a
        sparse matrix.  Doesn't actually create a polynomial; instead,
        overrides the unpickling routine to avoid the memory and CPU
        overhead of converting to a polynomial.

        The idea is call this method repeatedly to load a series of
        pickle files, each containing terms that when added together
        form the entire polynomial equation that we're trying to
        solve.  When we're done, we convert the entire sparse matrix
        from DOK to CSR form (not done in this method).

        The optional "section" arguments allow only a section of the
        row space to be loaded, to reduce the memory footprint.  The
        resulting matrices will either have to be stacked together
        to form a single matrix, or used with code designed to
        distribute the matrix operations over several machines.
        """

        assert type(section) == int
        assert type(total_sections) == int

        # Hardwired indices for Helium Ansatz 4

        #row_indices = (0,1,2, 283,284,285,286,287,288,289,290)
        #coeff_indices = tuple(range(3,283))
        #all_indices = tuple(range(0,291))

        R_coeff_vars = tuple(map(R, coeff_vars))

        from sage.rings.polynomial.polydict import ETuple

        class custom_unpickler:
            def __init__(sself, R):
                sself.R = R

            def __setitem__(sself, mon_tuple, coeff):

                #row_tuple = tuple([mon_tuple[i] for i in row_indices])
                #coeff_tuple = ETuple([(mon_tuple[i] if i in coeff_indices else 0) for i in all_indices])

                # Hardwired indices for Helium Ansatz 4
                (start, mid_plus_end) = mon_tuple.split(2)
                (coeff_tuple, end) = mid_plus_end.split(282)
                row_tuple = start.eadd(end)

                if hash(row_tuple) % total_sections == section:

                    if row_tuple not in self.rows:
                        self.rows[row_tuple] = len(self.rows)
                        self.dok.resize((len(self.rows), len(self.indices)))
                    row = self.rows[row_tuple]

                    while coeff_tuple not in self.indices:
                        self.max_degree += 1
                        self.indices = {next(pair[1]._iter_ETuples())[0] : pair[0] for pair in enumerate(self.generate_multi_vector(R_coeff_vars))}
                        # self.indices = {list(pair[1].dict().keys())[0] : pair[0] for pair in enumerate(self.generate_multi_vector(coeff_vars))}
                        self.dok.resize((len(self.rows), len(self.indices)))
                    index = self.indices[coeff_tuple]

                    self.dok[row, index] += coeff

            def new_MP(sself):
                return sself.R(0)

        sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular_unpickler = custom_unpickler

        f = open(fn, 'rb')
        up = pickle.Unpickler(f)
        up.load()
        f.close()

    def load_from_glob(self, globstr, section=0, total_sections=1):
        for fn in sorted(glob.glob(globstr)):
            print("Loading", fn)
            timefunc(self.load_from_pickle, fn, section, total_sections)
        print('Current RSS: {:6.1f} GB'.format(float(current_process.memory_info().rss/(2^30))))
        print("Converting to CSR")
        self.M = sp_unique(self.dok, axis=0, new_format='csr')
        self.dok = None

    @async_method
    def delete_rows(self, indices):
        logger.info('deleting %d rows', len(indices))
        self.M = delete_rows_csr(self.M, indices)

    # no longer works because I can't figure how to multiply numpy matrices with Sage Expressions in them
    def verify_matrix(self, vars):
        vec = self.generate_multi_vector(vars)
        set1 = set([eval(preparse(result)) for result in self.result])
        set2 = set(list(self.M.dot(vec)))
        return set1 == set2

    def dot(self, multivec):
        start_time = time.time()
        res = self.M.dot(multivec)
        self.dot_times += time.time() - start_time
        self.dot_calls += 1
        return res

    def get_eqn(self, row):
        r"""
        Recreate and retreive the Sage equation corresponding to a
        particular row in the matrix.  Only used for debugging
        """
        B = self.generate_multi_vector(coeff_vars)
        return sum(self.M[row].toarray() * B)

    @async_method
    def compute_partial_mdv(self, vec, shm_name, start, stop):
        shm = shared_memory.SharedMemory(shm_name)
        mdv_shape = (self.M.shape[1], len(coeff_vars))
        mdv = np.ndarray(mdv_shape, dtype=vec.dtype, buffer=shm.buf)
        for i in range(start, stop):
            mdv[:,i] = self.generate_multi_D_vector(vec, coeff_vars[i])
        shm.close()

    @async_result
    def eval_fns(self, vec):
        r"""
        Evaluate our polynomials

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars

        OUTPUT:

        - a vector of all of our polynomials, evaluated at ``vec``
        """
        multivec = self.generate_multi_vector(vec)
        return self.dot(multivec)

    def jac_fns(self, vec, shm_name=None):
        r"""
        Evaluate the Jacobian matrix (the matrix of first-order
        partial derivatives) of our polynomials

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars (ignored in shm_name is given)
        - ``shm_name`` -- a shared memory block with the pre-computed multi-D-vectors

        OUTPUT:

        - the Jacobian matrix, evaluated at ``vec``, as a numpy
        matrix wrapped in a proxy object
        """
        return self.autocreate('JacobianMatrix', self, vec, shm_name)

    def jac_fns_divA(self, vec):
        r"""
        Evaluate the Jacobian matrix (the matrix of first-order
        partial derivatives) of our polynomials divided by
        the sum of the squares of the A-vectors.

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars

        OUTPUT:

        - the Jacobian matrix, evaluated at ``vec``, as a numpy
        matrix wrapped in a proxy object
        """
        return self.autocreate('JacobianDivAMatrix', self, vec)

    @async_result
    def sum_of_squares(self, vec):
        r"""
        Compute the sum of squares of our polynomials

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars

        OUTPUT:

        - sum(p^2), evaluated at ``vec``
        """
        # logger.info('sum_of_squares %s %s', vec.dtype, vec)
        multivec = self.generate_multi_vector(vec)
        return sum(square(self.dot(multivec)))

    @async_result
    def jac_sum_of_squares(self, vec):
        r"""
        Compute the Jacobian vector (aka gradient) of the sum of
        squares of our polynomials

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars

        OUTPUT:

        - the vector of partial derivatives of sum(p^2) w.r.t each coeff_var,
          evaluated at ``vec``

        """
        # d(p^a s^b t^c)/ds = b(p^a s^(b-1) t^c),
        # so (p^a s^b t^c) should map to b(p^a s^(b-1) t^c)
        # dp/ds = 0         ds/ds = 1
        # d(s^2)/ds = 2s    d(ps)/ds = p     d(pt)/ds = 0
        # d(s^3)/ds = 3s^2  d(ps^2)/ds = 2ps   d(pst) = pt
        multivec = self.generate_multi_vector(vec)
        return [sum(2 * self.dot(multivec) * self.dot(self.generate_multi_D_vector(vec, var))) for var in coeff_vars]

# count up which equations have 1st, 2nd, 3rd degree monomials
# (`a` is a CollectionClass obtained from get_obj() in an rconsole session)
#
# ones = np.array(map(lambda l: int(any(l)), a.M.astype(bool).T[:17].T.toarray()))
# twos = np.array(map(lambda l: int(any(l)), a.M.astype(bool).T[17:170].T.toarray()))
# threes = np.array(map(lambda l: int(any(l)), a.M.astype(bool).T[170:].T.toarray()))

# convert first row of a CollectorClass array back to a Sage expression
#
# B = a.generate_multi_vector(coeff_vars)
# e = sum(a.M[0].toarray() * B)

def square(x):
    return x*x

# These functions are here to make it easier to work with long running
# worker processes without restarting them.  `get_objs` returns a
# dictionary containing all of the objects being managed by the
# multiprocessing server, and since we only create a single *Class
# object in each process, `get_obj` returns that object.  Thus,
# you can use the worker's `load` method to load a snippet of
# code that access the worker's *Class object via `get_obj`.
#
# These functions are here to make it easier to update the code of
# worker processes without restarting them.  You 'load' the new code
# into the worker, then 'rconsole' to the worker and bind the new
# class functions to the existing instance.  `get_objs` returns a
# dictionary containing all of the objects being managed by the
# multiprocessing server, and since we only create a single *Class
# object in each process, `get_obj` returns that object.
#
# The bind function isn't working; I cribbed it from stackexchange

def get_objs():
    return multiprocessing.current_process()._manager_server.id_to_obj
def get_obj(cls=Autoself):
    return (v[0] for v in get_objs().values() if isinstance(v[0], cls)).next()

def bind(instance, func, as_name=None):
    """
    Bind the function *func* to *instance*, with either provided name *as_name*
    or the existing name of *func*. The provided *func* should accept the
    instance as the first argument, i.e. "self".
    """
    if as_name is None:
        as_name = func.__name__
    bound_method = func.__get__(instance, instance.__class__)
    setattr(instance, as_name, bound_method)
    return bound_method



# bind = lambda instance, func, asname: setattr(instance, asname, func.__get__(instance, instance.__class__))

# Register all of these classes with BaseManager.  After this step,
# instantiating BaseManager will start a new process in which we can
# request the creation of ManagerClass, ExpanderClass, or CollectorClass
# objects and receive back proxy objects referring to them.  In fact,
# we'll only create a single *Class object in each process, and only a
# single process with a ManagerClass, as all of the other processes
# will be either workers or collectors.

work_queue = queue.Queue()

def get_work():
    return work_queue.get()
def add_work(l):
    for i in l:
        work_queue.put(i)

from multiprocessing.managers import BaseManager
BaseManager.register('ManagerClass', ManagerClass)
BaseManager.register('ExpanderClass', ExpanderClass)
BaseManager.register('CollectorClass', CollectorClass)
BaseManager.register('AsyncResult', AsyncResult)
BaseManager.register('JacobianMatrix', JacobianMatrix)
BaseManager.register('JacobianDivAMatrix', JacobianDivAMatrix)

BaseManager.register('mc', callable = lambda: mc)
BaseManager.register('get_work', get_work)
BaseManager.register('add_work', add_work)

def start_manager_process():
    global manager, mc
    manager = BaseManager(address=manager_address)
    manager.start()
    mc = manager.ManagerClass()
    mc.register_self()


# Look for solutions using an approximate numerical technique

# Standard operator overloading lets us use Sage's multi-precision
# floating point numbers for most numpy operations, but a few need to
# be overridden to use Sage tests.  Use some Python magic to achieve
# this.

import numpy as np

if 'np_isfinite' not in vars():
    np_isfinite = np.isfinite
def bwb_isfinite(x):
    if isinstance(x, SageObject) and isinstance(x.parent(), Field):
        return not x.is_infinity()
    else:
        return np_isfinite(x)
np.isfinite = bwb_isfinite

if 'np_isnan' not in vars():
    np_isnan = np.isnan
def bwb_isnan(x):
    if isinstance(x, SageObject) and isinstance(x.parent(), Field):
        return x.is_NaN()
    else:
        return np_isnan(x)
np.isnan = bwb_isnan


# What type should we use for our numerical approximation?

real_type = np.float64
#real_type = RR
#real_type = RealField(100)

last_time = 0

def get_eqn(row):
    for cc in ccs:
        if row > cc.nrows():
            row = row - cc.nrows()
        else:
            return cc.get_eqn(row)
    raise IndexError("row out of range")


def fn(v):
    r"""
    Evaluate our polynomials at a given coordinate.

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - a numpy vector of all of our polynomials, evaluated at ``vec``

    ALGORITHM:

    Call the `eval_fns` method for all CollectionClass's
    in parallel, then concatenate all of the results together.
    """

    res = np.hstack(tuple(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])))
    return res

def fns_divSqrtA(v):
    r"""
    The vector function we're trying to minimize: the polys that
    define the solution variety, divided by the square root of the
    zero variety we're trying to avoid.
    """

    # Save a copy of vector to aid in stopping and restarting the
    # calculation.  An explicit call to the copy method is required if
    # we're using the Fortran minpack code (i.e, scipy's optimize
    # package) because in that case, 'v' is only a pointer to a
    # Fortran array that gets deallocated once the Fortran code exits.
    # (I think - the copy's definitely needed, though)
    global last_v
    last_v = v.copy()

    res = np.hstack(list(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])))

    denom = mul(map(np.linalg.norm, (v * mask for mask in zero_variety_masks)))

    global last_time
    sum_of_squares = sum(square(res/denom))
    if last_time == 0:
        print(sum_of_squares)
    else:
        print("{:<30} {:10.2f} sec".format(sum_of_squares, time.time()-last_time))
    last_time = time.time()
    return res/denom

def fns_divExpA(v):
    r"""
    The vector function we're trying to minimize: the polynomials that
    define the solution variety, divided by (1-exp(-A)), where A is
    the `zero_variety` we're trying to avoid.
    """

    # Save a copy of vector to aid in stopping and restarting the
    # calculation.  An explicit call to the copy method is required if
    # we're using the Fortran minpack code (i.e, scipy's optimize
    # package) because in that case, 'v' is only a pointer to a
    # Fortran array that gets deallocated once the Fortran code exits.
    # (I think - the copy's definitely needed, though)
    global last_v
    last_v = v.copy()

    res = np.hstack(list(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])))

    sqrtA = mul(map(np.linalg.norm, (v * mask for mask in zero_variety_masks)))
    denom = 1 - math.exp(- sqrtA*sqrtA)

    global last_time
    sum_of_squares = sum(square(res/denom))
    if last_time == 0:
        print(sum_of_squares)
    else:
        print("{:<30} {:10.2f} sec".format(sum_of_squares, time.time()-last_time))
    last_time = time.time()
    return res/denom

def jac_fn(v):
    r"""
    Evaluate the Jacobian matrix (the matrix of first-order
    partial derivatives) of our polynomials

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - the Jacobian matrix, evaluated at ``vec``, as a numpy matrix,
    with as many rows as polynomials and as many columns as coeff_vars

    ALGORITHM:

    Call the `jac_fns` method for all CollectionClass's in
    parallel, then concatenate all of the results together
    (and transpose them).
    """

    res = np.vstack(map(lambda x: x.get(), [cc.jac_fns(v) for cc in ccs]))
    return res

def get_nparray(shms):
    r"""
    To be used in a map call.  Returns a function that retrieves a
    numpy array either as the image of shared memory or by actually
    copying the data.  If a shared memory region was used, append the
    SharedMemory object to the `shms` list to be closed and unlinked
    when the array is no longer needed.
    """
    def get(obj):
        (node, shm_name, shape, dtype) = obj.get_shm_info()
        if node == platform.node():
            shm = shared_memory.SharedMemory(shm_name)
            shms.append(shm)
            a = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
            return a
        else:
            return obj.get()
    return get

def jac_fns_divSqrtA(v):
    global N,dN,Av,Adenom
    N = np.hstack(list(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])))

    mdv_start_time = time.time()

    mdv_shm_size = len(coeff_vars) * ccs[0].ncols() * ccs[0].dtype().itemsize
    mdv_shm = shared_memory.SharedMemory(create=True, size=mdv_shm_size)
    coeff_distribution = [int(i*len(coeff_vars)/len(ccs)) for i in range(len(ccs) + 1)]

    for i in range(len(ccs)):
        ccs[i].compute_partial_mdv(v, mdv_shm.name, coeff_distribution[i], coeff_distribution[i+1])

    for cc in ccs:
        cc.join_threads()

    mdv_time = time.time()
    print("   Compute multi-D-vectors {:7.2f} sec".format(mdv_time - mdv_start_time))

    shms = []
    dN = np.vstack(list(map(get_nparray(shms), [cc.jac_fns(v, mdv_shm.name) for cc in ccs])))

    dN_time = time.time()
    print("   Compute dN              {:7.2f} sec".format(dN_time - mdv_time))

    zero_variety_factors = tuple(v * mask for mask in zero_variety_masks)
    denom = mul(map(np.linalg.norm, zero_variety_factors))
    res = dN/denom - sum(np.outer(N,f)/(denom*np.linalg.norm(f)^2) for f in zero_variety_factors)

    jac_time = time.time()
    print("   Compute Jacobian matrix {:7.2f} sec".format(jac_time - dN_time))

    mdv_shm.close()
    mdv_shm.unlink()

    for shm in shms:
        shm.close()
        shm.unlink()

    return res

def jac_fns_divExpA(v):
    global N,dN,Av,Adenom
    N = np.hstack(list(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])))

    mdv_start_time = time.time()

    mdv_shm_size = len(coeff_vars) * ccs[0].ncols() * ccs[0].dtype().itemsize
    mdv_shm = shared_memory.SharedMemory(create=True, size=mdv_shm_size)
    coeff_distribution = [int(i*len(coeff_vars)/len(ccs)) for i in range(len(ccs) + 1)]

    for i in range(len(ccs)):
        ccs[i].compute_partial_mdv(v, mdv_shm.name, coeff_distribution[i], coeff_distribution[i+1])

    for cc in ccs:
        cc.join_threads()

    mdv_time = time.time()
    print("   Compute multi-D-vectors {:7.2f} sec".format(mdv_time - mdv_start_time))

    shms = []
    dN = np.vstack(list(map(get_nparray(shms), [cc.jac_fns(v, mdv_shm.name) for cc in ccs])))

    dN_time = time.time()
    print("   Compute dN              {:7.2f} sec".format(dN_time - mdv_time))

    # this block is the only thing that's different between jac_fns_divSqrtA and jac_fns_divExpA
    zero_variety_factors = tuple(v * mask for mask in zero_variety_masks)
    sqrtA = mul(map(np.linalg.norm, zero_variety_factors))
    denom = 1 - math.exp(- sqrtA*sqrtA)
    res = dN/denom - 2*math.exp(-sqrtA)/denom^2*sum(np.outer(N,f)*sqrtA^2/np.linalg.norm(f)^2 for f in zero_variety_factors)

    jac_time = time.time()
    print("   Compute Jacobian matrix {:7.2f} sec".format(jac_time - dN_time))

    mdv_shm.close()
    mdv_shm.unlink()

    for shm in shms:
        shm.close()
        shm.unlink()

    return res

def sum_of_squares(v):
    return sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))

def minfunc(v):
    r"""
    The function we're trying to minimize: the sum of squares of the
    polys that define the solution variety, divided by the zero
    variety we're trying to avoid.
    """

    # Save a copy of vector to aid in stopping and restarting the calculation
    global last_v
    last_v = v.copy()

    d = dict(zip(coeff_vars, v))
    sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))

    # Compute zero_var on operands() to make the order of addition consistent
    # for testing purposes (otherwise we see variance in LSBs)
    #zero_var = v.dtype.type(zero_variety.subs(d))
    zero_var = v.dtype.type(sum(map(lambda bwb: bwb.subs(d), zero_variety.operands())))

    res = real_type(sum_of_squares / zero_var)

    global last_time
    if last_time == 0:
        print(res)
    else:
        print("{:<30} {:10.2f} sec".format(res, time.time()-last_time))
    last_time = time.time()
    return res

# jac_minfunc(v) - the Jacobian matrix of minfunc()
#
# For testing purposes, we can compute this either by letting the
# collectors compute the jacobian of the sum of the squares and then
# modifying that result to take the denominator into account, or
# we can have the collectors compute the jacobian matrix for the
# individual functions (taking the denominator into account - divA)
# and summing it up.  The result should be the same in either case,
# up to floating point inaccuracies.

use_matrix_jacobian = True

def jac_minfunc(v):
    global last_v
    last_v = v.copy()

    if use_matrix_jacobian:
        res = sum(map(lambda x: np.array(x.get().sum(axis=0)), [cc.jac_fns_divA(v) for cc in ccs]))
    else:
        d = dict(zip(coeff_vars, v))
        sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))
        jac_sum_of_squares = sum(map(lambda x: np.array(x.get()), [cc.jac_sum_of_squares(v) for cc in ccs]))

        # Compute zero_var on operands() to make the order of addition consistent
        # for testing purposes (otherwise we see variance in LSBs)
        #zero_var = v.dtype.type(zero_variety.subs(d))
        zero_var = v.dtype.type(sum(map(lambda bwb: bwb.subs(d), zero_variety.operands())))

        Av = v * zero_variety_mask
        res = ((jac_sum_of_squares*zero_var - 2*np.array(Av)*sum_of_squares)/zero_var^2)

    return res

import random

def LU_decomposition(matrices):
    (rows, cols) = matrices[0].shape()

    L = np.zeros((cols,cols))
    U = np.zeros((cols,cols))
    f = np.zeros(cols)
    diag_val = 0
    for j in range(cols):
        for i in range(j):
            U[i,j] = L[i,j] - sum(L[i,:].dot(U[:,j]))
        worker_results = map(lambda x: (x[0].get(), x[1]), [(m.LU_step(j, diag_val, U[:,j]), m) for m in matrices])
        ((selection_val, row, diag_val, f_val), submatrix) = max(worker_results)
        if diag_val == 0:
            diag_val = 1e-20
        # This is a second RPC exchange that could be avoided by
        # having LU_step() return the row and collapsing the remove
        # operation into the next call to LU_step().
        L[j] = submatrix.fetch_and_remove_row(row)
        U[j,j] = diag_val
        f[j] = f_val
    for j in range(cols):
        L[j,j] = 1
        L[j,j+1:] = 0
    for m in matrices: m.multiply_column(j, 1/diag_val)
    return (L, U, f)

def LU_test1(seed):
    np.random.seed(seed)
    M = np.random.rand(10,3)
    bwb = LUMatrix(M.copy())
    (L,U,f) = LU_decomposition([bwb])
    M1 = np.array(sorted(M.tolist()))
    M2 = np.array(sorted((np.vstack((L,bwb.M)).dot(U)).tolist()))
    return np.isclose(M1,M2).all()

# assert all([LU_test1(s) for s in range(100)])

def LU_test2(seed):
    np.random.seed(seed)
    M = np.random.rand(10,3)
    bwb1 = LUMatrix(M[0:5].copy())
    bwb2 = LUMatrix(M[5:].copy())
    (L,U,f) = LU_decomposition([bwb1, bwb2])
    M1 = np.array(sorted(M.tolist()))
    M2 = np.array(sorted((np.vstack((L,bwb1.M,bwb2.M)).dot(U)).tolist()))
    return np.isclose(M1,M2).all()

# assert all([LU_test2(s) for s in range(100)])

def optimize_step(vec):
    r"""x

    NOT CURRENTLY USED, as it doesn't converge very quickly near the
    solution, probably due to a poor choice of the direction vector.

    My attempt at an improved line search algorithm, which takes
    advantage that given a direction vector on an algebraic variety,
    the projection of the sum of squares onto that line is itself a
    polynomial, and given a degree bound on the variety's defining
    polynomials we can put a degree bound on the sum of squares
    polynomial, thus we can easily compute a polynomial that doesn't
    just fit the sum of squares, but describes it exactly.

    "Numerical Recipes" recommends against this type of procedure
    on p. 384, section 9.7, "Line Searches and Backtracking":

        Until the early 1970s, standard practice was to choose λ so
        that x new exactly minimizes f in the direction p. However, we
        now know that it is extremely wasteful of function evaluations
        to do so.

    In this case, NR's advice might not apply so readily.  Each step
    requires the calculation of the gradient vector, which for us is
    the most time consuming operation, requiring N(coeff_vars) matrix
    multiplications, while evaluating the function only requires a
    single matrix multiplication.  So for N(coeff_vars) large (more
    than seven, if the variety's degree bound is three), it makes
    sense to do a few more function evaluations at each step.

    """

    # solve J d = f to find a direction vector

    if use_scipy_lstsq:
        jacobian = jac_fns_divSqrtA(vec)
        (evalstep, *_) = scipy.linalg.lstsq(jacobian, fns_divSqrtA(vec))
    elif use_scipy_lu:
        jacobian = jac_fns_divSqrtA(vec)
        (p,l,u) = scipy.linalg.lu(jacobian)
        pb = scipy.linalg.inv(p).dot(fns_divSqrtA(vec))
        n = u.shape[0]
        lu = scipy.linalg.tril(l[0:n], -1) + u
        piv = np.array(range(0,n))
        evalstep = scipy.linalg.lu_solve((lu, piv), pb[0:n])
    else:
        # XXX - our collector class currently doesn't implement jac_fns_divSqrtA
        jacobians = [cc.jac_fns_divSqrtA(vec) for cc in ccs]
        (L, U, f) = LU_decomposition(jacobians)
        lu = np.tril(L, -1) + U
        piv = np.array(range(lu.shape[0]))
        evalstep = scipy.linalg.lu_solve((lu, piv), f)

    # use gradient descent on the sum of squares to find a direction vector

    # v0 = minfunc(vec)
    # gradient = jac_minfunc(vec)
    # norm = sum(square(gradient))
    # evalstep = gradient*v0/norm

    if use_exact_linesearch:
        # We want to sample at seven points to fit a sixth degree polynomial.
        # We expect a zero "close" to -1, so this will sample three points
        # on either size of it
        points = [vec + i*evalstep for i in [-4,-3,-2,-1,0,1,2]]
        values = list(map(sum_of_squares, points))

        global N,D
        # Now fit a polynomial to this data
        N = np.polynomial.polynomial.Polynomial.fit([-4,-3,-2,-1,0,1,2],values,6,[])

        # The denominator is the sum of (A_0 + lambda A_d)^2 for all As
        D = sum([square(np.polynomial.polynomial.Polynomial((vec[i], evalstep[i]))) for i in zero_variety_indices])

        # the numerator of the first derivative of the function we're trying to minimize
        f = (D * N.deriv() - N * D.deriv())

        # find the real roots of the first derivative
        all_roots = np.roots(list(reversed(list(f.coef))))
        roots = [c.real for c in all_roots if np.isclose(c.imag, 0)]

        # find the minimum and its location
        value_root_pairs = [(N(r)/D(r), r) for r in roots]
        value_root = min(value_root_pairs)

        # return the computed step
        nextstep = vec + evalstep*value_root[1]

        # just do this to print the value
        minfunc(nextstep)

    elif use_scipy_line_search:

        (alpha, *_) = scipy.optimize.line_search(minfunc, jac_minfunc, vec, -evalstep)
        nextstep = vec - evalstep*alpha

    else:
        nextstep = vec - evalstep

        current_val = minfunc(vec)
        newton_val = minfunc(nextstep)
        gradient = jac_minfunc(vec)

        g_deriv = - gradient.dot(evalstep)

        if newton_val > current_val - abs(g_deriv)/1000:

            l = - g_deriv / (2*(newton_val - current_val - g_deriv))
            #if l < 0.1:
            #    l = 0.1

            nextstep = vec - l*evalstep

            print('newton step not acceptable - using l =', l)

            nextval = minfunc(nextstep)

            if nextval > current_val - abs(g_deriv)/1000:

                ls = [l, 1.0]
                gs = [nextval, newton_val]

                while nextval > current_val - abs(g_deriv)/1000:
                    # A = g(l) - g'(0) l - g(0)
                    A = [gs[i] - g_deriv * ls[i] - current_val for i in [0,1]]
                    a = (A[0]/ls[0]^2 - A[1]/ls[1]^2) / (ls[0] - ls[1])
                    b = (-A[0]*ls[1]/ls[0]^2 + A[1]*ls[0]/ls[1]^2) / (ls[0] - ls[1])

                    nextl = (- b + sqrt(b^2 - 3*a*g_deriv))/(3*a)
                    nextstep = vec - nextl*evalstep
                    nextval = minfunc(nextstep)

                    ls.insert(0,nextl)
                    gs.insert(0,nextval)

                print('second step not acceptable - using ls = ', ls)

    return nextstep

def make_iv(seed):

    nvars = len(coeff_vars)

    # even though we're using numpy, we don't need to set its PRNG
    # seed, (which would require calling numpy.random.seed()), since
    # the algorithm is deterministic after the iv is picked

    random.seed(seed)        # for random
    set_random_seed(seed)    # for RR.random_element()
    return np.array([random.random() for i in range(nvars)])


def printv(v):
    for pair in zip(coeff_vars, v): print( pair)

def random_numerical(seed=0, limit=None):

    iv = make_iv(seed)

    global SciMin

    # optimize.minimize searchs for minimums of a scalar-valued
    # function, in our case the sum of squares of the variety's
    # defining polynomials.  The problem, as explained in Numerical
    # Recipes at the end of section 9.6, in the subsection titled
    # "Newton’s Method versus Minimization", is that the sum of
    # squares has many local minima that aren't global minima, and any
    # of these minimization algorithms will tend to latch on to a
    # local minimum instead of the global minimum.  The solution is to
    # use an algorithm that searches for zeros of the vector-valued
    # function instead of the sum of its squares.

    if use_scipy_minimize:

        SciMin = scipy.optimize.minimize(minfunc, iv, jac=jac_minfunc, method='BFGS', options={'return_all':True})

        print()
        print()

        if SciMin.success:
            for pair in zip(coeff_vars, SciMin.x): print( pair)
        else:
            print( SciMin.message)

    elif use_scipy_root:

        # optimize.root methods:
        # 'hybr' (the default) requires same num of eqns as vars
        # 'lm' uses a QR factorization of the Jacobian, then the Levenberg–Marquardt line search algorithm
        # the others uses various approximations to the Jacobian

        if use_divExp:
            SciMin = scipy.optimize.root(fns_divExpA, iv, jac=jac_fns_divExpA, method='lm')
        else:
            SciMin = scipy.optimize.root(fns_divSqrtA, iv, jac=jac_fns_divSqrtA, method='lm')

        print()
        print()

        if SciMin.success:
            printv(SciMin.x)
        else:
            print( SciMin.message)

    else:

        # Custom root-finding algorithm, based on Numerical Recipes
        #
        # - LU factorization of the Jacobian, computed in parallel
        # - exact-fit line search algorithm

        i = 0
        gtol = 1e-5
        gnorm = np.linalg.norm(jac_minfunc(iv))

        while (not limit or i < limit) and gnorm > gtol:
            #minfunc(iv)
            iv = optimize_step(iv)
            gnorm = np.linalg.norm(jac_minfunc(iv))
            i += 1

        global final_iv
        final_iv = iv
        for pair in zip(coeff_vars, iv): print(pair)

def random_numerical_ten(limit=10):
    for i in range(limit):
        print("Random seed", i)
        random_numerical(i)


def find_relation():

    # search for integer relations among the approximate solutions

    #ints = [ZZ(round(v)) for v in 2^26 / sqrt(sum(SciMin.x * SciMin.x)) * SciMin.x]
    #norm = 2^26 / sqrt(sum(SciMin.x * SciMin.x))
    norm = 1/min(abs(SciMin.x))
    ints = [ZZ(round(v)) for v in  norm * SciMin.x] + [ZZ(round(norm))]

    print(ints)

    L = matrix(list(matrix.identity(len(ints))) + [tuple(ints)]).transpose().LLL()

    print(L)

    for Lrow in L:

        rel = matrix(BWB.gens() + (1,)) * Lrow[0:-1]

        print(rel)

        V = Lrow[0:-1]

        # len(V)-1 so as to drop the "1" term at the end
        for i in range(len(V)-1):
            if V[i] == 1:
                print(i,Lrow)
                V[i] = 0
                ints[i] = - (matrix(V) * matrix(ints).transpose())[0,0]
                break

        if Lrow[-1] != 0:
            break

def convert_matrix(index):
    dump_fn = 'csr0-' + str(index)
    collector = CollectorClass()
    collector.load_from_glob('newpickle[0246]*', index, int(16))
    f = open(dump_fn, 'wb')
    pickle.dump(collector, f)
    f.close()
    print("Dumped to", dump_fn)


def do_work():
    m = BaseManager(address=('c200-1.fios-router.home', 50000))
    m.connect()
    while True:
        work = m.get_work()._getvalue()
        print("Working on", work)
        timefunc(convert_matrix, work)

def prep_work():
    m = BaseManager(address=('c200-1.fios-router.home', 50000))
    m.connect()
    m.add_work(list(range(16)))

# If we're running on my development laptop, initialize a simple calculation when this file is loaded

import platform

if platform.node() == 'samsung' and not 'no_init' in globals():

    prep_hydrogen()
    multi_init()
    multi_expand()

if platform.node().startswith('c200') and not 'no_init' in globals() and not 'mc' in globals():
    timefunc(prep_helium)
    timefunc(create_polynomial_ring)

if platform.node() == 'c200-1' and not 'no_init' in globals() and not 'mc' in globals():
    start_manager_process()
    for i in range(12):
        start_collector()
    ccs = mc.get_collectors()
    for i in range(12):
        ccs[i].load_matrix('matrix-' + str(i) + '.pickle')
    for i in range(12):
        ccs[i].join_threads()

def load_v(fn):
    f = open(fn, 'rb')
    up = pickle.Unpickler(f)
    v = up.load()
    f.close()
    return v
