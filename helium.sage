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
# gradient descent algorithm in the hopes of finding an approximate
# point on that variety.  If successful, we look for algebraic
# relationships between the point's coordinates, hoping to find the
# algebraic equations that define the irreducible component that the
# point lies on.
#
# A big computational bottleneck is the space complexity of expanding
# out the differential equation.  To solve this problem, we use
# multiple processes to expand out the polynomial incrementally.
#
# by Brent Baccala
#
# first version - August 2019
# latest version - September 2019
#
# no rights reserved; you may freely copy, modify, or distribute this
# program
#
# TODO list:
# - collector should parse directly into matrix form
# - compile a better regular expression to parse terms
# - more easily turn multiprocessing on and off
# - better naming of functions and variables like bwb4
# - allow worker processes on different hosts
# - remove unused code like Rosenfeld-Groebner
# - allow second-order ODEs in trial form
# - check collected coefficient polynomials to see if they factor
# - compute bound on degree of coefficient monomials
# - automate finding polynomial relations
# - save checkpoints of optimization iterations
# - optimize scipy sparse matrix by vector multiplication
# - optimize creation of multi-vectors
# - allow workers to be added or removed on the fly

from itertools import *

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

def trial_polynomial(base, cvars, rvars, degree):
    # cterms are coefficient terms
    # rterms are radius terms
    cterms = flatten([combinations_with_replacement(cvars, d) for d in range(degree+1)])
    # use this 'terms' for real
    terms = list(map(mul, (product(map(mul, cterms), map(mul, powerset(rvars))))))
    # use this 'terms' for testing
    # terms = list(map(mul,cterms)) + list(rvars)
    coefficients = tuple(var(base+str(c)) for c in range(len(terms)))
    poly = sum([var(base+str(c))*v for c,v in enumerate(terms)])
    return (coefficients, poly)

var('E')

def Del(Psi,vars):
    return sum([diff(Psi,v,2) for v in vars])

def prep_hydrogen():
    global eq, H, Psi, coeff_vars, A, B, Avars, Bvars, cvars, rvars

    cvars = (x1,y1,z1)
    rvars = (r1,)
    (Avars, A) = trial_polynomial('a', cvars, rvars, 1)
    (Bvars, B) = trial_polynomial('b', cvars, rvars, 1)

    coeff_vars = (E,) + Avars + Bvars

    global Phi
    Phi = function('Phi')(*cvars)

    Psi = A*Phi

    def H(Psi):
        return - 1/2 * Del(Psi,[x1,y1,z1]) - (1/r1)*Psi

    eq = H(Psi) - E*Psi

    # Phi = e^B
    # diff(Phi,B) = Phi
    # diff(Phi,t) = diff(Phi,B) * diff(B,t)
    # diff(Phi,t) = diff(B,t) * e^B = diff(B,t) * Phi
    # diff(Phi,t,2) = diff(B,t,2) * Phi + diff(B,t) * diff(Phi,t) = diff(B,t,2) * Phi + diff(B,t)^2 * Phi

    # eq = eq.subs({diff(Phi, x1):    diff(B, x1) * Phi,
    #               diff(Phi, x1, 2): (diff(B,x1,2) + diff(B,x1)^2) * Phi,
    #               diff(Phi, y1):    diff(B, y1) * Phi,
    #               diff(Phi, y1, 2): (diff(B,y1,2) + diff(B,y1)^2) * Phi,
    #               diff(Phi, z1):    diff(B, z1) * Phi,
    #               diff(Phi, z1, 2): (diff(B,z1,2) + diff(B,z1)^2) * Phi})

    dict1 = {diff(Phi,v): diff(B,v)*Phi for v in cvars}
    dict2 = {diff(Phi,v,2): diff(dict1[diff(Phi,v)],v) for v in cvars}

    # replace Phi(x1,y1,z1) with Phi to reduce ginac's memory utilization
    eq = eq.subs(dict2).subs(dict1).subs({Phi: SR.var('Phi')})
    Phi = SR.var('Phi')

def prep_helium():
    global eq, H, Psi, coeff_vars, A, B, Avars, Bvars, cvars, rvars

    cvars = (x1,y1,z1, x2,y2,z2)
    rvars = (r1,r2,r12)
    (Avars, A) = trial_polynomial('a', cvars, rvars, 1)
    (Bvars, B) = trial_polynomial('b', cvars, rvars, 1)

    coeff_vars = (E,) + Avars + Bvars

    global Phi
    Phi = function('Phi')(*cvars)

    Psi = A*Phi

    def H(Psi):
        return - 1/2 * Del(Psi,[x1,y1,z1]) - 1/2 * Del(Psi,[x2,y2,z2]) - (2/r1)*Psi - (2/r2)*Psi + (1/r12)*Psi

    eq = H(Psi) - E*Psi

    dict1 = {diff(Phi,v): diff(B,v)*Phi for v in cvars}
    dict2 = {diff(Phi,v,2): diff(dict1[diff(Phi,v)],v) for v in cvars}

    # replace Phi(x1,y1,z1) with Phi to reduce ginac's memory utilization
    eq = eq.subs(dict2).subs(dict1).subs({Phi: SR.var('Phi')})
    Phi = SR.var('Phi')

# Now we want to replace all of the sqrt(...) factors with 'r',
# and we use a clever Python trick to build a dictionary
# that maps expressions to variable names.

def varName(var):
    for name,value in globals().items():
        if id(var) == id(value):
            return name
    return None

def mk_maps(rvars):
    return {v.operands()[0] : SR.var(varName(v)) for v in rvars}

# convert all (x^2+y^2+z^2)^(n/2) expressions to r^n
def bwb(expr):
    if isinstance(expr, Expression) and expr.operator():
       if expr.operator() == operator.pow and bool(expr.operands()[0] in maps):
           return maps[expr.operands()[0]]^(expr.operands()[1] * 2)
           #num = int(expr.operands()[1] * 2)
           #return ((expr.operands()[0]^(int(num/2))) * (maps[expr.operands()[0]]^(int(num%2))))
       else:
           return expr.operator()(*map(bwb, expr.operands()))
    else:
       return expr

def create_bwb4():
    global maps, bwb4, bwb4a
    # first, build the dictionary that maps expressions like (x1^2+y1^2+z1^2) to variables like r1
    maps = mk_maps(rvars)
    # next, convert all of the roots in the equation to use the r-variables
    bwb4a = bwb(eq)
    # Find the least common denominator of all of the terms, then
    # clear the denominators and expand out all of the powers.
    # This is faster than expand(bwb4a.numerator()).  One test
    # on helium ran in 61 sec, vs 337 sec for expand/numerator.
    lcm_denominator = lcm(map(denominator, bwb4a.operands()))
    bwb4 = expand(bwb4a*lcm_denominator)
    # bwb4 is now a polynomial, but it's got higher powers of r's in it
    # assert bwb4.numerator() == bwb4



# Next... convert powers of r's to x,y,z's, expand out powers, and
# collect like x,y,z's terms together to get a system of polynomials
#
# This is a slow step, so I've tried several different ways to do it.

# Runs at a reasonable speed after create_bwb4() has been called.
# 'bwb4a' still needs to be expanded.

def SR_expand2a():
    global bwb4a
    sdict = {SR.var(v)^d : (globals()[v]^d, SR.var(v)*globals()[v]^(d-1))[d%2] for d in range(2,8) for v in ('r1','r2','r12')}
    bwb4a = bwb4.subs(sdict)

# Perform the expansion step using multiple processes.
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
# create_bwb4()
# SR_expand2a()
# ops = bwb4a.operands()
# start_manager_process()
# mc.set_range(len(ops), 100)
# mc.start_worker()
# wc = mc.getq()[0]  (or wc = mc.get_workers()[0])
# wc.get_data()

# number of operands to process in each expander worker
blocksize = 100

# number of collector processes
num_collectors = 2

# number of simultaneous expander processes
num_expanders = 2

import timeit

def extract_ops():
    global ops
    ops = bwb4a.operands()

def multi_init():
    prep_hydrogen()
    t = timeit.timeit(lambda: create_bwb4(), number=1)
    print 'create_bwb4() : %s sec'%(t)
    t = timeit.timeit(lambda: SR_expand2a(), number=1)
    print 'SR_expand2a() : %s sec'%(t)
    t = timeit.timeit(lambda: extract_ops(), number=1)
    print 'extract_ops() : %s sec'%(t)
    start_manager_process()
    mc.set_range(len(ops), blocksize)

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

def start_collector():
    collector_manager = BaseManager()
    collector_manager.start()
    cc = collector_manager.CollectorClass()
    mc.register_collector(cc)
    managers[cc._token] = collector_manager
    return cc

def start_worker():
    worker_manager = BaseManager()
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
    fetched_results = map(lambda x: x.get(), results)
    sizes = [cc.nrows() for cc in ccs]
    value_to_count = dict(zip(*np.unique(np.hstack(fetched_results), return_counts=True)))
    dups = set((u for u,c in value_to_count.items() if c>1))
    generators = [((j,e) for j,e in enumerate(fetched_results[i]) if e in dups) for i in range(len(ccs))]
    indices = [[] for i in range(len(ccs))]
    remaining_iterations = sum(sizes) - len(value_to_count)
    while len(dups) > 0:
        if remaining_iterations % 1000 == 0:
            print remaining_iterations, 'iterations remaining'
        maxsize = max(sizes)
        try:
            i = (i for i,c in enumerate(sizes) if c == maxsize and generators[i] is not None).next()
        except StopIteration:
            break
        try:
            (index, value) = generators[i].next()
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

import multiprocessing, logging
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

from rfoo.utils import rconsole

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
        return None

    def autocreate(self, classname, *args):
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        (ident, exposed) = server.create(None, classname, *args)
        token = multiprocessing.managers.Token(typeid=classname, address=server.address, id=ident)
        proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
        server.decref(None, ident)
        return proxy

    # convenience functions for development
    def getpid(self):
        return os.getpid()
    def load(self, filename):
        load(filename)
    def rconsole(self, port=54321):
        rconsole.spawn_server(port=port)

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

    def register_manager(self, mc):
        self.mc = mc
    def register_collectors(self, collectors):
        self.collectors = collectors
        self.dicts = map(dict, [[]] * len(collectors))
    @async_method
    def shutdown(self):
        for i in range(len(self.collectors)):
            del self.collectors[-1]
    @async_method
    def start_expand(self, start, stop):
        # sleep for a fraction of a second here to let this method
        # return back to the manager before we begin work
        time.sleep(0.1)
        expr = expand(sum(islice(ops, start, stop)))
        for monomial in expr.operands():
            s = str(monomial)
            if s[0] is '-':
                sign='-'
                s=s[1:]
            else:
                sign='+'
            vs = s.split('*')
            key = '*'.join(ifilter(lambda x: any([x.startswith(c) for c in ('x', 'y', 'z', 'r', 'P')]), vs))
            value = '*'.join(ifilterfalse(lambda x: any([x.startswith(c) for c in ('x', 'y', 'z', 'r', 'P')]), vs))
            dictnum = hash(key) % len(self.collectors)
            self.dicts[dictnum][key] = self.dicts[dictnum].get(key, '') + sign + value
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
    """

    def __init__(self, matrix):
        self.M = matrix
        (self.rows, cols) = matrix.shape
        # XXX might be a bit slow
        self.scaling = np.amax(abs(matrix), axis=1)
        self.U = np.zeros((cols, cols))

    def get(self):
        self.join_threads()
        return self.M

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
        row = (abs_b / self.scaling).argmax()
        return (abs_b[row], row, 0)

    def fetch_and_remove_row(self, row):
        res = self.M[row]
        self.M = np.delete(self.M, row, axis=0)
        self.scaling = np.delete(self.scaling, row)
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
    def _start_calculation(self, vec):
        M = np.stack([self.collector.dot(self.collector.generate_multi_D_vector(vec, var)) for var in coeff_vars], axis=1)
        LUMatrix.__init__(self, M)

    def __init__(self, collector, vec):
        self.collector = collector
        self._start_calculation(vec)

class CollectorClass(Autoself):
    result = {}
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
        fp = open(fn, 'w')
        pickle.dump(self.M, fp)
        fp.close()
        logger.info('matrix dumped to %s', fn)
    @async_method
    def load_matrix(self, fn):
        fp = open(fn)
        if self.M is None:
            self.M = pickle.load(fp)
        else:
            self.M = sp_unique(sp.vstack((self.M, pickle.load(fp))))
        fp.close()
        logger.info('matrix loaded from %s', fn)
        logger.info(repr(self.M))

    def times(self):
        return (self.dot_calls, self.dot_times, self.multi_vec_calls, self.multi_vec_times, self.multi_D_vec_calls, self.multi_D_vec_times)
    def nrows(self):
        return self.M.shape[0]
    def len_result(self):
        return len(self.result)
    def len_values(self):
        return sum([len(v) for v in a.result.values()])
    def len_keys(self):
        return sum([len(v) for v in a.result.keys()])
    def nterms(self):
        return sum([v.count('+')+v.count('-') for v in self.result.values()])

    # Now we want to evaluate the sum of squares of polynomials p, as
    # well as the first derivatives of the sum of squares, which is
    # the sum of (2 p dp/dv).  Using Sage symbolic expressions is slow
    # and consumes much memory, so this version of the code uses numpy
    # matrices instead.

    # Given a vector, generate a multi-vector of the vector itself,
    # then all the biproducts of the vector elements, then all the
    # triproducts of the vector elements.
    #
    # generate_multi_vector([a,b,c])
    #   -> [a,b,c,a^2,a*b,a*c,b^2,b*c,c^2,a^3,a^2*b,a^2*c,a*b^2,a*b*c,a*c^2,a*b^2,b^2*c,b*c^2,c^3]

    def generate_multi_vector(self, v):
        start_time = time.time()
        npv = np.array(v)
        two_products = [npv[i] * npv[i:] for i in range(len(v))]
        three_products = [np.hstack(two_products[i:]) * npv[i] for i in range(len(coeff_vars))]
        res = np.hstack([npv] + two_products + three_products)
        self.multi_vec_times += time.time() - start_time
        self.multi_vec_calls += 1
        return res

    def term_to_vector(self, term):
        if term[0] == '-':
            sign = -1
            term = term[1:]
        else:
            sign = 1
            if term[0] == '+':
                term = term[1:]
        (coeff, monomial) = re.split('\*', term, 1)
        if not re.match('^[0-9]*$', coeff):
            monomial = coeff + '*' + monomial
            coeff = '1'
        coeff = sign * int(coeff)
        index = self.indices[monomial]
        return coeff * self.monomial_vectors[monomial]

    @async_method
    def convert_to_matrix(self):
        self.i = 0
        self.indices = {str(pair[1]) : pair[0] for pair in enumerate(self.generate_multi_vector(coeff_vars))}
        veclen = len(self.indices)
        self.dok = scipy.sparse.dok_matrix((len(self.result), veclen), np.int64)
        for value in self.result.values():
            terms = re.split('([+-][^+-]+)', value)
            for term in ifilter(bool, terms):
                if term[0] == '-':
                    sign = -1
                    term = term[1:]
                else:
                    sign = 1
                    if term[0] == '+':
                        term = term[1:]
                try:
                    (coeff, monomial) = re.split('\*', term, 1)
                except ValueError:
                    coeff = '1'
                    monomial = term
                if not re.match('^[0-9]*$', coeff):
                    monomial = coeff + '*' + monomial
                    coeff = '1'
                coeff = sign * int(coeff)
                index = self.indices[monomial]
                self.dok[self.i, index] += coeff
            self.i += 1

        #self.M = np.unique(self.M, axis=0)
        #self.M = self.dok.tocsr()
        self.M = sp_unique(self.dok, axis=0, new_format='csr')
        logger.info('convert_to_matrix done')

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

    def generate_multi_D_vector(self, v, var):
        start_time = time.time()
        ind = coeff_vars.index(var)

        npv = np.array(v)

        firsts = np.zeros(npv.size)
        firsts[ind] = int(1)

        two_products = [npv[i] * npv[i:] for i in range(npv.size)]
        two_products_D = [firsts[i] * npv[i:] + firsts[i:] * npv[i] for i in range(npv.size)]

        three_products = [np.hstack(two_products[i:]) * firsts[i] + np.hstack(two_products_D[i:]) * npv[i] for i in range(npv.size)]

        res = np.hstack([firsts] + two_products_D + three_products)
        self.multi_D_vec_times += time.time() - start_time
        self.multi_D_vec_calls += 1
        return res

    def verify_D_vector(self):
        all([all([bool(diff(e,v)==d) for e,d in zip(generate_multi_vector(coeff_vars), generate_multi_D_vector(coeff_vars, v))]) for v in coeff_vars])

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

    def jacobian_fns(self, vec):
        r"""
        Evaluate the Jacobian matrix (the matrix of first-order
        partial derivatives) of our polynomials

        INPUT:

        - ``vec`` -- a vector of real values for all coeff_vars

        OUTPUT:

        - the Jacobian matrix, evaluated at ``vec``, as a numpy
        matrix wrapped in a proxy object
        """
        return self.autocreate('JacobianMatrix', self, vec)

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
    def jacobian_sum_of_squares(self, vec):
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

from multiprocessing.managers import BaseManager
BaseManager.register('ManagerClass', ManagerClass)
BaseManager.register('ExpanderClass', ExpanderClass)
BaseManager.register('CollectorClass', CollectorClass)
BaseManager.register('AsyncResult', AsyncResult)
BaseManager.register('JacobianMatrix', JacobianMatrix)

def start_manager_process():
    global manager, mc
    manager = BaseManager()
    manager.start()
    mc = manager.ManagerClass()


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

# The function we're trying to minimize: the sum of squares of the polys
# that define the solution variety, divided by two factors we're trying
# to avoid: the norm of 'v' (avoid the origin), and zero_variety (which
# includes the origin).

last_time = 0

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

    res = np.hstack(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs]))
    return res

def fndivA(v):
    # Save a copy of vector to aid in stopping and restarting the calculation
    global last_v
    last_v = v

    res = np.hstack(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs]))
    Adenom = sqrt(sum([square(v[i]) for i in Aindices]))

    global last_time
    sum_of_squares = sum(square(res/Adenom))
    if last_time == 0:
        print sum_of_squares
    else:
        print "{:<30} {:20} sec".format(sum_of_squares, time.time()-last_time)
    last_time = time.time()
    return res/Adenom

def jacfn(v):
    r"""
    Evaluate the Jacobian matrix (the matrix of first-order
    partial derivatives) of our polynomials

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - the Jacobian matrix, evaluated at ``vec``, as a numpy matrix,
    with as many rows as polynomials and as many columns as coeff_vars

    ALGORITHM:

    Call the `jacobian_fns` method for all CollectionClass's in
    parallel, then concatenate all of the results together
    (and transpose them).
    """

    res = np.vstack(map(lambda x: x.get(), [cc.jacobian_fns(v) for cc in ccs]))
    return res

def jac_fndivA(v):
    global N,dN,Av,Adenom
    N = np.hstack(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs]))
    dN = np.vstack(map(lambda x: x.get(), [cc.jacobian_fns(v) for cc in ccs]))
    Av = v * np.array([c in Avars for c in coeff_vars])   # could form a global vector for this
    Adenomsq = sum([square(v[i]) for i in Aindices])
    Adenom = sqrt(Adenomsq)
    res = dN/Adenom - np.outer(N,Av)/Adenom/Adenomsq
    return res

def sum_of_squares(v):
    return sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))

def minfunc(v):
    # Save a copy of vector to aid in stopping and restarting the calculation
    global last_v
    last_v = v

    d = dict(zip(coeff_vars, v))
    sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))
    res = real_type(sum_of_squares / zero_variety.subs(d))
    global last_time
    if last_time == 0:
        print res
    else:
        print "{:<30} {:20} sec".format(res, time.time()-last_time)
    last_time = time.time()
    return res

def jac(v):
    global last_v
    last_v = v

    d = dict(zip(coeff_vars, v))
    sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))
    jacobian_sum_of_squares = sum(map(lambda x: np.array(x.get()), [cc.jacobian_sum_of_squares(v) for cc in ccs]))
    zero_var = v.dtype.type(zero_variety.subs(d))
    Av = v * np.array([c in Avars for c in coeff_vars])
    res = ((jacobian_sum_of_squares*zero_var - 2*np.array(Av)*sum_of_squares)/zero_var^2)
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
    jacobians = [cc.jacobian_fns(v) for cc in ccs]
    (L, U, f) = LU_decomposition(jacobians)
    direction = scipy.linalg.lu_solve((lu, piv), f)

    # pick a step size in the `direction`
    v0 = minfunc(vec)
    norm = sum(square(direction))
    evalstep = direction*v0/norm

    # We want to sample at seven points to fit a sixth degree polynomial.
    # We expect a zero "close" to -1, so this will sample three points
    # on either size of it
    points = [vec + i*evalstep for i in [-4,-3,-2,-1,0,1,2]]
    values = map(sum_of_squares, points)

    # Now fit a polynomial to this data
    N = np.polynomial.polynomial.Polynomial.fit([-4,-3,-2,-1,0,1,2],values,6,[])

    # The denominator is the sum of (A_0 + lambda A_d)^2 for all As
    D = sum([square(np.polynomial.polynomial.Polynomial((vec[i], evalstep[i]))) for i in Aindices])

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
    return nextstep

def random_numerical(iv=0):

    import scipy.optimize

    # eqns.append(E+1/4)
    # eqns.append(BWB.gen(1) - 1)

    nvars = len(coeff_vars)

    # even though we're using numpy, we don't need to set its PRNG
    # seed, (which would require calling numpy.random.seed()), since
    # the algorithm is deterministic after the iv is picked

    if isinstance(iv, int) or isinstance(iv, Integer):
        random.seed(iv)        # for random
        set_random_seed(iv)    # for RR.random_element()
        iv = np.array([random.random() for i in range(nvars)])

    # We know the zero variety (all Avar's zero, so Psi is zero) will be
    # a "solution", but we want to avoid it

    global zero_variety, Aindices
    zero_variety = sum(map(square, Avars))
    Aindices = [i for i,c in enumerate(coeff_vars) if c in Avars]

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

    # SciMin = scipy.optimize.minimize(minfunc, iv, jac=jac, method='BFGS', options={'return_all':True})

    # optimize.root methods:
    # 'hybr' (the default) requires same num of eqns as vars
    # 'lm' uses a QR factorization of the Jacobian, then the Levenberg–Marquardt line search algorithm
    # the others uses various approximations to the Jacobian

    SciMin = scipy.optimize.root(fndivA, iv, jac=jac_fndivA, method='lm')

    # Our root-finding algorithm:
    #
    # - LU factorization of the Jacobian, computed in parallel
    # - exact-fit line search algorithm

    #i = 0
    #while i < 50:
    #    iv = optimize_step(iv)
    #    i += 1

    print
    print

    if SciMin.success:
        for pair in zip(coeff_vars, SciMin.x): print pair
    else:
        print SciMin.message

def random_numerical_ten(limit=10):
    for i in range(limit):
        print "Random seed", i
        random_numerical(i)


def find_relation():

    # search for integer relations among the approximate solutions

    #ints = [ZZ(round(v)) for v in 2^26 / sqrt(sum(SciMin.x * SciMin.x)) * SciMin.x]
    #norm = 2^26 / sqrt(sum(SciMin.x * SciMin.x))
    norm = 1/min(abs(SciMin.x))
    ints = [ZZ(round(v)) for v in  norm * SciMin.x] + [ZZ(round(norm))]

    print ints

    L = matrix(list(matrix.identity(len(ints))) + [tuple(ints)]).transpose().LLL()

    print L

    for Lrow in L:

        rel = matrix(BWB.gens() + (1,)) * Lrow[0:-1]

        print rel

        V = Lrow[0:-1]

        # len(V)-1 so as to drop the "1" term at the end
        for i in range(len(V)-1):
            if V[i] == 1:
                print i,Lrow
	        V[i] = 0
	        ints[i] = - (matrix(V) * matrix(ints).transpose())[0,0]
	        break

        if Lrow[-1] != 0:
            break
