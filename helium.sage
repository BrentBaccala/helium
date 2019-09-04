#
# Python code to search for solutions of Hydrogen and Helium
# non-relativistic time-invariant Schrodinger equation
#
# by Brent Baccala
# August 2019

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

# Rosenfeld-Groebner runs very slowly, even on Hydrogen first excited state
#
# I don't use it anymore.

def rosenfeld_groebner():
    # If I specify derivations=[x,y,z], then I have to use all three
    # variables when constructing A(x,y,z), i.e, A(x,y) produces an error

    (x,y,z) = var('x,y,z')
    (A,B,expB,r,Psi) = function('A,B,expB,r,Psi')
    (a,b,c,d,e) = var('a,b,c,d,e')
    (f,g,h,i,j) = var('f,g,h,i,j')
    var('E')

    parameters = [a,b,c,d,e,  f,g,h,i,j,  E]

    from sage.calculus.DifferentialAlgebra import DifferentialRing, BaseFieldExtension
    DR = DifferentialRing(derivations = [x,y,z],
                          blocks = [[A,B,expB,Psi,r], parameters],
                          parameters = parameters)

    def H(Psi):
       return -diff(Psi,x,2)-diff(Psi,y,2)-diff(Psi,z,2)-(1/r(x,y,z))*Psi

    # timings for various trial functions:
    # all used [[A,B,expB,Psi,r], parameters] for block ordering
    #
    #  5 sec: Psi = e * exp(i*r)
    # 45 sec: Psi = (a*x + e) * exp(i*r)
    # 52 sec: Psi = (a*x + b*y + e) * exp(i*r)
    # >13 hr: Psi = (d*r + e) * exp(i*r)

    #rels = [A(x,y,z) == a*x + b*y + c*z + d*r(x,y,z) + e,
    #        B(x,y,z) == f*x + g*y + h*z + i*r(x,y,z) + j,
    #rels = [A(x,y,z) == d*r(x,y,z) + e,
    rels = [A(x,y,z) == e,
            B(x,y,z) == i*r(x,y,z),
            r(x,y,z)^2 == x^2 + y^2 + z^2,
            Psi(x,y,z) != 0,
            Psi(x,y,z) == A(x,y,z)*expB(x,y,z),
            diff(expB(x,y,z),x) == diff(B(x,y,z), x) * expB(x,y,z),
            diff(expB(x,y,z),y) == diff(B(x,y,z), y) * expB(x,y,z),
            diff(expB(x,y,z),z) == diff(B(x,y,z), z) * expB(x,y,z),
            H(Psi(x,y,z)) == E*Psi(x,y,z)]

    # expanding H and dividing out expB helps (a little):
    #
    #  8 sec: Psi = (a*x + e) * exp(i*r)
    #  8 sec: Psi = (a*x + b*y + e) * exp(i*r)
    # >1 min: Psi = (d*r + e) * exp(i*r)

    Psi = A(x,y,z)*exp(B(x,y,z))
    rels = [A(x,y,z) == (e),
            B(x,y,z) == i*r(x,y,z),
            r(x,y,z)^2 == x^2 + y^2 + z^2,
            expand(H(Psi)/Psi - E)]

    RDC = DR.RosenfeldGroebner(rels)

    #Field = BaseFieldExtension (generators = parameters)
    #RDC = DR.RosenfeldGroebner (rels, basefield = Field)

    #print RDC

    for id in RDC: print id



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

# This version (SR_expand) expands in the symbolic ring, by first converting
# higher powers of the rvars to their expansions, then recursively extracting
# coefficients from all of the cvars and rvars.

def bwb2(expr):
    if isinstance(expr, Expression) and expr.operator():
       if expr.operator() == operator.pow and bool(expr.operands()[0] in SRr_s):
           var = expr.operands()[0]
           num = int(expr.operands()[1])
           #return ((x1^2+y1^2+z1^2)^(int(num/2))) * (expr.operands()[0]^(int(num%2)))
           return (globals()[str(var)])^(2*int(num/2)) * (var^(int(num%2)))
       else:
           return expr.operator()(*map(bwb2, expr.operands()))
    else:
       return expr

# This version expands everything using Gynac (i.e, the symbolic ring), and collects
# all the coefficients together in a list.  It runs very slow.

def SR_expander(expr, vars):
    if isinstance(expr, Expression) and len(vars) > 0:
        l = map(lambda x: SR_expander(x, vars[1:]), expr.coefficients(vars[0], sparse=False))
        return list(flatten(l))
    else:
        return [expr]

# This version collects coefficients together into a dictionary.  It is also very slow.

def SRdict_expander(expr, vars_expanded, vars_to_expand):
    if isinstance(expr, Expression) and len(vars_to_expand) > 0:
        v = vars_to_expand[0]
	map(lambda l: SRdict_expander(l[0], vars_expanded * v^l[1], vars_to_expand[1:]), expr.coefficients(v))
    else:
        SRdict[vars_expanded] = SRdict.get(vars_expanded, 0) + expr

# This version converts each monomial into a string and then splits
# the strings apart using Python operations.  It's pretty slow, but
# is more amenable to parallelisation.

def SRdict_expander2(expr):
    assert expr.operator() is sage.symbolic.operators.add_vararg
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
        SRdict[key] = SRdict.get(key, '') + sign + value

# this version returns a dictionary instead of using a global

def SRdict_expander2a(expr):
    assert expr.operator() is sage.symbolic.operators.add_vararg
    SRdict = dict()
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
        SRdict[key] = SRdict.get(key, '') + sign + value
    return SRdict

def SR_expand():
    global eqns
    bwb4 = create_bwb4()
    eqns = set(SR_expander(bwb2(bwb4), cvars + SRr_s + (Phi,)))

# Runs at a reasonable speed after create_bwb4() has been called.
# 'bwb4a' still needs to be expanded.

def SR_expand2a():
    global bwb4a
    sdict = {SR.var(v)^d : (globals()[v]^d, SR.var(v)*globals()[v]^(d-1))[d%2] for d in range(2,8) for v in ('r1','r2','r12')}
    bwb4a = bwb4.subs(sdict)

# Runs out of memory.

def SR_expand2b():
    global bwb4b
    bwb4b = expand(bwb4a)

def SRexpand2c(start, stop):
    ops = bwb4a.operands()
    ex = []
    ex.append(expand(sum(islice(ops, start, stop))))

def SRdict_expander4():
    global thousands
    for thousands in range(0, len(ops), 1000):
        SRdict_expander2(expand(sum(islice(ops, thousands, thousands+1000))))

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
        start_collector().load_json(fn)
    global ccs
    ccs = mc.get_collectors()

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
# actual result is obtained by calling the `get` method on the proxy.

class AsyncResult:
    def __init__(self, outerself, method, *args):
        self._cond = threading.Condition()
        self._ready = False

        def target(outerself, method, *args):
            self._result = method(outerself, *args)
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

    def get(self):
        self._cond.acquire()
        if not self._ready:
            self._cond.wait()
        self._cond.release()
        return self._result

def async_result(method):
    def inner(self, *args):
        classname = 'AsyncResult'
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        (ident, exposed) = server.create(None, classname, *((self, method) + args))
        token = multiprocessing.managers.Token(typeid=classname, address=server.address, id=ident)
        proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
        server.decref(None, ident)
        return proxy

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

class Autoself:
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

    def join_threads(self):
        if hasattr(self, 'threads'):
            for th in self.threads:
                logger.debug('join %s', th)
                th.join()

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

class CollectorClass(Autoself):
    result = {}

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
    def dump_matrix(self):
        fn = '/tmp/' + str(os.getpid()) + '.pickle'
        fp = open(fn, 'w')
        pickle.dump(self.M, fp)
        fp.close()
        logger.info('matrix dumped to %s', fn)

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

    @staticmethod
    def generate_multi_vector(v):
        two_products = [map(mul, izip(v[i:], repeat(e))) for i,e in enumerate(v)]
        three_products = [map(mul, izip(flatten(two_products[i:]), repeat(e))) for i,e in enumerate(v)]
        return np.array(list(flatten([v, flatten(two_products), flatten(three_products)])))

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

    # no longer works because I can't figure how to multiply numpy matrices with Sage Expressions in them
    def verify_matrix(self, vars):
        vec = self.generate_multi_vector(vars)
        set1 = set([eval(preparse(result)) for result in self.result])
        set2 = set(list(self.M.dot(vec)))
        return set1 == set2

    def generate_multi_D_vector(self, v, var):
        ind = coeff_vars.index(var)
        firsts = [int(0)] * len(coeff_vars)
        firsts[ind] = int(1)

        two_products = [map(mul, izip(v[i:], repeat(e))) for i,e in enumerate(v)]
        two_products_D = [map(sum, zip(*(map(mul, izip(v[i:], repeat(firsts[i]))), map(mul, izip(firsts[i:], repeat(v[i])))))) for i in range(len(coeff_vars))]

        three_products = [map(sum, zip(*(map(mul, izip(flatten(two_products[i:]), repeat(firsts[i]))), map(mul, izip(flatten(two_products_D[i:]), repeat(v[i])))))) for i in range(len(coeff_vars))]

        return np.array(list(flatten([firsts, flatten(two_products_D), flatten(three_products)])))

    def verify_D_vector(self):
        all([all([bool(diff(e,v)==d) for e,d in zip(generate_multi_vector(coeff_vars), generate_multi_D_vector(coeff_vars, v))]) for v in coeff_vars])

    @async_result
    def sum_of_squares(self, vec):
        multivec = self.generate_multi_vector(vec)
        return sum(square(self.M.dot(multivec)))
    @async_result
    def D_sum_of_squares(self, vec, var):
        # d(p^a s^b t^c)/ds = b(p^a s^(b-1) t^c),
        # so (p^a s^b t^c) should map to b(p^a s^(b-1) t^c)
        # dp/ds = 0         ds/ds = 1
        # d(s^2)/ds = 2s    d(ps)/ds = p     d(pt)/ds = 0
        # d(s^3)/ds = 3s^2  d(ps^2)/ds = 2ps   d(pst) = pt
        multivec = self.generate_multi_vector(vec)
        multiDvec = self.generate_multi_D_vector(vec, var)
        return sum(2 * self.M.dot(multivec) * self.M.dot(multiDvec))


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

def start_manager_process():
    global manager, mc
    manager = BaseManager()
    manager.start()
    mc = manager.ManagerClass()





def PolynomialRing_expand():

    # This method creates a PolynomialRing QQ(E,As,Bs)[cvars, rvars],
    # then quotients it by the ideal generated by the relationships
    # between the rvars and the cvars, converts bwb4 into this new
    # ring, extracts all of the QQ(E,As,Bs) coefficients, and
    # forms the numerators into eqns.
    #
    # BWB : QQ(E,As,Bs)
    # BWB2: QQ(E,As,Bs)[cvars, rvars]
    # BWB3: QQ(E,As,Bs)[cvars, rvars]/(relations)

    global BWB, BWB2, BWB3

    # use custom term ordering to prioritize elimination of 'E' variable
    # if Sage's Groebner basis-based techniques are used
    order = TermOrder('deglex(1),degrevlex({})'.format(len(Avars)+len(Bvars)))
    BWB = PolynomialRing(QQ, (var('E'),) + Avars + Bvars, order=order)
    #BWB = PolynomialRing(QQ, (var('E'),) + Avars + Bvars)

    #BWB2.<x,y,z,r> = Frac(BWB)[]
    BWB2 = PolynomialRing(Frac(BWB), cvars + tuple(maps.get(v^2, v) for v in rvars))
    #BWB3 = BWB2.quo(r^2-(x^2+y^2+z^2))
    BWB3 = BWB2.quo([k - v^2 for k,v in maps.items()])

    bwb4 = create_bwb4()

    global bwb3
    global eqns
    #bwb3 = BWB3(numerator(expand(bwb(eq/exp(B)))))
    bwb3 = BWB3(bwb4)

    global eqns
    eqns = map(numerator, bwb3.lift().coefficients())

    #for poly in eqns:
    #    print poly

# PolynomialRing_expand() runs very slowly on helium, so I've tried to wrap my own version of it...

# probably doesn't work - more work has gone into numpy_expand
def bwb_expand():
  bwb4 = create_bwb4()
  for count,monomial in enumerate(bwb4.operands()):
    index = [monomial.degree(c) for c in cvars]
    map(operator.add, index, [1,1,1,0,0,0] * monomial.degree(SR.var('r1'))/2)
    map(operator.add, index, [0,0,0,1,1,1] * monomial.degree(SR.var('r2'))/2)
    map(operator.add, index, [1,1,1,-1,-1,-1] * monomial.degree(SR.var('r12'))/2)
    index = index + [monomial.degree(SR.var('r1')) % 2, monomial.degree(SR.var('r2')) % 2, monomial.degree(SR.var('r12')) % 2]
    print count, monomial, index, monomial.coefficient(map(operator.mul, v_s, index))

import numpy as np
term_expansion = dict()
equations = dict()

def numpy_expand():
  bwb4 = create_bwb4()
  for count,monomial in enumerate(bwb4.operands()):
    index = np.array([int(monomial.degree(c)) for c in cvars] + [0,0,0])
    index2 = np.array([int(monomial.degree(c)) for c in v_s])
    index3 = np.array([int(monomial.degree(c)) for c in SRr_s])
    try:
        terms = term_expansion[tuple(index3/2)]
    except KeyError:
        expansion = expand(mul(map(operator.pow, rvars, 2*(index3/2))))
        assert expansion.operator() is sage.symbolic.operators.add_vararg
        terms = [np.array([int(term.degree(c)) for c in cvars] + [0,0,0]) for term in expansion.operands()]
        term_expansion[tuple(index3/2)] = terms

    index = index + np.array([0,0,0,0,0,0, int(monomial.degree(SR.var('r1'))) % 2, int(monomial.degree(SR.var('r2'))) % 2, int(monomial.degree(SR.var('r12'))) % 2])
    #for term in terms:
    #    print count, monomial, index + term, monomial.coefficient(mul(map(operator.pow, v_s, index2)))
    multiplier = monomial.coefficient(mul(map(operator.pow, v_s, index2)))
    for term in terms:
        equations[tuple(index + term)] = equations.get(tuple(index + term), 0) + multiplier
    if count % 100 == 0: print count

#raise(None)

print
print

#Ss = solve (map(SR, eqns), *((E,) + Avars + Bvars), algorithm='sympy')
#for s in Ss: print s
#exit()

# Tried with a finite field.  Same performance issues
# (BWB's coefficient fields needs to be changed to ZZ)
#BWBff = PolynomialRing(GF(3), (var('E'),) + Avars + Bvars)
#BWB2ff.<x,y,z,r> = Frac(BWBff)[]
#hom = BWB.hom(BWBff)
#bwbIf = ideal(map(hom, bwbI.gens()))

# If we used PolynomialRing_expand(), then we can form an ideal from
# eqns and use Groebner basis techniques to factor the ideal.

def associated_primes():

    bwbI = ideal(eqns)
    primes = bwbI.associated_primes()

    for prime in primes:
        print prime._repr_short()


# Look for solutions using an approximate numerical technique

# Standard operator overloading lets us use Sage's multi-precision
# floating point numbers for most numpy operations, but a few need to
# be overridden to use Sage tests.  Use some Python magic to achieve
# this.

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

use_multiprocessing = True

if use_multiprocessing:

    def minfunc(v):
        d = dict(zip(coeff_vars, v))
        sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))
        res = real_type(sum_of_squares / zero_variety.subs(d))
        print res
        return res

    def jac(v):
        d = dict(zip(coeff_vars, v))
        sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))
        D_sum_of_squares = {var : sum(map(lambda x: x.get(), [cc.D_sum_of_squares(v,var) for cc in ccs])) for var in coeff_vars}
        zero_var = zero_variety.subs(d)
        dvar = {var : 0 for var in coeff_vars}
        dvar.update({var : d[var] for var in Avars})
        A = [real_type((D_sum_of_squares[var]*zero_var - 2*dvar[var]*sum_of_squares)/zero_var^2) for var in coeff_vars]
        res = np.array(A)
        return res

else:

    def minfunc(v):
        res = real_type(minpoly.subs(dict(zip(coeff_vars, v))) / zero_variety.subs(dict(zip(coeff_vars, v))))
        print res
        return res

    def jac(v):
        res = np.array([real_type(d.subs(dict(zip(coeff_vars, v)))) for d in minpoly_derivatives])
        return res

import random

def random_numerical(seed=0):

    import scipy.optimize

    # eqns.append(E+1/4)
    # eqns.append(BWB.gen(1) - 1)

    nvars = len(coeff_vars)

    # even though we're using numpy, we don't need to set its PRNG
    # seed, (which would require calling numpy.random.seed()), since
    # the algorithm is deterministic after the iv is picked
    random.seed(seed)        # for random
    set_random_seed(seed)    # for RR.random_element()

    global iv
    iv = [random.random() for i in range(nvars)]

    # We know the zero variety (all Avar's zero, so Psi is zero) will be
    # a "solution", but we want to avoid it

    global zero_variety, minpoly, minpoly_derivatives
    zero_variety = sum(map(square, Avars))
    if not use_multiprocessing:
        global minpoly, minpoly_derivatives
        minpoly = sum([poly*poly for poly in eqns])
        minpoly_derivatives = [diff(minpoly / zero_variety, v) for v in coeff_vars]

    global SciMin
    SciMin = scipy.optimize.minimize(minfunc, iv, jac=jac, method='BFGS', options={'return_all':True})

    #print SciMin

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

# print an input file for Bertini to numerically find all irreducible components

def bertini():
    print '''
CONFIG
  TrackType:1;
END;

INPUT

'''
    print 'variable_group ', ','.join(map(str, BWB.gens())), ';'

    print 'function ', ','.join(['f{}'.format(i) for i in range(len(eqns))]), ';'

    for i in range(len(eqns)):
        print 'f{} = '.format(i), eqns[i], ';'

    print 'END;'
