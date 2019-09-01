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
    global eq, H, Psi, A, B, Avars, Bvars, cvars, rvars

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
# serialization.
#
# BaseManager wraps its class instances with proxy objects that
# performs serialization and synchronization.  We use two classes -
# ManagerClass exists on the main process and WorkerClass exists
# on the spawned workers.
#
# It's possible to pass a proxy for the ManagerClass to each worker
# when the process is started, so there's only a single ManagerClass
# on the main process and proxies for it on each worker.
#
# It doesn't seem possible to pass proxy objects in method calls to
# other proxy objects, so the workers are identified by address, and a
# separate WorkerClass is created for each remote process that wants
# to make calls to that worker.  This allows the workers to pass data
# directly to other workers, avoiding turning the main process into a
# bottleneck.
#
# Since every connection to a worker process has a separate
# WorkerClass, I use global and class variables a lot in WorkerClass.
# The proxy to the main process's ManagerClass is stored in a global
# variable 'mc', as is the local server in a global variable
# 'local_server' (we need to know its address to identify ourselves to
# the main process).  The data being processes is stored in a class
# variable 'data' (an array).
#
# Right now, it works like this:
#
# load('helium.sage')
# prep_hydrogen()
# create_bwb4()
# SR_expand2a()
# ops = bwb4a.operands()
# start_manager_thread()
# mc.set_range(len(ops), 100)
# mc.start_worker()  (13 times)
# wc = mc.getq()[0]
# wc.get_data()

def multi_init():
    global ops
    prep_hydrogen()
    create_bwb4()
    SR_expand2a()
    ops = bwb4a.operands()
    start_manager_process()
    mc.set_range(len(ops), 100)
    #SRdict_background()



blocksize = 100

import queue

import multiprocessing, logging
logger = multiprocessing.get_logger()
logger.setLevel(logging.INFO)
if len(logger.handlers) == 0:
    multiprocessing.log_to_stderr()

running_workers = 2

# I want some of longer-running methods in WorkerClass to return
# immediately after arranging to start processing in background.  This
# decorator causes a method to run asynchronously in a background
# thread, as well as keeping a record of such threads.  WorkerClass
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

# Standard Python proxy objects from the multiprocessing library can't
# be used as keys in dictionaries because multiple proxies can point
# to the same object.  This code adjusts the hashing and equality
# methods to ensure that proxies to the same object test equal and can
# then be used as dictionary keys.
#
# Probably needs to be reported as a bug in the multiprocessing library.

def BaseProxy_hash(self):
    return hash((self._token.address, self._token.id))
def BaseProxy_eq(self, other):
    return hash(self) == hash(other)

from multiprocessing.managers import BaseProxy
BaseProxy.__hash__ = BaseProxy_hash
BaseProxy.__eq__ = BaseProxy_eq


class ManagerClass:
    localq = queue.Queue()
    workers = {}
    thousands = []
    worker_data_count = {}
#    def __init__(self):
#        self.localq = queue.Queue()
#        self.workers = {}
    def set_range(self, limit, blocksize):
        self.thousands.extend(range(0, limit, blocksize))
    # race conditions here
    def register_manager(self, mc):
        logger.debug('register_manager %s', mc)
        self.mc = mc
    def debug(self):
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        if server:
            logger.debug(server.id_to_obj)
    def autoself(self):
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        if server:
            for key,value in server.id_to_obj.items():
                if value[0] == self:
                    token = multiprocessing.managers.Token(typeid='ManagerClass', address=server.address, id=key)
                    proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
                    return proxy
        return None
    def start_worker(self):
        worker_manager = BaseManager()
        worker_manager.start()
        workerclass = worker_manager.WorkerClass()
        workerclass.register_manager(self.autoself())

        logger.debug('start_worker %s %s', workerclass, workerclass._token)
        #self.workers[workerclass] = worker_manager
        self.workers[workerclass] = (workerclass, worker_manager)
        #self.workers[0] = worker_manager

        global running_workers

        if self.localq.qsize() > 1:
            logger.debug('combine')
            self.worker_data_count[workerclass] = 2
            wc1 = self.localq.get()
            wc2 = self.localq.get()
            logger.debug('wc1 %s', wc1._token)
            logger.debug('wc2 %s', wc2._token)
            workerclass.get_data_from_worker(wc1)
            workerclass.get_data_from_worker(wc2)
            #del self.workers[wc1]
            #del self.workers[wc2]
            running_workers = running_workers - 2
        elif len(self.thousands) > 0:
            thousand = self.thousands.pop()
            logger.debug('%s expand (%d,%d)', self, thousand, thousand+blocksize)
            workerclass.start_expand(thousand, thousand+blocksize)
        else:
            running_workers = running_workers - 1
            self.shutdown_worker(workerclass)
            #workerclass.shutdown()
            #del self.workers[address]
    def shutdown_worker(self, wc):
        logger.debug('shutdown_worker %s %s', repr(wc), wc._token)
        #self.workers[wc][1].shutdown()
        wc.join_threads()
        del self.workers[wc]
    def worker_addresses(self):
        return self.workers.keys()
    def release_all_workers(self):
        for addr in self.workers.keys():
            logger.debug('droping %s', addr)
            del self.workers[addr]
    def notify_result_ready(self, worker):
        logger.debug('notify_result_ready %s', worker._token)
        global running_workers
        #running_workers = running_workers + 1
        self.localq.put(worker)
    def notify_data_received(self, worker):
        # thread safe because RPC messages are processed on a single thread
        logger.debug('notify_data_received %s', worker._token)
        self.worker_data_count[worker] = self.worker_data_count[worker] - 1
        if (self.worker_data_count[worker] == 0):
            # We do need to delete this dictionary item, otherwise the
            # reference associated with the key will prevent the
            # worker process from terminating.
            del self.worker_data_count[worker]
            worker.start_combine()
    def qsize(self):
        return self.localq.qsize()
    def getq(self):
        return self.localq.queue
    def thousands_len(self):
        #logger.debug('thousands_len')
        return len(self.thousands)
    def getpid(self):
        logger.debug('getpid')
        return os.getpid()

import time

class WorkerClass:
    data = []
    def getpid(self):
        return os.getpid()
    def join_threads(self):
        if hasattr(self, 'threads'):
            for th in self.threads:
                logger.debug('join %s', th)
                th.join()
    def autoself(self):
        server = getattr(multiprocessing.current_process(), '_manager_server', None)
        if server:
            for key,value in server.id_to_obj.items():
                if value[0] == self:
                    token = multiprocessing.managers.Token(typeid='WorkerClass', address=server.address, id=key)
                    proxy = multiprocessing.managers.AutoProxy(token, 'pickle', authkey=server.authkey)
                    logger.debug('autoself %s', token)
                    return proxy
        return None
    def register_manager(self, mc):
        self.mc = mc
    def get_mc(self):
        return self.mc
    @async_method
    def start_expand(self, start, stop):
        self.data.append(SRdict_expander2a(expand(sum(islice(ops, start, stop)))))
        self.mc.notify_result_ready(self.autoself())
    @async_method
    def send_result_to_worker(self, wc):
        wc.data_transfer(self.data[0])
    @async_method
    def get_data_from_worker(self, wc):
        logger.debug('get_data_from_worker %s', wc._token)
        self.data.extend(wc.get_data())
        self.mc.shutdown_worker(wc)
        self.mc.notify_data_received(self.autoself())
    def get_data(self):
        return self.data
    def data_transfer(self, datain):
        self.data.append(datain)
        self.mc.notify_data_received(self.autoself())
    @async_method
    def start_combine(self):
        logger.debug('start_combine')
        result = {}
        for SRd in self.data:
            for key,value in SRd.items():
                if isinstance(value, str):
                    value = eval(preparse(value))
                result[key] = result.get(key, 0) + value
        self.data = [result]
        self.mc.notify_result_ready(self.autoself())

# Not sure just how this works.  Basically cribbed from the Python
# multiprocessing docs.

from multiprocessing.managers import BaseManager
BaseManager.register('ManagerClass', ManagerClass)
BaseManager.register('WorkerClass', WorkerClass)

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


import random

def random_numerical(seed=0):

    import scipy.optimize

    # eqns.append(E+1/4)
    # eqns.append(BWB.gen(1) - 1)
    minpoly = sum([poly*poly for poly in eqns])

    nvars = len(coeff_vars)

    # even though we're using numpy, we don't need to set its PRNG
    # seed, (which would require calling numpy.random.seed()), since
    # the algorithm is deterministic after the iv is picked
    random.seed(seed)        # for random
    set_random_seed(seed)    # for RR.random_element()

    global iv
    iv = [random.random() for i in range(nvars)]

    def square(x):
        return x*x

    # We know the zero variety (all Avar's zero, so Psi is zero) will be
    # a "solution", but we want to avoid it

    zero_variety = sum(map(square, Avars))

    # The function we're trying to minimize: the sum of squares of the polys
    # that define the solution variety, divided by two factors we're trying
    # to avoid: the norm of 'v' (avoid the origin), and zero_variety (which
    # includes the origin).

    real_type = np.float64
    #real_type = RR
    #real_type = RealField(100)

    global minfunc
    def minfunc(v):
        #return minpoly.subs(dict(zip(BWB.gens(), v))) / sum(map(square, v)) / zero_variety.subs(dict(zip(BWB.gens(), v)))
        return real_type(minpoly.subs(dict(zip(coeff_vars, v))) / zero_variety.subs(dict(zip(coeff_vars, v))))
        #return minpoly.subs(dict(zip(BWB.gens(), v)))

    minpoly_derivatives = [diff(minpoly / zero_variety, v) for v in coeff_vars]
    global jac
    def jac(v):
        return np.array([real_type(d.subs(dict(zip(coeff_vars, v)))) for d in minpoly_derivatives])

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
