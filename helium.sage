# -*- mode: python -*-
#
# Python code to search for solutions of Hydrogen and Helium
# non-relativistic time-invariant Schrodinger equation
#
# INTERACTIVE USAGE:
#
# load('helium.sage')
# prep_hydrogen()
# multi_init()
# multi_expand()
# random_numerical()
#
# This will produce a solution to the hydrogen atom using its
# default ansatz (number 1).  You can also prep_helium(), and
# supply an ansatz number as argument, as in prep_helium(-7).
#
# Negative ansatzen use a spherically symmetric Hamiltonian that
# reduces the dimension of the problem and lets us use the faster
# FLINT implementation because there are no roots in the Hamiltonian
# and FLINT doesn't have a Groebner basis implementation, which is
# required to handle roots.
#
# random_numerical() also takes an optional argument, the initial
# random seed, along with an optional 'homogenize' argument.
#
# "Homogenization" (not a very good term) is used to force selected
# polynomials to be non-zero by forcing their coefficients to be 1,
# one at a time.  Different values of the 'homogenize' argument
# (starting at 0) force different coefficients to be 1, and an
# exception is thrown once all of the homogenization possibilities
# have been exhausted.  Specifying homogenize=-1 runs through
# all possible homogenizations (usually what you want).
#
# Homogenization can also be done when the trial solution is
# constructed (this is the 'homogenize' argument to the
# trial_polynomial function), which produces systems with fewer free
# coefficients, but a new system has to be constructed for every
# homogenization possibility, so I don't do it this way, and only plan
# to develop this option if computational complexity becomes an issue.
#
# CONCEPT:
#
# We have a differential equation that we're trying to solve and a
# trial solution with lots of free parameters that we try to adjust
# to get a solution.
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
# latest version - March 2023
#
# no rights reserved; you may freely copy, modify, or distribute this
# program
#
# TODO list:
# - collector should parse directly into matrix form
# - more easily turn multiprocessing on and off
# - allow worker processes on different hosts
# - check collected coefficient polynomials to see if they factor
# - automate finding polynomial relations
# - save checkpoints of optimization iterations
# - optimize scipy sparse matrix by vector multiplication
# - optimize creation of multi-vectors
# - allow workers to be added or removed on the fly

# If True, use TCP/IP connections to interconnect Python sub-processes,
# otherwise use UNIX sockets.  TCP/IP has much slower connection setups,
# but allows multiple hosts to be used.

use_tcpip_multiprocessing = False

import platform
import glob
import psutil
import random

current_process = psutil.Process(os.getpid())

import itertools
from itertools import *
import scipy.optimize

from sage.symbolic.operators import add_vararg, mul_vararg

from sage.rings.polynomial.polydict import ETuple

# from python docs
def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)

# from python docs
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# Python itertools docs contains equivalent code to combinations_with_replacement
#
# i'm modifying it to pair each item in iterable with a maximum count of times it can be used

def combinations_with_replacement_with_maxcount_BROKEN(iterable, r):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        i = 0
        while i < r:
            # If we're picking three from ABCDE and D and E are limited to one, then
            # we want ADC -> ADD (rejected) -> ADE -> AEE (rejected) -> BBB
            for i in reversed(range(r)):
                if indices[i] != n - 1:
                    break
            else:
                return

            #if indices[i] == n - 1:
                # break back to outer while loop
                # back i up even further than selected above
                # try again
                # how do we do all this?
            while (i < r) and (indices[i] < n - 1):
                this_items_max = maxcount[indices[i] + 1]
                if this_items_max >= (r - i):
                    indices[i:] = [indices[i] + 1] * (r - i)
                    i = r
                    yield tuple(pool[i] for i in indices)
                    break
                elif space_is_left:
                    indices[i:i+this_items_max] = [indices[i] + 1] * this_items_max
                    i += this_items_max
                else:
                    # indices[i] == n - 1 (or is about to be)
                    # no space is left
                    indices[i:] = [indices[i] + 1] * (r - i)
                    i = r
                    # reject

        yield tuple(pool[i] for i in indices)

def trial_polynomial(base, coordinates, roots, degree, homogenize=None, constant=True):
    """trial_polynomial(base, coordinates, roots, degree, homogenize=None, constant=True)
    Form a trial polynomial in the Symbolic Ring

    base is a string to which we append numbers to get our coefficient names; i.e, 'a' -> (a0,a1,a2,...)
    coordinates is a tuple of symbolic expressions (currently all symbols; i.e, x1.is_symbol() == True)
    roots is a tuple of symbolic expressions for our roots (currently all powers; i.e, r.operator() == pow)
    degree is maximum degree of the trial polynomial
    constant=None is optional and drops the constant term (essentially mindeg=1 instead of mindeg=0)
    homogenize=N is unused right now and is intended as a future performance optimization
    """

    # base is a string to which we append numbers to get our coefficient names; i.e, 'a' -> (a0,a1,a2,...)
    # cterms are coefficient terms
    # rterms are radius terms
    cterms = flatten([combinations_with_replacement(coordinates, d) for d in range(degree+1)])
    # use this 'terms' for real
    # the 'roots' are assumed to be square roots, so all we need is their powerset;
    #    any higher powers will be replaced with expansions
    terms = list(map(mul, (product(map(mul, cterms), map(mul, powerset(roots))))))
    # The first term is the constant 1, so if we don't want a constant term, drop it
    if not constant:
        terms = terms[1:]
    # use this 'terms' for testing
    # terms = list(map(mul,cterms)) + list(roots)
    coefficients = [var(base+str(c)) for c in range(len(terms))]
    poly_coefficients = list(coefficients)
    if homogenize != None:
        # homogenize: use 1 as the coefficient of the homogenize'th term
        #   and set all previous terms to 0.
        # The idea is to prevent a polynomial from being zero by running successive
        #   calculations running through all the terms of the polynomial, forcing them to be 1.
        homogenize_coefficient = homogenize % len(coefficients)
        if homogenize_coefficient > 0:
            poly_coefficients[0:homogenize_coefficient] = [0] * (homogenize_coefficient)
        poly_coefficients[homogenize_coefficient]= 1
        del coefficients[0:homogenize_coefficient+1]
    poly = sum([poly_coefficients[c]*v for c,v in enumerate(terms)])
    return (tuple(coefficients), poly)

# Energy constant in Schroedinger's equations
var('E')

def Del(Psi,vars):
    return sum([diff(Psi,v,2) for v in vars])

# Create an operator (DD) that takes derivatives of symbolic functions
# w.r.t. their arguments, i.e, `DD[0](Phi)(A)` is the derivative of
# Phi w.r.t. its 0-th argument, `A`.
#
# from https://groups.google.com/g/sage-devel/c/xBHw11qUARg/m/0eqj3eUFsFkJ
# referenced from https://trac.sagemath.org/ticket/17445

from sage.symbolic.operators import FDerivativeOperator
class Doperator:
  def __init__(self,vars=None):
    self.vars= [] if vars is None else vars

  def __call__(self,f):
    return FDerivativeOperator(f,self.vars)

  def __getitem__(self,i):
    if isinstance(i,tuple):
       newvars=self.vars+list(i)
    else:
       newvars=self.vars+[i]
    return Doperator(newvars)

DD=Doperator()

def finish_prep(ansatz):
    global eq, H, coeff_vars, ODE_vars, coordinates, roots
    global A,B,C,D,F,G,M,N
    global homogenize_groups

    (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
    (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
    (Cvars, C) = trial_polynomial('c', coordinates, roots, 1)
    (Dvars, D) = trial_polynomial('d', coordinates, roots, 1)
    (Fvars, F) = trial_polynomial('f', coordinates, roots, 1)
    (Gvars, G) = trial_polynomial('g', coordinates, roots, 1)

    coeff_vars = (E,) + Avars + Bvars + Cvars + Dvars + Fvars + Gvars

    SR_function = sage.symbolic.function_factory.function

    if ansatz == 1:
        # A linear polynomial times the exponential of a linear polynomial
        # Phi is an exponential of a linear polynomial
        # Phi = e^B, so diff(Phi,B) = Phi and diff(Phi,v) = diff(B,v)*Phi
        # A is a linear polynomial; the solution is A times Phi.
        # Homogenization forces A to be non-zero.
        #
        # An earlier verion of this ansatz was used extensively for testing Hydrogen
        Phi = SR_function('Phi')
        (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
        Psi = A * Phi(B)
        pre_subs = {DD[0](Phi)(B) : Phi(B), DD[0,0](Phi)(B) : Phi(B)}
        post_subs = {Phi(B) : SR.var('Phi')}
        homogenize_groups = (Avars,)
        coeff_vars = (E,) + Avars + Bvars
        ODE_vars = ('Phi', )

    elif ansatz == 2:
        # A linear polynomial times the logarithm of a linear polynomial
        # Xi is a logarithm; Xi = ln C, so diff(Xi,C) = 1/C and diff(Xi,v) = diff(C,v)/C
        # A is a linear polynomial; the solution is A times Xi.
        # Homogenization forces A and C to be non-zero
        #
        # I've used this ansatz very little.
        Xi = SR_function('Xi')
        Psi = A * Xi(C)
        pre_subs = {DD[0](Xi)(C) : 1/C, DD[0,0](Xi)(C) : -1/C^2}
        post_subs = {Xi(C) : SR.var('Xi')}
        homogenize_groups = (Avars, Cvars)
        coeff_vars = (E,) + Avars + Cvars
        ODE_vars = ('Xi', )

    elif ansatz == 3:
        # Chi is a weird second-order mess: C d^2 Chi/dB^2 - D dChi/dB - F Chi - G = 0
        # A is a linear polynomial; the solution is A times Chi.
        #
        # I declared Chi to be a function of B, but in retrospect, it's also a
        # function of C, D, F, and G, none of which are (necessarily) functions of B
        #
        # Homogenization forces A, B and C to be non-zero, and B is non-constant.
        #
        # I've used this ansatz very little; ansatz 4 is closely related
        Chi = SR_function('Chi')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1, constant=None)
        Psi = A*Chi(B)

        homogenize_groups = (Avars, Bvars, Cvars)
        coeff_vars = (E,) + Avars + Bvars + Cvars + Dvars + Fvars + Gvars

        pre_subs = {DD[0,0](Chi)(B) : (D/C * DD[0](Chi)(B) + F/C * Chi(B)) + G/C}
        post_subs = {Chi(B) : SR.var('Chi'), DD[0](Chi)(B) : SR.var('DChi')}
        ODE_vars = ('Chi', 'DChi')

    elif ansatz == 4:
        # Like ansatz 3, but without the polynomial A as a factor, and thus simplier
        #
        # Homogenization forces B and C to be non-zero, and B is non-constant.
        #
        # This ansatz is fairly well explored, but in an earlier version of the code
        # (pre-edbaa8 and pre-e74ded) whose ring and class structures aren't compatible.
        Chi = SR_function('Chi')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1, constant=None)
        Psi = Chi(B)

        homogenize_groups = (Bvars, Cvars)
        coeff_vars = (E,) + Bvars + Cvars + Dvars + Fvars + Gvars

        pre_subs = {DD[0,0](Chi)(B) : (D/C * DD[0](Chi)(B) + F/C * Chi(B)) + G/C}
        post_subs = {Chi(B) : SR.var('Chi'), DD[0](Chi)(B) : SR.var('DChi')}
        ODE_vars = ('Chi', 'DChi')

    elif ansatz == 5:
        # A second-order homogeneous ODE: D(B) d^2 Zeta/dB^2 - M(B) dZeta/dB - N(B) Zeta = 0
        # where D(B), M(B), and N(B) are linear polynomials in B, which is itself a linear polynomial
        #
        # Homogenization forces B and D to be non-zero; B is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1, constant=None)
        Psi = Zeta(B)
        (Dvars, D) = trial_polynomial('d', [B], [], 1)
        (Mvars, M) = trial_polynomial('m', [B], [], 1)
        (Nvars, N) = trial_polynomial('n', [B], [], 1)

        homogenize_groups = (Dvars, Bvars)

        coeff_vars = (E,) + Bvars + Dvars + Mvars + Nvars

        pre_subs = {DD[0,0](Zeta)(B) : (M * DD[0](Zeta)(B) + N * Zeta(B)) / D}
        post_subs = {Zeta(B) : SR.var('Zeta'), DD[0](Zeta)(B) : SR.var('DZeta')}
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 6:
        # A second-order homogeneous ODE: D(B/C) d^2 Zeta/d(B/C)^2 - M(B/C) dZeta/d(B/C) - N(B/C) Zeta = 0
        # where D(B/C), M(B/C), and N(B/C) are linear polynomials in B/C, a first-degree rational function
        Zeta = SR_function('Zeta')
        Psi = Zeta(B/C)
        (Dvars, D) = trial_polynomial('d', [B/C], [], 1)
        (Mvars, M) = trial_polynomial('m', [B/C], [], 1)
        (Nvars, N) = trial_polynomial('n', [B/C], [], 1)

        coeff_vars = (E,) + Bvars + Cvars + Dvars + Mvars + Nvars

        pre_subs = {DD[0,0](Zeta)(B/C) : (M * DD[0](Zeta)(B/C) + N * Zeta(B/C)) / D}
        post_subs = {Zeta(B/C) : SR.var('Zeta'), DD[0](Zeta)(B/C) : SR.var('DZeta')}
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 7:
        # A second-order homogeneous ODE: D(B/C) d^2 Zeta/d(B/C)^2 - M(B/C) dZeta/d(B/C) - N(B/C) Zeta = 0
        # where D(B/C), M(B/C), and N(B/C) are second-degree polynomials in B/C, a second-degree rational function
        Zeta = SR_function('Zeta')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 2, homogenize=0)
        (Cvars, C) = trial_polynomial('c', coordinates, roots, 2, homogenize=1)
        Psi = Zeta(B/C)
        (Dvars, D) = trial_polynomial('d', [B/C], [], 2, homogenize=0)
        (Mvars, M) = trial_polynomial('m', [B/C], [], 2)
        (Nvars, N) = trial_polynomial('n', [B/C], [], 2)

        coeff_vars = (E,) + Bvars + Cvars + Dvars + Mvars + Nvars

        pre_subs = {DD[0,0](Zeta)(B/C) : (M * DD[0](Zeta)(B/C) + N * Zeta(B/C)) / D}
        post_subs = {Zeta(B/C) : SR.var('Zeta'), DD[0](Zeta)(B/C) : SR.var('DZeta')}
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 8:
        # A first-order homogeneous ODE: M(B) dZeta/dB - N(B) Zeta = 0
        # where M(B) and N(B) are linear polynomials in B, which is itself a linear polynomial
        # B can not be constant; neither B or M can be zero (homogenization)
        #
        # Logically it's a step backwards from ansatz 5, but I want to see it work.

        Zeta = SR_function('Zeta')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1, constant=None)
        Psi = Zeta(B)
        (Mvars, M) = trial_polynomial('m', [B], [], 1)
        (Nvars, N) = trial_polynomial('n', [B], [], 1)

        homogenize_groups = (Mvars, Bvars)
        coeff_vars = (E,) + Bvars + Mvars + Nvars

        # A limitation of the program is that I have to manually calculate DD[0,0](Zeta)(B) here
        #  DD[0,0](Zeta)(B)  = d^2 Zeta / dB^2 = d/dB (N(B) * Zeta(B) / M(B))
        #     = (dN/dB * Zeta(B) * M(B) + N(B) * DD[0](Zeta)(B) * M(B) - N(B) * Zeta(B) * dM/dB ) / M^2(B)
        #     = (n1 * Zeta(B) * M(B) + N(B) * DD[0](Zeta)(B) * M(B) - N(B) * Zeta(B) * m1 ) / M^2(B)
        #
        # Can't write diff(M,B) because B is a polynomial and diff only accepts a symbol as its second argument.
        # Yet we know that M = m1*B + m0, so diff(M,B)=m1
        m1 = Mvars[1]
        n1 = Nvars[1]
        pre_subs = {DD[0](Zeta)(B) : (N * Zeta(B)) / M,
                    DD[0,0](Zeta)(B) : (n1 * Zeta(B) * M + N * N * Zeta(B) - N * Zeta(B) * m1 ) / (M*M) }
        post_subs = {Zeta(B) : SR.var('Zeta')}
        ODE_vars = ('Zeta', )

    elif ansatz == 9:
        # A first-order homogeneous ODE: dZeta/dB - n0 Zeta = 0
        # where n0 is a constant and B is a linear polynomial
        #
        # Logically it's a further step backwards from ansatz 8, but ansatz 8 has too many free variables
        # for scipy.optimize.root to work on 1-dim hydrogen if we use homogenization
        Zeta = SR_function('Zeta')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1, constant=None)
        Psi = Zeta(B)
        (Nvars, N) = trial_polynomial('n', [B], [], 0)

        homogenize_groups = (Bvars, )
        coeff_vars = (E,) + Bvars + Nvars

        # A limitation of the program is that I have to manually calculate DD[0,0](Zeta)(B) here
        #  DD[0](Zeta)(B)  = n0 Zeta(B)
        #  DD[0,0](Zeta)(B)  = n0^2 Zeta(B)

        n0 = Nvars[0]
        pre_subs = {DD[0](Zeta)(B) : n0 * Zeta(B),
                    DD[0,0](Zeta)(B) : n0^2 * Zeta(B)}
        post_subs = {Zeta(B) : SR.var('Zeta')}
        ODE_vars = ('Zeta', )

    elif ansatz == 10:
        # A second-order homogeneous ODE: D(B) d^2 Zeta/dB^2 - M(B) dZeta/dB - N(B) Zeta = 0
        # where D(B), M(B), and N(B) are quadratic polynomials in B, which is itself a quadratic polynomial
        #
        # Homogenization forces B and D to be non-zero; B is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 2, constant=None)
        Psi = Zeta(B)
        (Dvars, D) = trial_polynomial('d', [B], [], 2)
        (Mvars, M) = trial_polynomial('m', [B], [], 2)
        (Nvars, N) = trial_polynomial('n', [B], [], 2)

        homogenize_groups = (Dvars, Bvars)

        coeff_vars = (E,) + Bvars + Dvars + Mvars + Nvars

        pre_subs = {DD[0,0](Zeta)(B) : (M * DD[0](Zeta)(B) + N * Zeta(B)) / D}
        post_subs = {Zeta(B) : SR.var('Zeta'), DD[0](Zeta)(B) : SR.var('DZeta')}
        ODE_vars = ('Zeta', 'DZeta')

    else:
        raise 'Bad ansatz'

    eq = H(Psi) - E*Psi

    eq = eq.subs(pre_subs).subs(post_subs)

    # reduce coeff_vars to those which actually appear in the equation
    coeff_vars = tuple(sorted(set(eq.free_variables()).intersection(coeff_vars), key=lambda x:str(x)))
    coeff_vec = np.array(coeff_vars)

def prep_hydrogen(ansatz=1):
    global H, coordinates, roots

    if ansatz < 0:

        var('r')
        coordinates = (r,)
        roots = tuple()

        def H(Psi):
            return - 1/2 * (1/r^2 * diff(r^2 * diff(Psi,r), r)) - (1/r)*Psi

    else:
        var('x1,y1,z1')
        coordinates = (x1,y1,z1)

        global r1
        r1 = sqrt(x1^2+y1^2+z1^2)
        roots = (r1,)

        def H(Psi):
            return - 1/2 * Del(Psi,[x1,y1,z1]) - (1/r1)*Psi

    finish_prep(ansatz=abs(ansatz))

def prep_helium(ansatz=6):
    global H, coordinates, roots

    if ansatz < 0:
        # According to Nakatsuji, we can write the helium Hamiltonian
        # for S states (no angular momentum) in a r1/r2/r12 coordinate system.
        #
        # That's what we do for a negative ansatz

        var('R1,R2,R12')
        coordinates = (R1,R2,R12)
        roots = tuple()

        def H(Psi):
            return - 1/2 *sum(diff(Psi, Ri, 2) + 2/Ri*diff(Psi,Ri) for Ri in [R1,R2])  \
                   - (diff(Psi, R12, 2) + 2/R12*diff(Psi,R12))                          \
                   - (R1^2 + R12^2 - R2^2)/(2*R1*R12) * diff(diff(Psi,R12),R1)         \
                   - (R2^2 + R12^2 - R1^2)/(2*R2*R12) * diff(diff(Psi,R12),R2)         \
                   - sum(2/Ri for Ri in [R1,R2])*Psi + 1/R12*Psi

    else:

        var('x1,y1,z1,x2,y2,z2')

        global r1, r2, r12
        r1 = sqrt(x1^2+y1^2+z1^2)
        r2 = sqrt(x2^2+y2^2+z2^2)
        r12 = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)

        coordinates = (x1,y1,z1, x2,y2,z2)
        roots = (r1,r2,r12)

        def H(Psi):
            return - 1/2 * Del(Psi,[x1,y1,z1]) - 1/2 * Del(Psi,[x2,y2,z2]) - (2/r1)*Psi - (2/r2)*Psi + (1/r12)*Psi

    finish_prep(ansatz=abs(ansatz))


# Now we want to replace all of the sqrt(...) factors with 'r',
# and we use a clever Python trick to build a dictionary
# that maps expressions to variable names.

def varName(var):
    for name,value in globals().items():
        if id(var) == id(value):
            return name
    return None

def mk_maps(roots):
    return {v.operands()[0] : SR.var(varName(v)) for v in roots}

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
    maps = mk_maps(roots)
    # next, convert all of the roots in the equation to use the r-variables
    eq_a = roots_to_rs(eq)

# Print "equation A" in the form parsed by GNU EMACS's outline mode.
# This function is only used for debugging.

def analyze_eq_a(eq, depth=1, print_depth=2, file=sys.stdout):
    if eq.operator():
       print('*' * depth, eq.operator(), len(eq.operands()), file=file)
       if depth >= print_depth: print(eq, file=file)
       for o in eq.operands():
           analyze_eq_a(o, depth+1, print_depth=print_depth, file=file)


# Create a polynomial ring to hold our expressions.
#
# Sage does this using Singular, which stores polynomials internally in standard
# form (i.e, fully expanded).  Singular is the default, while FLINT is available
# as an option, which I tend to use because I've worked on the FLINT code and
# found it more accessible than Singular.
#
# The variables are listed from most significant to least significant.  We use lexicographic
# ordering to group terms together conveniently.
#
# I want 'roots' to be first in the ordering, because they're going to be substituted for,
# so making them the most significant variables groups like 'roots' terms together.
#
# I want the coeff_vars to be last in the ordering, because the system we're going to solve
# will be grouped by coordinates and ODE_vars.
#
# 'encoding' is an option that I've added to the Sage/FLINT implementation to describe
# how the polynomial terms will be written out to disk.  dexlex64(N) encodes N variables
# using 64-bit deglex encoding (see M. Gastineau, Storage of Multivariate Polynomials,
# Advanced School on Specific Algebraic Manipulators, 2007); sint64 encodes the coefficient
# as a signed 64-bit integer, so this encoding encodes each polynomial term using
# three 64-bit words.  If things don't fit, it throws an exception.

def create_polynomial_ring():
    global R,F,num_rvars,num_cvars
    roots_names = list(map(varName, roots))
    num_rvars = len(roots_names) + len(ODE_vars) + len(coordinates)
    num_cvars = len(coeff_vars)
    encoding = 'deglex64({}),deglex64({}),sint64'.format(num_rvars, num_cvars)
    if len(roots) > 0:
        # FLINT multivariates can't handle reduction modulo an ideal, so use Singular multivariates instead
        print('Using Singular implementation')
        R = PolynomialRing(ZZ, names=tuple(flatten((roots_names, ODE_vars, coordinates, coeff_vars))), order='lex')
    else:
        print('Using FLINT implementation')
        R = PolynomialRing(ZZ, names=tuple(flatten((roots_names, ODE_vars, coordinates, coeff_vars))),
                           implementation="FLINT", order='lex', encoding=encoding)
    F = Frac(R)

def mk_ideal(R, roots):
    "Given a list or tuple of roots, return a ideal of ring R that reduces the global variable names of those roots"
    ideal_generators = []
    for v in roots:
        assert v.operator() is operator.pow
        Rname = R(varName(v))
        Rexpr = R(v.operands()[0])
        power = int(1/v.operands()[1])
        ideal_generators.append(Rname^power - Rexpr)
    return ideal(ideal_generators)

def convert_eq_a():
    global F_eq_a, F_eq_a_n, F_eq_a_d
    create_polynomial_ring()
    # If we write this as 'F_eq_a = F(eq_a)', Sage will attempt to construct F_eq_a by calling
    # eq_a.numerator() and eq_a.denominator(), which will perform lots of rational function
    # math in the Symbolic Ring, which is very slow and memory intensive.  Calling it
    # like 'eq_a.polynomial(ring=F)' recurses through the expression tree and builds the
    # expression from the bottom up using polynomial ring operations, which are much more efficient.
    #
    # This trick (currently) only works on my development Sage, so try it and fall back on the slower way.
    try:
        F_eq_a = eq_a.polynomial(ring=F)
    except:
        F_eq_a = F(eq_a)

    # clear higher powers of roots
    if len(roots) > 0:
        F_eq_a_n = F_eq_a.numerator().mod(mk_ideal(R, roots))
        F_eq_a_d = F_eq_a.denominator().mod(mk_ideal(R, roots))
    else:
        F_eq_a_n = F_eq_a.numerator()
        F_eq_a_d = F_eq_a.denominator()
    print('F_eq_a: numerator', F_eq_a_n.number_of_terms(), 'terms; denominator', F_eq_a_d.number_of_terms(), 'terms')

import time

def timefunc(func, *args):
    start_time = time.perf_counter()
    func(*args)
    end_time = time.perf_counter()
    print('{:30} {:10.2f} sec'.format(func.__name__, end_time - start_time))

def multi_init():
    timefunc(create_eq_a)
    timefunc(convert_eq_a)

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
# until we've calculated the polynomial that we want to expand.
# This avoids having to serialize that polynomial, since each
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
blocksize = 100000

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
    mc.set_range(F_eq_a_n.number_of_terms(), blocksize)

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
#   - handle CNTL-C in interactive session without killing subprocesses
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
        self.lists = [[]] * len(collectors)
    @async_method
    def shutdown(self):
        for i in range(len(self.collectors)):
            del self.collectors[-1]
    @async_method
    def start_expand(self, start, stop):
        # sleep for a fraction of a second here to let this method
        # return back to the manager before we begin work
        time.sleep(float(0.1))
        # Each term is presented as a tuple of ETuple and cofficient
        # iterator_exp_coeff was introduced later than Sage 9.0, so it doesn't work on Ubuntu 20's stock Sage
        try:
            iterator = F_eq_a_n.iterator_exp_coeff()
        except:
            iterator = F_eq_a_n.dict().items()
        for etuple, coeff in islice(iterator, start, stop):
            # the first rvars form the key
            key = etuple[:num_rvars]
            # the remaining variables form the value
            value = etuple[num_rvars:]
            # hash the key to figure which collector it's assigned to
            listnum = hash(key) % len(self.collectors)
            # save in a list to be sent to the collectors
            self.lists[listnum].append((key, coeff, value))
        # expansion finished; send the expanded monomials to the collector processes
        for i in range(len(self.collectors)):
            self.collectors[i].combine_data(self.lists[i])
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
    def combine_data(self, lst):
        for key,coeff,value in lst:
            # key is an ETuple of rvars; coeff is an integer (not sure what type); value is an ETuple of the cvars
            if key not in self.result:
                self.result[key] = []
            self.result[key].append((coeff, value))

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
        return sum([len(v) for v in self.result.values()])

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
    # [a,b,c]
    # [a,b,c] * a, [b,c] * b, [c] * c = [a^2, a*b, a*c, b^2, b*c, c^2]
    #
    # Graded lexicographic order
    #
    # generate_multi_vector([a,b,c])
    #   -> [a,b,c,a^2,a*b,a*c,b^2,b*c,c^2,a^3,a^2*b,a^2*c,a*b^2,a*b*c,a*c^2,a*b^2,b^2*c,b*c^2,c^3]

    # XXX don't really need this in a method

    def generate_multi_vector(self, v):
        start_time = time.time()
        npv = np.array(v)
        stack = [[1], npv]
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
        D_stack = [[0], firsts]
        stack = [[1], npv]

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
        # so start at -1 and reshape the matrix every time we hit a monomial
        # that's too big.
        self.max_degree = -1
        self.dok = scipy.sparse.dok_matrix((len(self.result), 0), np.int64)

        # we don't care about the keys (roots, ODE vars, coordinates) because they were just used to accumulate the values
        # the values include both the constant coefficient and the value (the coefficient variable names)
        for l in self.result.values():
            #print(len(l))
            for coeff, value in l:
                if value.unweighted_degree() > self.max_degree:
                    # increase max_degree and rebuild indices
                    self.max_degree = value.unweighted_degree()
                    # index of the smallest tuple of the next higher degree
                    veclen=encode_deglex([self.max_degree + 1] + [0]*(len(value) - 1))
                    self.dok.resize((len(self.result), veclen))
                index = encode_deglex(value)
                #if decode_deglex(index + 1, len(value)) != list(value):
                #    print(decode_deglex(index + 1, len(value)), list(value))
                self.dok[self.i, index] += coeff
            self.i += 1
            #if (self.i % 1000 == 0):
            print(self.i, "of", len(self.result), "done")

        #self.M = np.unique(self.M, axis=0)
        #self.M = self.dok.tocsr()
        self.M = sp_unique(self.dok, axis=0, new_format='csr')
        logger.info('convert_to_matrix done')

    def load_from_polynomial(self, poly, section=int(0), total_sections=int(1)):
        r"""
        Loads a polynomial and stores it in a sparse matrix.

        The optional "section" arguments allow only a section of the
        row space to be loaded, to reduce the memory footprint.  The
        resulting matrices will either have to be stacked together
        to form a single matrix, or used with code designed to
        distribute the matrix operations over several machines.
        """

        assert type(section) == int
        assert type(total_sections) == int

        polyR = poly.parent()
        polyR_coeff_vars = tuple(map(polyR, coeff_vars))

        coeff_indices = tuple(i for i in range(len(coeff_vars)) if polyR.gens()[i] in polyR_coeff_vars)

        # To speed the code, require the coefficient variables to be the first set of generators
        len_coeff_indices = len(coeff_indices)
        assert coeff_indices == tuple(range(0, len_coeff_indices))

        for mon_tuple, coeff in poly._iter_ETuples():

            # Much slower - use split method if available
            #coeff_tuple = ETuple([(mon_tuple[i] if i in coeff_indices else 0) for i in all_indices])
            #row_tuple = ETuple([(mon_tuple[i] if i not in coeff_indices else 0) for i in all_indices])

            (coeff_tuple, row_tuple) = mon_tuple.split(len_coeff_indices - 1)

            if hash(row_tuple) % total_sections == section:

                if row_tuple not in self.rows:
                    self.rows[row_tuple] = len(self.rows)
                    self.dok.resize((len(self.rows), len(self.indices)))
                row = self.rows[row_tuple]

                while coeff_tuple not in self.indices:
                    self.max_degree += 1
                    self.indices = {next(pair[1]._iter_ETuples())[0] : pair[0] for pair in enumerate(self.generate_multi_vector(polyR_coeff_vars))}
                    # self.indices = {list(pair[1].dict().keys())[0] : pair[0] for pair in enumerate(self.generate_multi_vector(coeff_vars))}
                    self.dok.resize((len(self.rows), len(self.indices)))
                index = self.indices[coeff_tuple]

                self.dok[row, index] += coeff

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


def fns(v):
    r"""
    Evaluate our polynomials at a given coordinate.

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - a numpy vector of all of our polynomials, evaluated at ``vec``

    ALGORITHM:

    Call the `eval_fns` method for all CollectionClass's
    in parallel, then concatenate all of the results together,
    adding additional polynomials to enfore "homogenize" conditions.
    """

    # Save a copy of vector to aid in stopping and restarting the
    # calculation.  An explicit call to the copy method is required if
    # we're using the Fortran minpack code (i.e, scipy's optimize
    # package) because in that case, 'v' is only a pointer to a
    # Fortran array that gets deallocated once the Fortran code exits.
    # (I think - the copy's definitely needed, though)
    global last_v
    last_v = v.copy()

    homogenize_terms = [v[i] for i in homogenize_zero_indices] + [v[i] - 1 for i in homogenize_one_indices]
    res = np.hstack(list(map(lambda x: x.get(), [cc.eval_fns(v) for cc in ccs])) + homogenize_terms)

    global last_time
    sum_of_squares = sum(square(res))
    if last_time == 0:
        print(sum_of_squares)
    else:
        print("{:<30} {:10.2f} sec".format(sum_of_squares, time.time()-last_time))
    last_time = time.time()
    return res

def jac_fns(v):
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

    dN = np.vstack(list(map(lambda x: x.get(), [cc.jac_fns(v) for cc in ccs])) + [homogenize_derivatives])
    return dN

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

def jac_fns_shm(v):
    r"""
    Same as jac_fns method, but intended to be faster, uses shared
    memory to distribute the computation over multiple CollectorClass's
    running in different processes.
    """

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
    dN = np.vstack(list(map(get_nparray(shms), [cc.jac_fns(v, mdv_shm.name) for cc in ccs])) + [homogenize_derivatives])

    jac_time = time.time()
    print("   Compute Jacobian matrix {:7.2f} sec".format(jac_time - mdv_time))

    mdv_shm.close()
    mdv_shm.unlink()

    for shm in shms:
        shm.close()
        shm.unlink()

    return dN

def sum_of_squares(v):
    return sum(map(lambda x: x.get(), [cc.sum_of_squares(v) for cc in ccs]))

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

def random_numerical(seed=0, homogenize=None, limit=None):

    iv = make_iv(seed)

    global SciMin

    global homogenize_zeros, homogenize_ones
    global homogenize_zero_indices, homogenize_one_indices
    global homogenize_derivatives

    if homogenize != None:

        if homogenize == -1:
            # This is "auto-homogenization", where we run through all the homogenization possibilities
            # until we get the ValueError below
            hom = 0
            while True:
                random_numerical(homogenize=hom, seed=seed)
                hom = hom + 1

        # We perform "homogenization" (more like de-homogenization) by adding some simple
        # equations to our set of equations to be solved.  This is done in fns();
        # here we only need to identify which variables are to be set to zero
        # and which variables are to be set to one.

        homogenize_zeros = []
        homogenize_ones = []

        for grp in homogenize_groups:
            homogenize_coefficient = homogenize % len(grp)
            homogenize = int(homogenize / len(grp))
            for i in range(homogenize_coefficient):
                homogenize_zeros.append(grp[i])
            homogenize_ones.append(grp[homogenize_coefficient])

        if homogenize > 0:
            raise ValueError('homogenize greater than maximum number of homogenization options')

        print('Homogenize zeros:', homogenize_zeros, 'ones:', homogenize_ones)

        homogenize_zero_indices = tuple(coeff_vars.index(var) for var in homogenize_zeros if var in coeff_vars)
        homogenize_one_indices = tuple(coeff_vars.index(var) for var in homogenize_ones if var in coeff_vars)

        # These are the Jacobian matrices of the first derivatives of the extra functions.
        # The equations are simple (either v=0 or v=1), so the derivatives have a very simple (and constant) form.
        # They get precomputed here and will be appended to the Jacobian matrix in jac_fns

        zero_terms = [np.array([1 if i==j else 0 for i in range(len(coeff_vars))], ndmin=2) for j in homogenize_zero_indices]
        one_terms = [np.array([1 if i==j else 0 for i in range(len(coeff_vars))], ndmin=2) for j in homogenize_one_indices]
        homogenize_derivatives = np.vstack(zero_terms + one_terms)
    else:
        homogenize_zero_indices = tuple()
        homogenize_one_indices = tuple()
        homogenize_derivatives = np.empty(shape=(0, len(coeff_vars)))

    # optimize.root methods:
    # 'hybr' (the default) requires same num of eqns as vars
    # 'lm' uses a QR factorization of the Jacobian, then the Levenberg–Marquardt line search algorithm
    # the others uses various approximations to the Jacobian

    SciMin = scipy.optimize.root(fns, iv, jac=jac_fns, method='lm')

    print()
    print()

    if SciMin.success:
        sum_of_squares = sum(map(lambda x: x.get(), [cc.sum_of_squares(SciMin.x) for cc in ccs]))
        print("Sum of squares", sum_of_squares)
        if sum_of_squares < 1e-10:
            print()
            printv(SciMin.x)
            raise ValueError('solution found')
    else:
        print( SciMin.message)

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

#if platform.node() == 'samsung' and not 'no_init' in globals():
#
#    prep_hydrogen()
#    multi_init()
#    multi_expand()

#if not 'no_init' in globals() and not 'mc' in globals():
#    timefunc(prep_helium)
#    timefunc(create_polynomial_ring)

if platform.node() == 'c200-1' and not 'no_init' in globals() and not 'mc' in globals() and len(glob.glob('matrix-*.pickle')) > 0:
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


# Functions to encode and decode deglex exponents

from sage.functions.other import binomial
def choose_with_replacement(setsize,num):
    return binomial(setsize + num - 1, num)

def encode_deglex(exps):
    # modified Gastineau algorithm to produce graded lexicographic order
    #    order within each graded block is reversed from Gastineau's paper
    delta = sum(exps)
    retval = sum(choose_with_replacement(len(exps), j) for j in range(0,delta+1)) - 1
    d = delta
    for i in range(0,len(exps)-1):
        retval -= sum(choose_with_replacement(len(exps)-i-1, d-j) for j in range(0,exps[i]))
        d = d - exps[i]
    return retval

def decode_deglex(ind, len_exps):
    total_degree = 0
    exps = []
    while ind >= choose_with_replacement(len_exps, total_degree):
        ind -= choose_with_replacement(len_exps, total_degree)
        total_degree += 1
    d = total_degree
    ind -= choose_with_replacement(len_exps, total_degree)
    for i in range(0, len_exps-1):
        this_exp = 0
        while ind < -choose_with_replacement(len_exps-i-1, d-this_exp):
            ind += choose_with_replacement(len_exps-i-1, d-this_exp)
            this_exp += 1
        exps.append(this_exp)
        d -= this_exp
    exps.append(d)
    return exps

# Function to compress files being written out uncompressed by the substitution code

import os
import time
import threading
import subprocess

compressor_limit = 2
compressors_available = threading.Semaphore(compressor_limit)
compressors_running = ['radii-526340.out']
compressors_keep_running = True

def compress_file(fn):
    subprocess.run(['gzip', fn])
    compressors_available.release()

def run_compressors():
    while compressors_keep_running:
        candidates = [(os.stat(fn).st_ctime, fn) for fn in os.listdir() if fn.startswith('radii-') and fn.endswith('.out') and fn not in compressors_running]
        if len(candidates) > 12:
            candidates.sort()
            next_candidate = candidates[0][1]
            compressors_available.acquire()
            print("Compressing", next_candidate, "file size", os.stat(next_candidate).st_size, file=sys.stderr)
            compressors_running.append(next_candidate)
            th = threading.Thread(target = compress_file, args = (next_candidate,))
            th.start()
        else:
            print('Sleeping 60 seconds', file=sys.stderr)
            time.sleep(60)

def run_compressors_bg():
    global main_compressor_th
    main_compressor_th = threading.Thread(target = run_compressors)
    main_compressor_th.start()

# Function to run functions in background.
#
# They need to release the GIL for this to be useful.

bg_threads = []

def bg(*args, **kwargs):
    th = threading.Thread(target=args[0], args=args[1:], kwargs=kwargs)
    bg_threads.append(th)
    th.start()
