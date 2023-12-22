# -*- mode: python -*-
#
# Python code to search for solutions of Hydrogen and Helium
# non-relativistic time-invariant Schrodinger equation
#
# INTERACTIVE USAGE:
#
# load('helium.sage')
# prep_hydrogen(5)
# init()
# random_numerical(homogenize=-1)
#
# This will produce a solution to the hydrogen atom using ansatz 5.
# You can also prep_helium(), and supply an ansatz number as argument,
# as in prep_helium(-7).
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
# After init(), instead of calling random_numerical(), you can also
# try an exact calculation:
#
# ideal(eqns_RQQ).radical().primary_decomposition()
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
# by Brent Baccala
#
# first version - August 2019
# latest version - April 2023
#
# no rights reserved; you may freely copy, modify, or distribute this
# program
#
# TODO list:
# - check collected coefficient polynomials to see if they factor
# - automate finding polynomial relations
# - save checkpoints of optimization iterations

import random

import itertools
from itertools import *

import numpy as np

import scipy.optimize

import threading

from sage.symbolic.operators import add_vararg, mul_vararg

from sage.rings.polynomial.polydict import ETuple

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

# from python docs
def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)

# from python docs
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def trial_polynomial(base, coordinates, roots, degree, homogenize=None, constant=True):
    """trial_polynomial(base, coordinates, roots, degree, homogenize=None, constant=True)
    Form a trial polynomial in the Symbolic Ring

    base is a string to which we append numbers to get our coefficient names; i.e, 'a' -> (a0,a1,a2,...)
    coordinates is a tuple of symbolic expressions (currently all symbols; i.e, x1.is_symbol() == True)
    roots is a tuple of symbolic expressions for our roots (currently all powers; i.e, r.operator() == pow) (and all square roots)
    degree is maximum degree of the trial polynomial
    constant=None is optional and drops the constant term (essentially mindeg=1 instead of mindeg=0)
    homogenize=N is unused right now and is intended as a future performance optimization
    """

    if not constant:
        mindegree = 1
    else:
        mindegree = 0
    terms = []
    for deg in range(mindegree, degree+1):
        terms += list(map(mul, (x for x in combinations_with_replacement(roots + coordinates, deg) if all(x.count(r) < 2 for r in roots))))

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
    global A,B,C,D,F,G,M,N,V
    global homogenize_groups
    global alg_exts
    global gamma

    (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
    (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
    (Cvars, C) = trial_polynomial('c', coordinates, roots, 1)
    (Dvars, D) = trial_polynomial('d', coordinates, roots, 1)
    (Fvars, F) = trial_polynomial('f', coordinates, roots, 1)
    (Gvars, G) = trial_polynomial('g', coordinates, roots, 1)

    SR_function = sage.symbolic.function_factory.function

    alg_exts = tuple()

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
        # subs is a list of dictionaries defining substitutions.  They are executed in order.
        subs = [{DD[0](Phi)(B) : Phi(B), DD[0,0](Phi)(B) : Phi(B)},
                {Phi(B) : SR.var('Phi')}
        ]
        homogenize_groups = (Avars, Bvars)
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
        subs = [{DD[0](Xi)(C) : 1/C, DD[0,0](Xi)(C) : -1/C^2},
                {Xi(C) : SR.var('Xi')}
        ]
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

        subs = [{DD[0,0](Chi)(B) : (D/C * DD[0](Chi)(B) + F/C * Chi(B)) + G/C},
                {Chi(B) : SR.var('Chi'), DD[0](Chi)(B) : SR.var('DChi')}
        ]
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

        subs = [{DD[0,0](Chi)(B) : (D/C * DD[0](Chi)(B) + F/C * Chi(B)) + G/C},
                {Chi(B) : SR.var('Chi'), DD[0](Chi)(B) : SR.var('DChi')}
        ]
        ODE_vars = ('Chi', 'DChi')

    elif ansatz == 5:
        # A second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, 1, constant=None)
        Psi = Zeta(V)
        (Avars, A) = trial_polynomial('a', [V], [], 1)
        (Bvars, B) = trial_polynomial('b', [V], [], 1)
        (Cvars, C) = trial_polynomial('c', [V], [], 1)

        homogenize_groups = (Avars, Vvars)

        coeff_vars = (E,) + Vvars + Avars + Bvars + Cvars

        subs = [{DD[0,0](Zeta)(V) : (B * DD[0](Zeta)(V) + C * Zeta(V)) / A},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 5.1:
        # A second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is a quadratic polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, 2, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], 1)
        (Mvars, M) = trial_polynomial('m', [V], [], 1)
        (Nvars, N) = trial_polynomial('n', [V], [], 1)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 5.2:
        # A second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are quadratic polynomials in V, which is a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, 1, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], 2)
        (Mvars, M) = trial_polynomial('m', [V], [], 2)
        (Nvars, N) = trial_polynomial('n', [V], [], 2)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 5.3:
        # A second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are quadratic polynomials in V, which is also a quadratic polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, 2, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], 2)
        (Mvars, M) = trial_polynomial('m', [V], [], 2)
        (Nvars, N) = trial_polynomial('n', [V], [], 2)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
        ]
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

        subs = [{DD[0,0](Zeta)(B/C) : (M * DD[0](Zeta)(B/C) + N * Zeta(B/C)) / D},
                {Zeta(B/C) : SR.var('Zeta'), DD[0](Zeta)(B/C) : SR.var('DZeta')}
        ]
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

        subs = [{DD[0,0](Zeta)(B/C) : (M * DD[0](Zeta)(B/C) + N * Zeta(B/C)) / D},
                {Zeta(B/C) : SR.var('Zeta'), DD[0](Zeta)(B/C) : SR.var('DZeta')}
        ]
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
        subs = [{DD[0](Zeta)(B) : (N * Zeta(B)) / M,
                 DD[0,0](Zeta)(B) : (n1 * Zeta(B) * M + N * N * Zeta(B) - N * Zeta(B) * m1 ) / (M*M) },
                {Zeta(B) : SR.var('Zeta')}
        ]
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
        subs = [{DD[0](Zeta)(B) : n0 * Zeta(B),
                 DD[0,0](Zeta)(B) : n0^2 * Zeta(B)},
                {Zeta(B) : SR.var('Zeta')}
        ]
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

        subs = [{DD[0,0](Zeta)(B) : (M * DD[0](Zeta)(B) + N * Zeta(B)) / D},
                {Zeta(B) : SR.var('Zeta'), DD[0](Zeta)(B) : SR.var('DZeta')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

    elif ansatz == 11:
        # A second-degree algebraic extension (linear coeffs) followed by
        # a second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant

        (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
        (Cvars, C) = trial_polynomial('c', coordinates, roots, 1)
        def deriv(self, *args,**kwds):
            #print("{} {} {}".format(self, args, kwds))
            wrt = args[kwds['diff_param']]
            return -(diff(A, wrt)*self(*coordinates)^2+diff(B,wrt)*self(*coordinates)+diff(C,wrt)/(2*A*self(*coordinates)+B))
        # anything that isn't constant w.r.t. coordinates is an SR_function
        gamma = SR_function('g', nargs=len(coordinates), derivative_func=deriv)

        # We can construct derivatives like this, too:
        # sage: DD[0](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), x1)
        # sage: DD[1](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), y1)
        # sage: DD[1,1](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), y1, y1)

        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots + (gamma(*coordinates),), 1, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], 1)
        (Mvars, M) = trial_polynomial('m', [V], [], 1)
        (Nvars, N) = trial_polynomial('n', [V], [], 1)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars + Avars + Bvars + Cvars
        print(coeff_vars)

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')},
                {gamma(*coordinates) : SR.var('g')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

        alg_exts = (('g', A*gamma(*coordinates)^2 + B*gamma(*coordinates) + C, subs[2]),)

    elif int(ansatz) == 12:
        # A second-degree homogenous ODE (linear coeffs) followed by another
        # second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant

        # trial_polynomial returns a tuple: the coefficient variables used, and the polynomial itself

        if ansatz == 12:
            maxdeg_v = 1
            maxdeg_u = 1
            maxdeg_ode_v = 1
            maxdeg_ode_u = 1
        elif ansatz == 12.1:
            maxdeg_v = 2
            maxdeg_u = 1
            maxdeg_ode_v = 1
            maxdeg_ode_u = 1
        elif ansatz == 12.2:
            maxdeg_v = 1
            maxdeg_u = 2
            maxdeg_ode_v = 1
            maxdeg_ode_u = 1
        elif ansatz == 12.3:
            maxdeg_v = 1
            maxdeg_u = 1
            maxdeg_ode_v = 2
            maxdeg_ode_u = 1
        elif ansatz == 12.4:
            maxdeg_v = 1
            maxdeg_u = 1
            maxdeg_ode_v = 1
            maxdeg_ode_u = 2

        Theta = SR_function('Theta')
        (Uvars, U) = trial_polynomial('u', coordinates, roots, maxdeg_u, constant=None)
        (Avars, A) = trial_polynomial('a', [U], [], maxdeg_ode_u)
        (Bvars, B) = trial_polynomial('b', [U], [], maxdeg_ode_u)
        (Cvars, C) = trial_polynomial('c', [U], [], maxdeg_ode_u)

        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates + (Theta(U),), roots, maxdeg_v, constant=None)

        (Dvars, D) = trial_polynomial('d', [V], [], maxdeg_ode_v)
        (Mvars, M) = trial_polynomial('m', [V], [], maxdeg_ode_v)
        (Nvars, N) = trial_polynomial('n', [V], [], maxdeg_ode_v)

        # Psi is the solution to the PDE
        Psi = Zeta(V)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Uvars + Avars + Bvars + Cvars + Vvars + Dvars + Mvars + Nvars
        print(coeff_vars)

        # Make sure Zeta(V)->SR(Zeta) comes before Theta(U)->SR(Theta) because
        # Zeta(V) will have Theta(U) in its arguments and the Theta(U) sub screws up that match
        subs = [{DD[0,0](Theta)(U) : (B * DD[0](Theta)(U) + C * Theta(U)) / A},
                {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')},
                {Theta(U) : SR.var('Theta'), DD[0](Theta)(U) : SR.var('DTheta')}
        ]

        # It does need a list of the additional variables for when it builds the polynomial ring
        # Probably it could just infer this information from the variables present in the equation
        ODE_vars = ('Theta', 'DTheta', 'Zeta', 'DZeta')

    elif int(ansatz) == 13:
        # A second-degree algebraic extension (linear coeffs) used as the coefficient ring
        # in a second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V (a linear polynomial) and gamma.
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant

        if ansatz == 13:
            maxdeg_v = 1
            maxdeg_alg = 1
            maxdeg_ode = 1
        elif ansatz == 13.1:
            maxdeg_v = 2
            maxdeg_ode = 1
            maxdeg_alg = 1
        elif ansatz == 13.2:
            maxdeg_v = 1
            maxdeg_ode = 2
            maxdeg_alg = 1
        elif ansatz == 13.3:
            maxdeg_v = 1
            maxdeg_ode = 1
            maxdeg_alg = 2
        elif ansatz == 13.4:
            maxdeg_v = 2
            maxdeg_ode = 2
            maxdeg_alg = 1
        elif ansatz == 13.5:
            maxdeg_v = 2
            maxdeg_ode = 1
            maxdeg_alg = 2
        elif ansatz == 13.5:
            maxdeg_v = 1
            maxdeg_ode = 2
            maxdeg_alg = 2
        elif ansatz == 13.6:
            maxdeg_v = 2
            maxdeg_ode = 2
            maxdeg_alg = 2

        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, maxdeg_v, constant=None)

        (Avars, A) = trial_polynomial('a', [V], [], maxdeg_alg)
        (Bvars, B) = trial_polynomial('b', [V], [], maxdeg_alg)
        (Cvars, C) = trial_polynomial('c', [V], [], maxdeg_alg)
        def deriv(self, *args,**kwds):
            #print("{} {} {}".format(self, args, kwds))
            return -(diff(A, V)*self(*coordinates)^2+diff(B,V)*self(*coordinates)+diff(C,V)/(2*A*self(*coordinates)+B))
        # anything that isn't constant w.r.t. coordinates is an SR_function
        gamma = SR_function('g', nargs=1, derivative_func=deriv)

        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [gamma(V)], maxdeg_ode)
        (Mvars, M) = trial_polynomial('m', [V], [gamma(V)], maxdeg_ode)
        (Nvars, N) = trial_polynomial('n', [V], [gamma(V)], maxdeg_ode)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars + Avars + Bvars + Cvars
        print(coeff_vars)

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')},
                {gamma(V) : SR.var('g')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

        #alg_exts = (('g', A*gamma(V)^2 + B*gamma(V) + C, post2_subs),)
        alg_exts = (('g', A*SR.var('g')^2 + B*SR.var('g') + C, subs[2]),)

    elif int(ansatz) == 16:
        # A second-degree algebraic extension (root of a linear polynomial) followed by
        # a second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant

        if ansatz == 16:
            maxdeg_v = 1
            maxdeg_alg = 1
            maxdeg_ode = 1
        elif ansatz == 16.1:
            maxdeg_v = 2
            maxdeg_ode = 1
            maxdeg_alg = 1
        elif ansatz == 16.2:
            maxdeg_v = 1
            maxdeg_ode = 2
            maxdeg_alg = 1
        elif ansatz == 16.3:
            maxdeg_v = 1
            maxdeg_ode = 1
            maxdeg_alg = 2
        elif ansatz == 16.4:
            maxdeg_v = 2
            maxdeg_ode = 2
            maxdeg_alg = 1
        elif ansatz == 16.5:
            maxdeg_v = 2
            maxdeg_ode = 1
            maxdeg_alg = 2
        elif ansatz == 16.5:
            maxdeg_v = 1
            maxdeg_ode = 2
            maxdeg_alg = 2
        elif ansatz == 16.6:
            maxdeg_v = 2
            maxdeg_ode = 2
            maxdeg_alg = 2

        (Avars, A) = trial_polynomial('a', coordinates, roots, maxdeg_alg)
        def deriv(self, *args,**kwds):
            wrt = args[kwds['diff_param']]
            return diff(A, wrt)/(2*A*self(*coordinates))
        # anything that isn't constant w.r.t. coordinates is an SR_function
        gamma = SR_function('g', nargs=len(coordinates), derivative_func=deriv)

        # We can construct derivatives like this, too:
        # sage: DD[0](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), x1)
        # sage: DD[1](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), y1)
        # sage: DD[1,1](gamma)(x1,y1,z1)
        # diff(g(x1, y1, z1), y1, y1)

        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots + (gamma(*coordinates),), maxdeg_v, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], maxdeg_ode)
        (Mvars, M) = trial_polynomial('m', [V], [], maxdeg_ode)
        (Nvars, N) = trial_polynomial('n', [V], [], maxdeg_ode)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars + Avars
        print(coeff_vars)

        subs = [{DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D},
                {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')},
                {gamma(*coordinates) : SR.var('g')}
        ]
        ODE_vars = ('Zeta', 'DZeta')

        # alg_exts is a list of tuples
        # each tuple is (name, minimal polynomial, substitution)
        # minimal polynomial has its roots converted to rs, then subsitution is applied, then append to ideal to mod out by
        alg_exts = (('g', gamma(*coordinates)^2 - A, subs[2]),)

    else:
        raise 'Bad ansatz'

    eq = H(Psi) - E*Psi

    for sub in subs:
        eq = eq.subs(sub)

    # reduce coeff_vars to those which actually appear in the equation
    # let's not do this, in case we've got algebraic extension elements (like ansatz 11)
    # coeff_vars = tuple(sorted(set(eq.free_variables()).intersection(coeff_vars), key=lambda x:str(x)))

    # I used to do this in convert_eq_a(), but that function can be slow, and this is pretty quick,
    # so let's put it in the "prep()" step instead of the "init()" step
    create_polynomial_rings(alg_exts)
    create_eq_a()

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

        # eq (5) in Nakashima and Nakatusji, Solving the Schrodinger equation for helium...
        # THE JOURNAL OF CHEMICAL PHYSICS 127, 224104 2007
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

def create_polynomial_rings(alg_exts):
    global convertRing,idealRing,reduceRing,RQQ,RZZflint,R32003,convertField,num_rvars,num_cvars
    # we need to add gamma to this to make ansatz 11 (algebraic extension) work
    roots_names = list(map(varName, roots))
    alg_exts_names = [p[0] for p in alg_exts]
    num_rvars = len(alg_exts_names) + len(roots_names) + len(ODE_vars) + len(coordinates)
    num_cvars = len(coeff_vars)
    # encoding isn't used except for writing FLINT polynomials to disk with my custom code
    encoding = 'deglex64({}),deglex64({}),sint64'.format(num_rvars, num_cvars)

    # the ordering here is intended to make reduction mod alg_exts and roots easy
    Rsingular = PolynomialRing(QQ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))),
                         order=f'lex({len(alg_exts) + len(roots_names)}), degrevlex({len(ODE_vars) + len(coordinates) + len(coeff_vars)})')

    # used with my custom option to set disk encoding
    #R = PolynomialRing(ZZ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))),
    #                   implementation="FLINT", order='lex', encoding=encoding)
    try:
        Rflint = PolynomialRing(ZZ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))),
                                implementation="FLINT", order='lex')
    except Exception as ex:
        print(ex)
        print('multivariate FLINT rings unavailable')
        Rflint = None

    # not only might FLINT be unavailable, it doesn't implement Groebner bases, so can't be used for reduction

    if Rflint:
        print('Using FLINT rings for convertion')
        convertRing = Rflint
    else:
        print('Using Singular rings for convertion')
        convertRing = Rsingular

    print('Using Singular rings for reduction ideal')
    idealRing = Rsingular

    #if len(roots) > 0 or len(alg_exts) > 0:
    if False:
        print('Using Singular rings for reduction')
        reduceRing = Rsingular

        # This doesn't work - Singular interface doesn't support splitting variables with the coeff field like this
        # (according to the first few lines of multi_polynomial_libsingular.pyx)
        # I'm not sure about Singular proper.  It has some kind of support for "transcendental extension of Q"
        #Rsingular1 = PolynomialRing(QQ, names=tuple(flatten((ODE_vars, coordinates, coeff_vars))))
        #Rsingular2 = PolynomialRing(FractionField(Rsingular1), names=tuple(flatten((alg_exts_names, roots_names))), order='lex')
    else:
        print('Using convertion ring for reduction')
        reduceRing = convertRing

    convertField = Frac(convertRing)

    # These are the rings used for the system of equations in the coefficients
    RQQ = PolynomialRing(QQ, names=coeff_vars)
    RZZflint = PolynomialRing(ZZ, names=coeff_vars, implementation='FLINT')
    R32003 = PolynomialRing(GF(32003), names=coeff_vars)

# we need to add gamma to this to make ansatz 11 (algebraic extension) work
def mk_ideal(R, roots, alg_exts):
    "Given a list or tuple of roots, return a ideal of ring R that reduces the global variable names of those roots"
    global reductionIdeal
    # We expect a tuple of pow's in the Symbolic Ring, so we can easily construct the minimal polynomials
    #
    # We can also take a pair of (varName, minpoly) where varName is a var in the Symbolic Ring
    # and minpoly is a polynomial in the Symbolic Ring that can be converted to R
    ideal_generators = []
    for v in roots:
        assert v.operator() is operator.pow
        Rname = R(varName(v))
        Rexpr = R(v.operands()[0])
        power = int(1/v.operands()[1])
        ideal_generators.append(Rname^power - Rexpr)
    for v,e,postsub in alg_exts:
        assert e in SR
        ideal_generators.append(R(roots_to_rs(e).subs(postsub)))
    reductionIdeal = ideal(ideal_generators)

def convert_eq_a():
    global eq_a_convertField
    # If we write this as 'eq_a_convertField = convertField(eq_a)', Sage will attempt to construct eq_a_convertField by calling
    # eq_a.numerator() and eq_a.denominator(), which will perform lots of rational function
    # math in the Symbolic Ring, which is very slow and memory intensive.  Calling it
    # like 'eq_a.polynomial(ring=convertField)' recurses through the expression tree and builds the
    # expression from the bottom up using polynomial ring operations, which are much more efficient.
    #
    # This trick (currently) only works on my development Sage, so try it and fall back on the slower way.
    try:
        eq_a_convertField = eq_a.polynomial(ring=convertField)
    except:
        print('WARNING: converting eq_a using the Symbolic Ring (this is slow)')
        eq_a_convertField = convertField(eq_a)
    print('eq_a_convertField numerator:', eq_a_convertField.numerator().number_of_terms(), 'terms')

def convertRing_to_reduceRing(element):
    # this doesn't work with my current development sage:
    #   eq_a_reduceRing_n = reduceRing(eq_a_convertField.numerator()).mod(I)
    #   eq_a_reduceRing_d = reduceRing(eq_a_convertField.denominator()).mod(I)
    # I tried converting to a string, but that hits "RecursionError: maximum recursion depth exceeded during compilation"
    #   eq_a_reduceRing_n = reduceRing(str(eq_a_convertField.numerator())).mod(I)
    #   eq_a_reduceRing_d = reduceRing(str(eq_a_convertField.denominator())).mod(I)
    # go this way instead: (works on a simple reduceRing)
    if convertRing != reduceRing:
        return reduceRing(element.dict())
    else:
        return element
    # if reduceRing is a ring over a field with convertRing variables split between the two, we need something else
    # don't bother with this code, as Singular can't support splitting variables between the ring and the coeff field
    #
    #baseRing = reduceRing.base_ring()
    #rest_term_sub = {convertRing(v):1 for v in reduceRing.variable_names()}
    #result = reduceRing(0)
    #for coeff, monomial in element:
    #    coeff_term = monomial.subs(rest_term_sub)
    #    ring_term = reduceRing(monomial / coeff_term)
    #    result += coeff * baseRing(coeff_term) * ring_term
    #return result

def reduce_mod_ideal(element, I=None):
    if I:
        # This is slower if convertRing is FLINT and reduceRing is Singular; it does the reduction in Singular:
        #    return convertRing_to_reduceRing(element).mod(I)
        # This way does the reduction in FLINT:
        # if type(reduceRing) == sage.rings.polynomial.multi_polynomial_flint.MPolynomialRing_flint:
        for p in I.groebner_basis():
            element %= convertRing(str(p))
        return convertRing_to_reduceRing(element)
    else:
        return convertRing_to_reduceRing(element)

def reduce_numerator(I=None):
    global eq_a_reduceRing_n
    eq_a_reduceRing_n = reduce_mod_ideal(eq_a_convertField.numerator(), I)
    print('eq_a_reduceRing_n:', eq_a_reduceRing_n.number_of_terms(), 'terms')

def reduce_denominator(I=None):
    global eq_a_reduceRing_d
    eq_a_reduceRing_d = reduce_mod_ideal(eq_a_convertField.denominator(), I)
    print('eq_a_convertField: denominator', eq_a_reduceRing_d.number_of_terms(), 'terms')

import time

def timefunc(func, *args, **kwargs):
    start_time = time.perf_counter()
    retval = func(*args, **kwargs)
    end_time = time.perf_counter()
    print('{:30} {:10.2f} sec'.format(func.__name__, end_time - start_time))
    return retval

# Now expand out powers, and collect like x,y,z's terms together to
# get a system of polynomials
#
# This is a slow step, so I've tried several different ways to do it.

# Look for solutions using an approximate numerical technique

last_time = 0

# We're going from reduceRing to whatever ring is specified, or reduceRing (if not specified)

def build_system_of_equations(ring=None):
    result = dict()
    # for speed, build this tuple here instead of letting the subs method do it in monomial.subs
    # if my custom __evaluate function is available, use it, it's yet faster
    if '__evaluate' in dir(eq_a_reduceRing_n):
        FLINT_evaluate = tuple((n,1) for n in range(reduceRing.ngens()) if reduceRing.gen(n) in coeff_vars)
        print('Using FLINT __evaluate method')
    else:
        FLINT_evaluate = None
        non_coeff_sub = tuple(1 if reduceRing.gen(n) in coeff_vars else reduceRing.gen(n) for n in range(reduceRing.ngens()))
    pb = ProgressBar(label='build_system_of_equations ', expected_size=eq_a_reduceRing_n.number_of_terms())
    # this loop works on Singular or FLINT elements, but not other things like rings with variables in their coeff field
    for i, (coeff, monomial) in enumerate(eq_a_reduceRing_n):
        if i%100 == 99:
            pb.show(i+1)
        if FLINT_evaluate:
            non_coeff_part = monomial.__evaluate(FLINT_evaluate)
        else:
            non_coeff_part = monomial(non_coeff_sub)
        if ring:
            coeff_part = ring(monomial / non_coeff_part)
        else:
            # this cast needs to be here because otherwise the division (even though it's exact) takes us to the fraction field
            coeff_part = reduceRing(monomial / non_coeff_part)
        if (non_coeff_part) in result:
            result[non_coeff_part] += coeff * coeff_part
        else:
            result[non_coeff_part] = coeff * coeff_part
    pb.show(i+1)
    pb.done()
    return tuple(set(result.values()))

def create_eqns_RQQ():
    global eqns_RQQ, jac_eqns_RQQ
    eqns_RQQ = timefunc(build_system_of_equations, RQQ)

def create_jac_eqns_RQQ():
    global eqns_RQQ, jac_eqns_RQQ
    jac_eqns_RQQ = [[diff(eqn, RQQ(v)) for eqn in eqns_RQQ] for v in coeff_vars]

def create_eqns_R32003():
    global eqns_R32003
    eqns_R32003 = tuple(map(lambda arg: arg.map_coefficients(GF(32003), GF(32003)), eqns_RQQ))

def init():
    global reductionIdeal
    # convert_eq_a is the first really time consuming step
    timefunc(convert_eq_a)
    if len(roots) > 0 or len(alg_exts) > 0:
        timefunc(mk_ideal, idealRing, roots, alg_exts)
    else:
        reductionIdeal = None
    timefunc(reduce_numerator, reductionIdeal)
    # We don't use the denominator for anything, currently
    # timefunc(reduce_denominator, reductionIdeal)
    timefunc(create_eqns_RQQ)
    timefunc(create_jac_eqns_RQQ)
    timefunc(create_eqns_R32003)

def bertini(eqns=None, file=sys.stdout):
    if not eqns:
        eqns = eqns(RQQ)
    print("""
CONFIG

TRACKTYPE: 1;

END;

INPUT
    """, file=file)
    print('variable_group ', ', '.join( [str(var) for var in set().union(*[eqn.variables() for eqn in eqns])]), ';', file=file)
    print("function ", ', '.join([f"f{i}" for i in range(1,len(eqns)+1)]), ';', file=file)
    for i,eqn in enumerate(eqns):
        print(f"f{i+1} = ", eqn, ";", file=file)
    print('END;', file=file)

def cocoa_dump(fn, eqns, elim=None):
    with open(fn, 'w') as f:
        print('use R ::= QQ[', ','.join(str(g) for g in eqns[0].parent().gens()), ']', ',' + elim if elim else '', ';', file=f)
        print('I := ideal(', ','.join(str(e) for e in eqns if e != 0), ');', file=f)

# Euclidean Distance polynomials are used to construct a masking function to drive our
# numerical solution algorithm away from known varieties.

def construct_ED_polynomial(I):
    """Construct and return the Euclidean Distance polynomial for a given ideal"""
    uvars = tuple('u_' + str(var) for var in coeff_vars)
    global RED, REDflint
    RED = PolynomialRing(QQ, names=tuple(flatten((coeff_vars, uvars, ('t',)))),
                         order=f'degrevlex({len(coeff_vars)}), degrevlex({len(coeff_vars) + 1})')
    REDflint = PolynomialRing(ZZ, names=tuple(flatten((coeff_vars, uvars, ('t',)))), implementation='FLINT')
    t = RED('t')
    codim = RQQ.ngens() - I.dimension()
    # compute singular locus, generated by all codim-minors of Jac(I)
    # this assert is the easy case.  do it first
    assert codim == I.ngens()
    #Jacobian = matrix([[g.derivative(RQQ(cv)) for g in I.gens()] for cv in coeff_vars])
    singular_locus_generators = [RED(0), RED(1)]
    # binomial(39,9) = 211915132, so let's just do it by hand right now
    # this is slow because it computes the derivatives in the inner loop, when they don't change over combinations
    #for subset_coeff_vars in combinations(coeff_vars, codim):
    #    M = matrix([[g.derivative(RQQ(cv)) for g in I.gens()] for cv in subset_coeff_vars])
    #    singular_locus_generators.append(M.determinant())
    singular_locus = ideal(singular_locus_generators)

    if max(g.degree() for g in I.gens()) == 1:
        # RED(g) doesn't work right because the variable names are different; need RED(str(g))
        # this version of critical_ideal_generators only works if degree(I) == 1
        critical_ideal_generators = [RED(str(g)) for g in I.gens()] + [RED(cv) - RED('u_'+str(cv)) for cv in coeff_vars if RQQ(cv) not in I]
    else:
        Jacobian = matrix([[g.derivative(RQQ(cv)) for g in I.gens()] for cv in coeff_vars])
        JacobianFLINT = Jacobian.apply_map(lambda x: REDflint(x) if x.is_constant() else REDflint(str(x)))
        extendedJacobianFLINT = JacobianFLINT.transpose().stack(matrix([REDflint('u_' + str(cv)) - REDflint(cv) for cv in coeff_vars]))
        # This next step can be quite slow, as there are binomial(nvars, codim+1) minors.  Use FLINT for speed.
        minors = set(extendedJacobianFLINT.minors(codim+1))
        # convert to Singular (RED) because we're about to do a Groebner basis calculation
        critical_ideal_generators = [RED(str(g)) for g in chain(minors, I.gens())]

    return critical_ideal_generators

    distance_function = t^2 - sum((RED(cv) - RED('u_'+str(cv)))^2 for cv in coeff_vars)
    elimination_ideal = ideal(critical_ideal_generators + [distance_function])
    return elimination_ideal.groebner_basis()[-1]

ED_ideals = []
ED_polynomials = []

#with open('giac.script', 'w') as f:
#    print('greduce(', str(distance_function), ',gbasis(',str(bwb), ',', list(RED.gens()), '))', file=f)

#bwb = construct_ED_polynomial(I)
#with open('giac.script', 'w') as f:
#    print('gbasis(',str(list(bwb + [distance_function])), ',', list(RED.gens()), ',plex)[-1]', file=f)

def masking_function(x):
    if x < 2:
        return exp((1-x)/(2*x-x^2))
    else:
        return 0

def masking_function_derivative(x):
    if x < 2:
        return (-x^2+2*x-2)/(2*x-x^2)^2 * exp((1-x)/(2*x-x^2))
    else:
        return 0

def ED_function(v, EDpoly):
    ED_eqn = EDpoly.subs(dict(zip(map(lambda x: RED('u_' + str(x)), coeff_vars), v)))
    solset_symbolic = solve([SR(ED_eqn) == 0], SR('t'), solution_dict=True)
    solset_numerical = [sol[SR('t')].n() for sol in solset_symbolic]
    sol = min(filter(lambda x: x>=0, solset_numerical))
    return sol

# Now we want to evaluate possibly lots of polynomials, as well as
# their first derivatives.  For speed, we use scipy sparse matrices.
# Each row in the matrix is a polynomial and each column corresponds
# to a monomial.  The matrix entry is that monomial's coefficient.  We
# can then evaluate all the polynomials at once by forming a vector of
# all the monomials (up to a given maximum degree) and multiplying it
# by the matrix.

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
#
# Array of n-th degree terms: [terms involving a, terms involving b but not a, terms involving only c]
#
# To generate array of n-th degree terms:
#    - multiply a by all (n-1)-th degree terms to get all n-th degree terms involving a
#    - multiply b by all (n-1)-th degree terms not involving a to get all n-th degree terms involving b but not a
#      etc
#
# stack is an array of arrays.  The second array is special, it's actually a numpy vector,
# not a python array, but the logic still works.
#
# second array (first degree terms):
#   -> [a,b,c]
#
# third array (second degree terms):
#   -> [[a^2,a*b,a*c],[b^2,b*c],[c^2]]
#   -> [a*[a,b,c], b*[b,c], c*[c]]
#
# fourth array (third order terms):
#   -> [[a^3,a^2*b,a^2*c,a*b^2,a*b*c,a*c^2,a*b^2,b^2*c,b*c^2,c^3]
#   -> [a*[a^2,a*b,a*c,b^2,b*c,c^2], b*[b^2,b*c,c^2], c*[c^2]]

def generate_multi_vector(max_degree, v):
    npv = np.array(v)
    stack = [[1], npv]
    for d in range(1, max_degree):
        stack.append([np.hstack(stack[-1][i:]) * npv[i] for i in range(len(v))])
    res = np.hstack(tuple(flatten(stack)))
    return res

# same idea, but for a first derivative

def generate_multi_D_vector(max_degree, v, var):
    ind = coeff_vars.index(var)

    npv = np.array(v)

    firsts = np.zeros(npv.size)
    firsts[ind] = int(1)
    D_stack = [[0], firsts]
    stack = [[1], npv]

    for d in range(1, max_degree):
        D_stack.append([np.hstack(stack[-1][i:]) * firsts[i] + np.hstack(D_stack[-1][i:]) * npv[i]
                        for i in range(len(v))])
        stack.append([np.hstack(stack[-1][i:]) * npv[i] for i in range(len(v))])

    res = np.hstack(tuple(flatten(D_stack)))

    return res

from sage.functions.other import binomial
def choose_with_replacement(setsize,num):
    return binomial(setsize + num - 1, num)

def encode_deglex(exps):
    # Needs to produce indices in the same ordering as generate_multi_vector()
    #
    # modified Gastineau algorithm to produce graded lexicographic order
    #    order within each graded block is reversed from Gastineau's paper
    delta = sum(exps)
    retval = sum(choose_with_replacement(len(exps), j) for j in range(0,delta+1)) - 1
    d = delta
    for i in range(0,len(exps)-1):
        retval -= sum(choose_with_replacement(len(exps)-i-1, d-j) for j in range(0,exps[i]))
        d = d - exps[i]
    return retval

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

def convert_to_matrix(system_of_equations):
    global matrix_RQQ
    # I haven't calculated in advance the maximum degree of the monomials,
    # so start at -1 and reshape the matrix every time we hit a monomial
    # that's too big.
    max_degree = -1
    dok = scipy.sparse.dok_matrix((len(system_of_equations), 0), np.int64)

    pb = ProgressBar(label='convert_to_matrix ', expected_size=len(system_of_equations))
    for i,eqn in enumerate(system_of_equations):
        if eqn.degree() > max_degree:
            # increase max_degree and rebuild indices
            max_degree = eqn.degree()
            # index of the smallest tuple of the next higher degree
            veclen=encode_deglex([max_degree + 1] + [0]*(eqn.parent().ngens() - 1))
            dok.resize((len(system_of_equations), veclen))
        for etuple, coeff in eqn.iterator_exp_coeff():
            index = encode_deglex(etuple)
            dok[i, index] += coeff
        pb.show(i+1)

    pb.show(i+1)
    pb.done()

    matrix_RQQ = sp_unique(dok, axis=0, new_format='csr')
    matrix_RQQ.max_degree = max_degree

use_matrix_code = False

def fns(v):
    r"""
    Evaluate our polynomials at a given coordinate.

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - a numpy vector of all of our polynomials, evaluated at ``vec``

    ALGORITHM:

    - evaluate the polynomials in eqns_RQQ and
      add additional polynomials to enfore "homogenize" conditions.
    """

    start_time = time.time()

    # Save a copy of vector to aid in stopping and restarting the
    # calculation.  An explicit call to the copy method is required if
    # we're using the Fortran minpack code (i.e, scipy's optimize
    # package) because in that case, 'v' is only a pointer to a
    # Fortran array that gets deallocated once the Fortran code exits.
    # (I think - the copy's definitely needed, though)
    global last_v
    last_v = v.copy()

    # not really used anymore; older "homogenization" code
    homogenize_terms = [v[i] for i in homogenize_zero_indices] + [v[i] - 1 for i in homogenize_one_indices]

    ED_terms = []
    for EDpoly in ED_polynomials:
        ED_eqn = EDpoly.subs(dict(zip(map(lambda x: RED('u_' + str(x)), coeff_vars), v)))
        #print(ED_eqns)
        soleqns = solve(SR(ED_eqn) == 0, SR('t'), solution_dict=True)
        #for sol in solve([SR(EDpoly) == 0 for EDpoly in ED_eqns], SR('t'), solution_dict=True): print(sol[SR('t')].n())

        solset = [sol[SR('t')].n() for sol in soleqns]
        sol = min(filter(lambda x: x>=0, solset))
        #print(sol)
        #raise ValueError
        ED_terms.append(masking_function(sol))

    if use_matrix_code:
        multivec = generate_multi_vector(matrix_RQQ.max_degree, v)
        res = np.hstack((matrix_RQQ.dot(multivec), homogenize_terms, ED_terms))
    else:
        # optimized code (only construct dict once):
        substitution = dict(zip(map(RQQ, coeff_vars), v))
        #res = np.hstack([eqn.subs(substitution) for eqn in eqns_RQQ] + homogenize_terms)

        # further optimized code (use call instead of subs):
        call_substitution = tuple(substitution.get(RQQ.gen(n), 0) for n in range(RQQ.ngens()))
        res = np.hstack([eqn(call_substitution) for eqn in eqns_RQQ] + homogenize_terms + ED_terms)

    global last_time
    sum_of_squares = sum(res*res)
    if last_time == 0:
        print(sum_of_squares)
    else:
        #if time.time() - last_time > 10:
        if time.time() - start_time > 10:
            printv(res)
        #print("{:<30} {:10.2f} sec".format(sum_of_squares, time.time()-last_time))
        print("{:<30} {:10.2f} sec".format(sum_of_squares, time.time()-start_time))
    last_time = time.time()
    return res

def jac_fns(vec):
    r"""
    Evaluate the Jacobian matrix (the matrix of first-order
    partial derivatives) of our polynomials

    INPUT:

    - ``vec`` -- a vector of real values for all coeff_vars

    OUTPUT:

    - the Jacobian matrix, evaluated at ``vec``, as a numpy matrix,
    with as many rows as polynomials (plus homogenization conditions)
    and as many columns as coeff_vars
    """

    # compute sol and masking_function(sol)

    global ED_derivatives
    ED_derivatives = [np.zeros((0, len(coeff_vars)))]
    for EDpoly in ED_polynomials:
        global RED_mapper, ED_derivs, t_deriv, coeff_derivs, coeff_derivs_eval

        sol = ED_function(vec, EDpoly)
        if masking_function(sol) == 0:
            ED_derivatives.append(np.array([0 for cv in coeff_vars]))
        else:
            RED_mapper = dict(zip(map(lambda x: RED('u_' + str(x)), coeff_vars), vec))
            RED_mapper[RED('t')] = sol
            ED_derivs = [EDpoly.derivative(RED('u_' + str(cv))).subs(RED_mapper) for cv in coeff_vars]
            t_deriv = sum(c*m.derivative(RED('t')) for c,m in EDpoly)
            #d0_deriv = sum(c*m.derivative(RED('u_d0')) for c,m in ED_polynomials[0])
            coeff_derivs = [-sum(c*m.derivative(RED('u_' + str(cv)))/t_deriv for c,m in EDpoly) for cv in coeff_vars]
            coeff_derivs_eval = [coeff_deriv.subs(RED_mapper) for coeff_deriv in coeff_derivs]
            ED_derivatives.append(np.array([masking_function_derivative(sol)*cde for cde in coeff_derivs_eval]))
    ED_derivatives = np.vstack(ED_derivatives)

    global mdv, dN
    if use_matrix_code:
        pb = ProgressBar(label='jac_fns', expected_size=len(coeff_vars))
        # each dot product is neqns in size; this is a tuple of ngens vectors, each neqns in size
        def dot_and_show(i,var):
            res = matrix_RQQ.dot(generate_multi_D_vector(matrix_RQQ.max_degree, vec, var))
            pb.show(i+1)
            return res
        mdv = tuple(dot_and_show(i,var) for i,var in enumerate(coeff_vars))
        # eqns on the vertical axis, coeffs on the horizontal
        dN = np.vstack((np.stack(mdv, axis=1), homogenize_derivatives, ED_derivatives))
    else:
        mapper = dict(zip(map(RQQ, coeff_vars), vec))
        call_substitution = tuple(mapper.get(RQQ.gen(n), 0) for n in range(RQQ.ngens()))
        dN = np.vstack([np.hstack([np.vstack([eqn(call_substitution) for eqn in eqnblock]) for eqnblock in jac_eqns_RQQ]), homogenize_derivatives, ED_derivatives])
    return dN

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
        v = SciMin.x
        homogenize_terms = [v[i] for i in homogenize_zero_indices] + [v[i] - 1 for i in homogenize_one_indices]
        res = np.hstack([eqn.subs(dict(zip(map(RQQ, coeff_vars), v))) for eqn in eqns_RQQ] + homogenize_terms)
        sum_of_squares = sum(res*res)
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

# Factor all of the polynomials in the system of equations and build a set of systems,
# all with irreducible polynomials, that generates the same variety.  From each factored
# set that arises from a polynomial in the original system, we want one factor.
#
# Currently works, but builds way more systems that are really needed.  871 on hydrogen-5,
# when we know we'll only get 5 irreducible varieties.

# commented out so it won't run on load, but needs to be done in build_systems():

def init_build_systems():
    global factored_eqns, factors, factor_dict, index_dict, indices, findices
    factored_eqns = [factor(eq) for eq in eqns_RQQ]
    factors=list(set(f for factored_eqn in factored_eqns for f,m in factored_eqn))
    factor_dict = {f:i for i,f in enumerate(factors)}
    index_dict = {i:f for i,f in enumerate(factors)}
    # indices is a tuple of tuples, each sub-tuple is a set of factor indexes for a single equation in the original system
    indices = tuple(tuple(factor_dict[f] for f,m in factored_eqn) for factored_eqn in factored_eqns)
    # findices is a tuples of tuples of factors, with the multiplicities dropped
    findices = tuple(tuple(f for f,m in factored_eqn) for factored_eqn in factored_eqns)

def add_factor(f):
    if f not in factors:
        factors.append(f)
        factor_dict[f] = len(factors)-1
        index_dict[len(factors)-1] = f

def build_systems():
    global systems
    systems = set()
    working_ideal = set()
    tracking_info = list()
    substitutions = dict()
    last_i = -1
    while True:
        for i in range(last_i+1, len(findices)):
            # if any index in the working ideal is in this equation, skip it, as it's already satisfied
            if not working_ideal.isdisjoint(findices[i]):
                continue
            working_ideal.add(findices[i][0])
            tracking_info.append((i, 0))
            #print('adding', i, 0)
        # yield working_ideal
        systems.add(tuple(sorted(working_ideal)))
        while True:
            try:
                last_i, last_index = tracking_info.pop()
            except IndexError:
                return
            # remove will raise KeyError if this doesn't work
            try:
                working_ideal.remove(findices[last_i][last_index])
            except KeyError:
                print('KeyError removing from', working_ideal)
                raise
            #print('removing', last_i, last_index)
            if last_index < len(indices[last_i])-1:
                break
        working_ideal.add(findices[last_i][last_index + 1])
        tracking_info.append((last_i, last_index + 1))
        #print('adding', last_i, last_index+1)
