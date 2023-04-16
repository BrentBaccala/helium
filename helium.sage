# -*- mode: python -*-
#
# Python code to search for solutions of Hydrogen and Helium
# non-relativistic time-invariant Schrodinger equation
#
# INTERACTIVE USAGE:
#
# load('helium.sage')
# prep_hydrogen()
# init()
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

    (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
    (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
    (Cvars, C) = trial_polynomial('c', coordinates, roots, 1)
    (Dvars, D) = trial_polynomial('d', coordinates, roots, 1)
    (Fvars, F) = trial_polynomial('f', coordinates, roots, 1)
    (Gvars, G) = trial_polynomial('g', coordinates, roots, 1)

    SR_function = sage.symbolic.function_factory.function

    alg_exts = tuple()

    post2_subs = dict()

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
        # A second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant
        Zeta = SR_function('Zeta')
        (Vvars, V) = trial_polynomial('v', coordinates, roots, 1, constant=None)
        Psi = Zeta(V)
        (Dvars, D) = trial_polynomial('d', [V], [], 1)
        (Mvars, M) = trial_polynomial('m', [V], [], 1)
        (Nvars, N) = trial_polynomial('n', [V], [], 1)

        homogenize_groups = (Dvars, Vvars)

        coeff_vars = (E,) + Vvars + Dvars + Mvars + Nvars

        pre_subs = {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D}
        post_subs = {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
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

        pre_subs = {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D}
        post_subs = {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
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

        pre_subs = {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D}
        post_subs = {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
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

        pre_subs = {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D}
        post_subs = {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
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

    elif ansatz == 11:
        # A second-degree algebraic extension (linear coeffs) followed by
        # a second-order homogeneous ODE: D(V) d^2 Zeta/dV^2 - M(V) dZeta/dV - N(V) Zeta = 0
        # where D(V), M(V), and N(V) are linear polynomials in V, which is itself a linear polynomial
        #
        # Homogenization forces V and D to be non-zero; V is also forced to be non-constant

        global gamma
        (Avars, A) = trial_polynomial('a', coordinates, roots, 1)
        (Bvars, B) = trial_polynomial('b', coordinates, roots, 1)
        (Cvars, C) = trial_polynomial('c', coordinates, roots, 1)
        def deriv(self, *args,**kwds):
            #print("{} {} {}".format(self, args, kwds))
            wrt = args[kwds['diff_param']]
            return -(diff(A, wrt)*self(*coordinates)^2+diff(B,wrt)*self(*coordinates)+diff(C,wrt)/(2*A*self(*coordinates)+B))
        # anything that isn't constant w.r.t. coordinates is an SR_function
        gamma = SR_function('g', nargs=3, derivative_func=deriv)

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

        pre_subs = {DD[0,0](Zeta)(V) : (M * DD[0](Zeta)(V) + N * Zeta(V)) / D}
        post_subs = {Zeta(V) : SR.var('Zeta'), DD[0](Zeta)(V) : SR.var('DZeta')}
        post2_subs = {gamma(*coordinates) : SR.var('g')}
        ODE_vars = ('Zeta', 'DZeta')

        alg_exts = (('g', A*gamma(*coordinates)^2 + B*gamma(*coordinates) + C, post2_subs),)

    else:
        raise 'Bad ansatz'

    eq = H(Psi) - E*Psi

    eq = eq.subs(pre_subs).subs(post_subs).subs(post2_subs)

    # reduce coeff_vars to those which actually appear in the equation
    # let's not do this, in case we've got algebraic extension elements (like ansatz 11)
    # coeff_vars = tuple(sorted(set(eq.free_variables()).intersection(coeff_vars), key=lambda x:str(x)))

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

def create_polynomial_ring(alg_exts):
    global R,RQQ,F,num_rvars,num_cvars
    # we need to add gamma to this to make ansatz 11 (algebraic extension) work
    roots_names = list(map(varName, roots))
    alg_exts_names = [p[0] for p in alg_exts]
    num_rvars = len(alg_exts_names) + len(roots_names) + len(ODE_vars) + len(coordinates)
    num_cvars = len(coeff_vars)
    encoding = 'deglex64({}),deglex64({}),sint64'.format(num_rvars, num_cvars)
    if len(roots) > 0 or len(alg_exts) > 0:
        # FLINT multivariates can't handle reduction modulo an ideal, so use Singular multivariates instead
        print('Using Singular implementation')
        R = PolynomialRing(QQ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))), order='lex')
    else:
        print('Using FLINT implementation')
        R = PolynomialRing(ZZ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))),
                           implementation="FLINT", order='lex', encoding=encoding)
    # I don't want order=lex because this is the ring I'll use for Groebner basis calculations
    RQQ = PolynomialRing(QQ, names=tuple(flatten((alg_exts_names, roots_names, ODE_vars, coordinates, coeff_vars))),
                         order=f'degrevlex({len(alg_exts) + len(roots_names) + len(ODE_vars) + len(coordinates)}), degrevlex({len(coeff_vars)})')
    F = Frac(R)

# we need to add gamma to this to make ansatz 11 (algebraic extension) work
def mk_ideal(R, roots, alg_exts):
    "Given a list or tuple of roots, return a ideal of ring R that reduces the global variable names of those roots"
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
    return ideal(ideal_generators)

def convert_eq_a():
    global F_eq_a, F_eq_a_n, F_eq_a_d
    create_polynomial_ring(alg_exts)
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
        F_eq_a_n = F_eq_a.numerator().mod(mk_ideal(R, roots, alg_exts))
        F_eq_a_d = F_eq_a.denominator().mod(mk_ideal(R, roots, alg_exts))
    else:
        F_eq_a_n = F_eq_a.numerator()
        F_eq_a_d = F_eq_a.denominator()
    print('F_eq_a: numerator', F_eq_a_n.number_of_terms(), 'terms; denominator', F_eq_a_d.number_of_terms(), 'terms')

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

global idealnum
idealnum = 0

def eqns_from_eq_a(ring=None):
    result = dict()
    for coeff, monomial in F_eq_a_n:
        rest_term = monomial.subs({R(v):1 for v in coeff_vars})
        coeff_term = monomial / rest_term
        if (rest_term) in result:
            if ring:
                result[rest_term] += ring(coeff * coeff_term)
            else:
                result[rest_term] += coeff * coeff_term
        else:
            if ring:
                result[rest_term] = ring(coeff * coeff_term)
            else:
                result[rest_term] = coeff * coeff_term
    global eqns
    eqns = tuple(result.values())

    # "Fan out" an ideal by considering all possible states (=0 or !=0)
    # of its generators.  The idea is to divide out by a prime ideal
    # and try to find a witness point for the remaining ideals.

    gennames = 'd0,d1,m0,m1,n0,n1'.split(',')
    # probably this has to be a Groebner basis
    gens = tuple(map(RQQ, gennames))
    def fn2(x,y):
        Q,R = x.quo_rem(y)
        print(Q, '*', y, ' + ', end='')
        return R
    for eqn in eqns:
        print(eqn, ' = ', end='')
        reduce(fn2, (eqn,) + gens)
        print('')
    global zerogens, nonzerogens, nonzeroprod
    zerogens = tuple(gen for i,gen in enumerate(gens) if (idealnum & 2^i) == 0)
    nonzerogens = tuple(gen for i,gen in enumerate(gens) if (idealnum & 2^i) != 0)
    nonzeroprod = mul(gen for i,gen in enumerate(gens) if (idealnum & 2^i) != 0)
    I = ideal(zerogens)
    #return [eqn.mod(I) for eqn in eqns] + list(zerogens)
    #return [eqn.mod(I)/nonzeroprod for eqn in eqns] + list(zerogens)
    return [reduce(fn2, (eqn,) + zerogens)/nonzeroprod for eqn in eqns] + list(zerogens)

def create_eqns_RQQ():
    global eqns_RQQ, jac_eqns_RQQ
    eqns_RQQ = tuple(set(timefunc(eqns_from_eq_a, RQQ)))

def create_jac_eqns_RQQ():
    global eqns_RQQ, jac_eqns_RQQ
    jac_eqns_RQQ = [[diff(eqn, RQQ(v)) for eqn in eqns_RQQ] for v in coeff_vars]

def init():
    timefunc(create_eq_a)
    timefunc(convert_eq_a)
    timefunc(create_eqns_RQQ)
    timefunc(create_jac_eqns_RQQ)

def bertini(eqns=None):
    if not eqns:
        eqns = eqns(RQQ)
    print("""
CONFIG

TRACKTYPE: 1;

END;

INPUT
    """)
    print('variable_group ', ', '.join( [str(var) for var in set().union(*[eqn.variables() for eqn in eqns])]), ';')
    print("function ", ', '.join([f"f{i}" for i in range(1,len(eqns)+1)]), ';')
    for i,eqn in enumerate(eqns):
        print(f"f{i+1} = ", eqn, ";")
    print('END;')

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

    # Save a copy of vector to aid in stopping and restarting the
    # calculation.  An explicit call to the copy method is required if
    # we're using the Fortran minpack code (i.e, scipy's optimize
    # package) because in that case, 'v' is only a pointer to a
    # Fortran array that gets deallocated once the Fortran code exits.
    # (I think - the copy's definitely needed, though)
    global last_v
    last_v = v.copy()

    homogenize_terms = [v[i] for i in homogenize_zero_indices] + [v[i] - 1 for i in homogenize_one_indices]
    res = np.hstack([eqn.subs(dict(zip(map(RQQ, coeff_vars), v))) for eqn in eqns_RQQ] + homogenize_terms)

    global last_time
    sum_of_squares = sum(res*res)
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
    with as many rows as polynomials (plus homogenization conditions)
    and as many columns as coeff_vars
    """

    mapper = dict(zip(map(RQQ, coeff_vars), v))
    # dN = np.vstack(list(map(lambda x: x.get(), [cc.jac_fns(v) for cc in ccs])) + [homogenize_derivatives])
    #dN = np.hstack([np.vstack([eqn.subs(mapper) for eqn in eqnblock] + [homogenize_derivatives]) for eqnblock in jac_eqns_RQQ])
    dN = np.vstack([np.hstack([np.vstack([eqn.subs(mapper) for eqn in eqnblock]) for eqnblock in jac_eqns_RQQ]), homogenize_derivatives])
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
