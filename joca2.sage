# -*- mode: python -*-
#
# Sage script to do the computations for my paper in the Journal of Computational Algebra.
#
# Author: Brent Baccala
# Date: December 12, 2025
#
# Tested on Ubuntu 24 with Sage 10.6
#
# François Boulier's Differential Algebra package is required to run this script.
# To install it, run this command from the sage prompt:
#
# %pip install DifferentialAlgebra
#
# DifferentialAlgebra (sympy based) is used to perform the differential algebra reduction.
# Sage is used to compute a primary decomposition of the resulting system of equations.
# Some fiddling is required to juggle back and forth between the two.

import sympy
import time

try:
    import DifferentialAlgebra
except ModuleNotFoundError as ex:
    raise ModuleNotFoundError(ex.msg + "\nInstall it with '%pip install DifferentialAlgebra'")

# Claude Sonnet 4's solution to printing "DPsi" as "\Psi'" in LaTeX

def patch_latex_varify():
    """
    Patch sage.misc.latex.latex_varify to handle 'DPsi' specially.
    """
    from sage.misc.latex import latex_varify
    import sage.misc.latex

    # Save reference to the original function
    original_latex_varify = latex_varify

    # Define the custom function
    def custom_latex_varify(a, is_fname=False):
        if a == "DPsi":
            return r"\Psi'"
        else:
            return original_latex_varify(a, is_fname=is_fname)

    # Replace the function in the module
    sage.misc.latex.latex_varify = custom_latex_varify

patch_latex_varify()

# Wrap a long-running function in a timer and print the time

def timefunc(func, *args, **kwargs):
    start_time = time.perf_counter()
    retval = func(*args, **kwargs)
    end_time = time.perf_counter()
    print('{:30} {:10.2f} sec'.format(func.__name__, end_time - start_time))
    return retval

# Declare our independent variables
R1,R2,R12 = sympy.var('R1,R2,R12')

# Declare our constants
E = sympy.var('E')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')
g0,g1,g2,g3,g4,g5,g6,g7,g8,g9 = sympy.var('g0,g1,g2,g3,g4,g5,g6,g7,g8,g9')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E,v1,v2,v3,v4,a0,a1,b0,b1,c0,c1,g0,g1,g2,g3,g4,g5,g6,g7,g8,g9]

# Declare our dependent variables
#
# `indexedbase` is DifferentialAlgebra's suggested way of declaring dependent variables
# if we want to use jet notation (i.e, v[x] for dv/dx) to write their derivatives.

Psi,DPsi,DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')
gamma = DifferentialAlgebra.indexedbase('gamma')

# Create the DifferentialRing and set the ranking on the variables

# gamma has be lower than v, otherwise we get things like v[R1] in the reduction
# because              v - (v1*R1 + v2*R2 + v3*R12 + v4*gamma)
# differentiates to    v[R1] - v1 + v4 gamma[R1]
# if gamma ranked greater than v, then gamma[R1] would be the leader of this differential polynomial
# and it wouldn't be used to reduce v[R1]
#
# also                 gamma^2 - (g0 + g1*R1 + g2*R2 + g3*R12 + g4*R1^2 + g5*R2^2 + g6*R12^2 + g7*R1*R2 + g8*R1*R12 + g9*R2*R12)]
# differentiates to    2 gamma gamma[R1] - g1 - 2 g4 R1 - g7 R2 - g8 R12

DiffRing = DifferentialAlgebra.DifferentialRing (derivations = [R1,R2,R12],
                                                 blocks = [[DDPsi,DPsi,Psi,v,r,gamma], constants],
                                                 parameters = constants,
                                                 notation = 'jet')

# Define the PDE we're trying to solve
# sympy can't handle a Sage Integer, so use casts to make these Python integers

# eq (5) in Nakashima and Nakatusji, Solving the Schrodinger equation for helium...
# THE JOURNAL OF CHEMICAL PHYSICS 127, 224104 2007

PDE = - int(1)/int(2) *sum(Psi[Ri, Ri] + int(2)/Ri*Psi[Ri] for Ri in [R1,R2])           \
        - (Psi[R12, R12] + int(2)/R12*Psi[R12])                                         \
        - (R1^2 + R12^2 - R2^2)/(int(2)*R1*R12) * DiffRing.differentiate(Psi[R12],R1)   \
        - (R2^2 + R12^2 - R1^2)/(int(2)*R2*R12) * DiffRing.differentiate(Psi[R12],R2)   \
        - sum(int(2)/Ri for Ri in [R1,R2])*Psi + int(1)/R12*Psi                         \
        - E*Psi

print("PDE:", PDE)

# Define the ansatz, the parameterized function space in which we're looking for solutions

ansatz = [Psi[R1] - DPsi * v[R1],
          Psi[R2] - DPsi * v[R2],
          Psi[R12] - DPsi * v[R12],
          DPsi[R1] - DDPsi * v[R1],
          DPsi[R2] - DDPsi * v[R2],
          DPsi[R12] - DDPsi * v[R12],
          (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
          v - (v1*R1 + v2*R2 + v3*R12 + v4*gamma),
          gamma^2 - (g0 + g1*R1 + g2*R2 + g3*R12 + g4*R1^2 + g5*R2^2 + g6*R12^2 + g7*R1*R2 + g8*R1*R12 + g9*R2*R12)]

# DifferentialAlgebra can't handle the parenthesized expressions directly, so expand them

ansatz = list(map(sympy.expand, ansatz))

print("\nAnsatz:", *ansatz, sep='\n')

# Reduce the PDE modulo the ansatz using Ritt's reduction algorithm

numerator, denominator = sympy.expand(PDE).as_numer_denom()

h,r = timefunc(DiffRing.differential_prem, numerator, ansatz)

print(f"\nRemainder ({len(sympy.expand(r).args)} terms):", r)

# Convert the remainder to Sage

PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))
# len(r.args) is 4 and the first three don't involve constants
PolyRing_r = PolyRing(r.args[3])

# Given an equation and a list of constants, factor each term into constant and non-constant factors,
# then group together all terms with identical non-constant factors and return the resulting
# list of equations (which will only involve constants).

def build_system_of_equations(eqn, constants):
    ring = eqn.parent()
    system_of_like_terms = dict()
    non_constant_sub = tuple(1 if ring.gen(n) in constants else ring.gen(n) for n in range(ring.ngens()))
    for coeff, monomial in eqn:
        non_constant_part = monomial(non_constant_sub)
        constant_part = coeff * monomial // non_constant_part
        if (non_constant_part) in system_of_like_terms:
            system_of_like_terms[non_constant_part] += constant_part
        else:
            system_of_like_terms[non_constant_part] = constant_part
    return tuple(set(system_of_like_terms.values()))

eqns = timefunc(build_system_of_equations, PolyRing_r, PolyRing_constants)

print(f"\nSystem of equations ({len(eqns)} equations):", *eqns, sep='\n')

raise 'hi'

# Build a polynomial ideal from the system of equations and construct its prime decomposition

I = ideal(eqns)
prime_decomposition = I.minimal_associated_primes()

# Sort this result (so it prints in the same order as in the JOCA paper) and print it

prime_decomposition.sort(key=lambda x:str(x))
print("\nMinimal associated prime ideals:", *prime_decomposition, sep='\n')
