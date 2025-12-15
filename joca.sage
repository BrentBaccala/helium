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

# Declare our independent variables
x,y,z = sympy.var('x,y,z')

# Declare our constants
E = sympy.var('E')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
constants = [E,v1,v2,v3,v4,a0,a1,b0,b1,c0,c1]

# Declare our dependent variables
#
# `indexedbase` is DifferentialAlgebra's suggested way of declaring dependent variables
# if we want to use jet notation (i.e, v[x] for dv/dx) to write their derivatives.

Psi,DPsi,DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

# Create the DifferentialRing and set the ranking on the variables

DiffRing = DifferentialAlgebra.DifferentialRing (derivations = [x,y,z],
                                                 blocks = [[DDPsi,DPsi,Psi,v,r], constants],
                                                 parameters = constants,
                                                 notation = 'jet')

# Define the PDE we're trying to solve
# sympy can't handle a Sage Integer, so use casts to make these Python integers

PDE = -int(1)/int(2)*(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - Psi - E*r*Psi

print("PDE:", PDE)

# Define the ansatz, the parameterized function space in which we're looking for solutions

ansatz = [Psi[x] - DPsi * v[x],
          Psi[y] - DPsi * v[y],
          Psi[z] - DPsi * v[z],
          DPsi[x] - DDPsi * v[x],
          DPsi[y] - DDPsi * v[y],
          DPsi[z] - DDPsi * v[z],
          (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
          v - (v1*x + v2*y + v3*z + v4*r),
          r**2 - x**2 - y**2 - z**2]

# DifferentialAlgebra can't handle the parenthesized expressions directly, so expand them

ansatz = list(map(sympy.expand, ansatz))

print("\nAnsatz:", *ansatz, sep='\n')

# Reduce the PDE modulo the ansatz using Ritt's reduction algorithm

h,r = DiffRing.differential_prem(PDE, ansatz)

print("\nRemainder:", r)

# Convert the remainder to Sage

PolyRing = PolynomialRing(QQ, names=[str(indet) for indet in DiffRing.indets(selection='all')])
PolyRing_constants = list(map(PolyRing, constants))
PolyRing_r = PolyRing(r)

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

eqns = build_system_of_equations(PolyRing_r, PolyRing_constants)

print("\nSystem of equations:", *eqns, sep='\n')

# Build a polynomial ideal from the system of equations and construct its prime decomposition

I = ideal(eqns)
prime_decomposition = I.minimal_associated_primes()

# Sort this result (so it prints in the same order as in the JOCA paper) and print it

prime_decomposition.sort(key=lambda x:str(x))
print("\nMinimal associated prime ideals:", *prime_decomposition, sep='\n')
