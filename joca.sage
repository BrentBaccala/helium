# -*- mode: python -*-
#
# Sage script to do the computations for my paper in the Journal of Computational Algebra.
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
import DifferentialAlgebra

# Declare our independent variables
x,y,z = sympy.var('x,y,z')

# Declare our constants
E = sympy.var('E')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')
constants = [E,a0,a1,b0,b1,c0,c1,v1,v2,v3,v4]

# Declare our dependent variables
#
# `indexedbase` is DifferentialAlgebra's suggested way of declaring dependent variables
# if we want to use jet notation (i.e, v[x] for dv/dx) to write their derivatives.

Psi,DPsi,DDPsi = DifferentialAlgebra.indexedbase('Psi,DPsi,DDPsi')
v = DifferentialAlgebra.indexedbase('v')
r = DifferentialAlgebra.indexedbase('r')

# Create the DifferentialRing and set the ranking on the variables

R = DifferentialAlgebra.DifferentialRing (derivations = [x,y,z],
                                          blocks = [[DDPsi,DPsi,Psi],[v,r], constants],
                                          parameters = constants,
                                          notation='jet')

# Define the PDE we're trying to solve
# sympy can't handle a Sage Integer, so use a cast to make these Python integers

PDE = -(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - int(2) * Psi - int(2)*E*r*Psi

# Define the ansatz, the parameterized function space in which we're looking for solutions

syst = [Psi[x] - DPsi * v[x],
        Psi[y] - DPsi * v[y],
        Psi[z] - DPsi * v[z],
        DPsi[x] - DDPsi * v[x],
        DPsi[y] - DDPsi * v[y],
        DPsi[z] - DDPsi * v[z],
#       DifferentialAlgebra can't parse this parenthesized expression, so expand it out "by hand"
#       (a0 + a1*v)*DDPsi + (b0 + b1*v)*DPsi + (c0 + c1*v)*Psi,
        a0*DDPsi + a1*v*DDPsi + b0*DPsi + b1*v*DPsi + c0*Psi + c1*v*Psi,
        v - (v1*x + v2*y + v3*z + v4*r),
        r**2 - x**2 - y**2 - z**2]

# Reduce the PDE modulo the ansatz using Ritt's reduction algorithm

h,r = R.differential_prem(PDE, syst)

# Convert the remainder to Sage (since sympy can't compute prime decompositions of ideals)

Rsage = PolynomialRing(QQ, names=[str(indet) for indet in R.indets(r, selection='all')])
sage_constants = list(map(Rsage, constants))
sage_r = Rsage(r)

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

eqns = build_system_of_equations(sage_r, sage_constants)

# Build a polynomial ideal from the system of equations and construct its prime decomposition

I = ideal(eqns)
prime_decomposition = I.minimal_associated_primes()

# Print this result

print(*prime_decomposition, sep='\n')
