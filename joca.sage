# -*- mode: python -*-
#
# François Boulier's Differential Algebra package is required to run this script.
# To install it, run this command from the sage prompt:
#
# %pip install DifferentialAlgebra

import sympy
from DifferentialAlgebra import *

x,y,z = sympy.var('x,y,z')
E = sympy.var('E')
a0,a1,b0,b1,c0,c1 = sympy.var('a0,a1,b0,b1,c0,c1')
v1,v2,v3,v4 = sympy.var('v1,v2,v3,v4')

# `indexedbase` is DifferentialAlgebra's suggested way of declaring dependent variables
# if we want to use jet notation to write their derivatives (see DifferentialRing's docstring)

Psi,DPsi,DDPsi = indexedbase('Psi,DPsi,DDPsi')
v = indexedbase('v')
r = indexedbase('r')

params = [E,a0,a1,b0,b1,c0,c1,v1,v2,v3,v4]

R = DifferentialRing (derivations = [x,y,z], blocks = [[DDPsi,DPsi,Psi],[v,r], params], parameters = params, notation='jet')

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

# sympy can't handle a Sage Integer, so use a cast to make this Python intergers
PDE = -(Psi[x,x] + Psi[y,y] + Psi[z,z])*r - int(2) * Psi - int(2)*E*r*Psi

bwb = R.differential_prem(PDE, syst)

Rsage = PolynomialRing(QQ, names=[str(indet) for indet in R.indets(bwb[1].args[1], selection='all')])

def build_system_of_equations(eqn, constants):
    global system_of_like_terms
    ring = eqn.parent()
    system_of_like_terms = dict()
    global non_coeff_sub, non_coeff_part, monomial
    # cast int(1) is here because otherwise we get a sympy.core.numbers.One; don't ask me why
    non_coeff_sub = tuple(int(1) if ring.gen(n) in constants else ring.gen(n) for n in range(ring.ngens()))
    for i, (coeff, monomial) in enumerate(eqn):
        non_coeff_part = monomial(non_coeff_sub)
        # this cast needs to be here because otherwise the division (even though it's exact) takes us to the fraction field
        coeff_part = ring(monomial / non_coeff_part)
        if (non_coeff_part) in system_of_like_terms:
            system_of_like_terms[non_coeff_part] += coeff * coeff_part
        else:
            system_of_like_terms[non_coeff_part] = coeff * coeff_part
    return tuple(set(system_of_like_terms.values()))

eqn = Rsage(bwb[1].args[1])
eqns = build_system_of_equations(eqn, list(map(Rsage, params)))

I = ideal(eqns)
I.minimal_associated_primes()

print(*I.minimal_associated_primes(), sep='\n')
