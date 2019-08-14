

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

#rosenfeld_groebner()
#print
#print


var('x,y,z')

r = sqrt(x^2+y^2+z^2)

from itertools import *

def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)

def trial_polynomial(base, vars, degree):
    terms = tuple(flatten([combinations_with_replacement(vars, d) for d in range(degree+1)]))
    coefficients = tuple(var(base+str(c)) for c in range(len(terms)))
    poly = sum([var(base+str(c))*mul(v) for c,v in enumerate(terms)])
    return (coefficients, poly)

(Avars, A) = trial_polynomial('a', [x,y,z,r], 1)
(Bvars, B) = trial_polynomial('b', [x,y,z,r], 1)

Psi = A*exp(B)

var('E')

def H(Psi):
   return -diff(Psi,x,2)-diff(Psi,y,2)-diff(Psi,z,2)-(1/r)*Psi

eq = H(Psi) - E*Psi

# Now we want to replace all of the sqrt(...) factors with 'r',
# and we use a clever Python trick to build a dictionary
# that maps expressions to variable names.

import inspect

# need to use [3] instead of [2] because we're using this in a dict comprehension (I guess)
def varName(var):
    lcls = inspect.stack()[3][0].f_locals
    for name in lcls:
        if id(var) == id(lcls[name]):
            return name
    return None

def mk_maps(*vars):
    return {v.operands()[0] : var(varName(v)) for v in vars}

def bwb(expr):
    if isinstance(expr, Expression) and expr.operator():
       if expr.operator() == operator.pow and bool(expr.operands()[0] in maps):
           return maps[expr.operands()[0]]^(expr.operands()[1] * 2)
       else:
           return expr.operator()(*map(bwb, expr.operands()))
    else:
       return expr

maps = mk_maps(r)

# print maps

# print numerator(expand(bwb(eq)/exp(B)))

# use custom term ordering to prioritize elimination of 'E' variable
# if Sage's Groebner basis-based techniques are used
order = TermOrder('deglex(1),degrevlex({})'.format(len(Avars)+len(Bvars)))
BWB = PolynomialRing(QQ, (var('E'),) + Avars + Bvars, order=order)
#BWB = PolynomialRing(QQ, (var('E'),) + Avars + Bvars)

BWB2.<x,y,z,r> = Frac(BWB)[]
BWB3 = BWB2.quo(r^2-(x^2+y^2+z^2))

bwb3 = BWB3(numerator(expand(bwb(eq/exp(B)))))
eqns = map(numerator, bwb3.lift().coefficients())
bwbI = ideal(eqns)

for poly in eqns:
    print poly

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

def associated_primes():

    primes = bwbI.associated_primes()

    for prime in primes:
        print prime._repr_short()


# Look for solutions using an approximate numerical technique

def random_numerical():

    import scipy.optimize

    # eqns.append(E+1/4)
    # eqns.append(BWB.gen(1) - 1)
    minpoly = sum([poly*poly for poly in eqns])

    nvars = len(BWB.gens())

    iv = [random() for i in range(nvars)]

    def square(x):
        return x*x

    def minfunc(v):
        return sqrt(minpoly.subs(dict(zip(BWB.gens(), v))) / sum(map(square, v)))

    SciMin = scipy.optimize.minimize(minfunc, iv, method='Nelder-Mead')

    #print SciMin

    print
    print

    if SciMin.success:
        for pair in zip(BWB.gens(), SciMin.x): print pair
    else:
        print SciMin.message


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
