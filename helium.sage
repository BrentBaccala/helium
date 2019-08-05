
# If I specify derivations=[x,y,z], then I have to use all three
# variables when constructing A(x,y,z), i.e, A(x,y) produces an error

(x,y,z) = var('x,y,z')
(A,B,r,Psi) = function('A,B,r,Psi')
(a,b,c,d,e) = var('a,b,c,d,e')
(f,g,h,i,j) = var('f,g,h,i,j')
var('E')

parameters = [a,b,c,d,e,  f,g,h,i,j,  E]

from sage.calculus.DifferentialAlgebra import DifferentialRing
DR = DifferentialRing(derivations = [x,y,z],
                      blocks = [[A,B,Psi,r], parameters],
	              parameters = parameters)

def H(Psi):
   return -diff(Psi,x,2)-diff(Psi,y,2)-diff(Psi,z,2)-(1/r(x,y,z))*Psi

rels = [A(x,y,z) == a*x + b*y + c*z + d*r(x,y,z) + e,
        B(x,y,z) == f*x + g*y + h*z + i*r(x,y,z) + j,
	r(x,y,z)^2 == x^2 + y^2 + z^2,
	H(Psi(x,y,z)) == E*Psi(x,y,z)]

RDC = DR.RosenfeldGroebner(rels)

print RDC

#exit()


var('x,y,z')

r = sqrt(x^2+y^2+z^2)

var('a,b,c,d,u')
A = a*x+b*y+c*z+u*r+d

var('e,f,g,h,v')
B = e*x+f*y+g*z+v*r+h

Psi = A*exp(B)
#Psi = a*exp(b*r)

var('E')

eq = -diff(Psi,x,2)-diff(Psi,y,2)-diff(Psi,z,2)-(1/r)*Psi-E*Psi

def bwb(expr):
    if isinstance(expr, Expression) and expr.operator():
       if expr.operator() == operator.pow and bool(expr.operands()[0] == (x^2+y^2+z^2)):
           return var('r')^(expr.operands()[1] * 2)
       else:
           return expr.operator()(*map(bwb, expr.operands()))
    else:
       return expr

# print numerator(expand(bwb(eq)/exp(B)))

BWB = QQ['a,b,c,d,u, e,f,g,h,v, E']
BWB2 = BWB['x,y,z,r']

#bwb2 = BWB2(numerator(expand(bwb(eq/exp(B)))))
bwb2 = BWB2(numerator(expand(bwb(eq/exp(B)))).subs({r^2:x^2+y^2+z^2, r^3:r*(x^2+y^2+z^2)}))

#print bwb2.coefficients()

# result = solve([SR(i) == 0 for i in bwb2.coefficients()], [E,a,b,c,d,u,e,f,g,h,v])
# print result

primes = ideal(bwb2.coefficients()).associated_primes()

for prime in primes:
    print prime._repr_short()
