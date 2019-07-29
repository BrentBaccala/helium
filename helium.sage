
var('x,y,z')

r = sqrt(x^2+y^2+z^2)

var('a,b,c,d,u')
A = a*x+b*y+c*z+u*r+d

var('e,f,g,h,v')
B = e*x+f*y+g*z+v*r+h

#Psi = A*exp(B)
Psi = a*exp(b*r)

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
bwb2 = BWB2(numerator(expand(bwb(eq/Psi))).subs({r^2:x^2+y^2+z^2, r^3:r*(x^2+y^2+z^2)}))

# print bwb2.coefficients()

result = solve([SR(i) == 0 for i in bwb2.coefficients()], [E,b])

print result
