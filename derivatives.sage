
# Sage script to verify the substitutions in finish_prep() to express
# the second derivatives of Chi w.r.t. various variables in terms
# of Chi and DChi, where Chi is a second-order ODE w.r.t. B, i.e:
#
# DChi = diff(Chi,B)
# DDChi = diff(Chi,B,2)
#
# C * DDChi + D * DChi + F * Chi + G = 0
#
# (since E is already used for Energy)

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

var('v')
function('B C D F G Chi')

Bv = B(v)
Cv = C(v)
Dv = D(v)
Fv = F(v)
Gv = G(v)

ChiB = Chi(Bv)
ChiC = Chi(Cv)
ChiD = Chi(Dv)
ChiF = Chi(Fv)
ChiG = Chi(Gv)

# bwb = diff(ChiB,v,2).subs({DD[0](Chi)(Bv) : SR.var('DChi'), DD[0,0](Chi)(Bv): Dv/Cv*SR.var('DChi') + Fv/Cv*SR.var('Chi') + Gv/Cv})

bwb = diff(ChiB,v,2).subs({DD[0,0](Chi)(Bv) : Dv/Cv*SR.var('DChi') + Fv/Cv*SR.var('Chi') + Gv/Cv}).subs({DD[0](Chi)(Bv) : SR.var('DChi')})

print(bwb)

# compute diff(Chi(B(v)), v, 2) and substitute for...

# Cv * DD[0,0](Chi)(Bv) == Dv * DD[0](Chi)(Bv) + Fv * Chi(Bv) + Gv

# DD[0](Phi)(B) == Phi
# DD[0](Xi)(C) == 1/C
# DD[0,0](Zeta)(B) == F * DD[0](Zeta)(B) + G * Zeta
#
# B = trial_polynomial('b', coordinates, radii, 1)
# C = trial_polynomial('c', coordinates, radii, 1)
# F = trial_polynomial('f', B, [], 1)

# is DD[0](Phi)(B) the first derivative of Phi, with B as an argument?
# no, we really want the first derivative of Phi w.r.t B: DD[B](Phi) ?

# DD[B](Phi) == Phi
# DD[B](Xi) == 1/C
# DD[B,B](Zeta) == F * DD[B](Zeta) + G * Zeta

# diff(Phi(*coordinates),x1) == DD[0](Phi)(*coordinates)       (currently implemented; True)

# DD[B](Phi) == DD[x1](Phi) * diff(x1, B)

# DD[x1](Phi) == DD[B](Phi) * diff(B, x1)

# DD[x1](Phi(B(*coordinates))) == DD[B](Phi(B(*coordinates))) * diff(B, x1)

# Currently works:
#
# sage: diff(Phi(B(*coordinates)),x1)
# diff(B(x1, y1, z1), x1)*D[0](Phi)(B(x1, y1, z1))

# Doesn't work:
# sage: diff(Phi(B(*coordinates)),x1).subs({diff(Phi(B(*coordinates)), B(*coordinates)) : Phi(B(*coordinates))})
# TypeError: argument symb must be a symbol

# Works:
#
# sage: diff(Phi(B(*coordinates)),x1).subs({DD[0](Phi)(B(*coordinates)) : Phi(B(*coordinates))})
# Phi(B(x1, y1, z1))*diff(B(x1, y1, z1), x1)

# Works:
#
# sage: trial_B = trial_polynomial('b', coordinates, [], 1)[1]
# sage: trial_B
# b1*x1 + b2*y1 + b3*z1 + b0
# sage: diff(Phi(trial_B), x1)
# b1*D[0](Phi)(b1*x1 + b2*y1 + b3*z1 + b0)
# sage: diff(Phi(trial_B), x1).subs({DD[0](Phi)(trial_B) : Phi(trial_B)})
# b1*Phi(b1*x1 + b2*y1 + b3*z1 + b0)
