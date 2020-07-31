
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
