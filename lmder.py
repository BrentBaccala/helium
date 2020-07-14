
#       i is an integer input variable set to 1, 2, or 3 which
#         selects the desired machine parameter. If the machine has
#         t base b digits and its smallest and largest exponents are
#         emin and emax, respectively, then these parameters are
#
#         dpmpar[0] = b**(1 - t), the machine precision,
#
#         dpmpar[1] = b**(emin - 1), the smallest magnitude,
#
#         dpmpar[2] = b**emax*(1 - b**(-t)), the largest magnitude.
#
#       Original code used 1-based numbering


#     Machine constants for IEEE machines.
#

dpmpar = (2.22044604926e-16, 2.22507385852e-308, 1.79769313485e+308)

(one,p1,p5,p25,p75,p0001,zero) = (1.0e0, 1.0e-1, 5.0e-1, 2.5e-1, 7.5e-1, 1.0e-4, 0.0e0)

#     **********
#
#     subroutine lmder
#
#     the purpose of lmder is to minimize the sum of the squares of
#     m nonlinear functions in n variables by a modification of
#     the levenberg-marquardt algorithm. the user must provide a
#     subroutine which calculates the functions and the jacobian.
#
#     the subroutine statement is
#
#       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
#                        maxfev,diag,mode,factor,nprint,info,nfev,
#                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
#
#     where
#
#       fcn is the name of the user-supplied subroutine which
#         calculates the functions and the jacobian. fcn must
#         be declared in an external statement in the user
#         calling program, and should be written as follows.
#
#         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
#         integer m,n,ldfjac,iflag
#         double precision x(n),fvec(m),fjac(ldfjac,n)
#         ----------
#         if iflag = 1 calculate the functions at x and
#         return this vector in fvec. do not alter fjac.
#         if iflag = 2 calculate the jacobian at x and
#         return this matrix in fjac. do not alter fvec.
#         ----------
#         return
#         end
#
#         the value of iflag should not be changed by fcn unless
#         the user wants to terminate execution of lmder.
#         in this case set iflag to a negative integer.
#
#       m is a positive integer input variable set to the number
#         of functions.
#
#       n is a positive integer input variable set to the number
#         of variables. n must not exceed m.
#
#       x is an array of length n. on input x must contain
#         an initial estimate of the solution vector. on output x
#         contains the final estimate of the solution vector.
#
#       fvec is an output array of length m which contains
#         the functions evaluated at the output x.
#
#       fjac is an output m by n array. the upper n by n submatrix
#         of fjac contains an upper triangular matrix r with
#         diagonal elements of nonincreasing magnitude such that
#
#                t     t           t
#               p *(jac *jac)*p = r *r,
#
#         where p is a permutation matrix and jac is the final
#         calculated jacobian. column j of p is column ipvt(j)
#         (see below) of the identity matrix. the lower trapezoidal
#         part of fjac contains information generated during
#         the computation of r.
#
#       ldfjac is a positive integer input variable not less than m
#         which specifies the leading dimension of the array fjac.
#
#       ftol is a nonnegative input variable. termination
#         occurs when both the actual and predicted relative
#         reductions in the sum of squares are at most ftol.
#         therefore, ftol measures the relative error desired
#         in the sum of squares.
#
#       xtol is a nonnegative input variable. termination
#         occurs when the relative error between two consecutive
#         iterates is at most xtol. therefore, xtol measures the
#         relative error desired in the approximate solution.
#
#       gtol is a nonnegative input variable. termination
#         occurs when the cosine of the angle between fvec and
#         any column of the jacobian is at most gtol in absolute
#         value. therefore, gtol measures the orthogonality
#         desired between the function vector and the columns
#         of the jacobian.
#
#       maxfev is a positive integer input variable. termination
#         occurs when the number of calls to fcn with iflag = 1
#         has reached maxfev.
#
#       diag is an array of length n. if mode = 1 (see
#         below), diag is internally set. if mode = 2, diag
#         must contain positive entries that serve as
#         multiplicative scale factors for the variables.
#
#       mode is an integer input variable. if mode = 1, the
#         variables will be scaled internally. if mode = 2,
#         the scaling is specified by the input diag. other
#         values of mode are equivalent to mode = 1.
#
#       factor is a positive input variable used in determining the
#         initial step bound. this bound is set to the product of
#         factor and the euclidean norm of diag*x if nonzero, or else
#         to factor itself. in most cases factor should lie in the
#         interval (.1,100.).100. is a generally recommended value.
#
#       nprint is an integer input variable that enables controlled
#         printing of iterates if it is positive. in this case,
#         fcn is called with iflag = 0 at the beginning of the first
#         iteration and every nprint iterations thereafter and
#         immediately prior to return, with x, fvec, and fjac
#         available for printing. fvec and fjac should not be
#         altered. if nprint is not positive, no special calls
#         of fcn with iflag = 0 are made.
#
#       info is an integer output variable. if the user has
#         terminated execution, info is set to the (negative)
#         value of iflag. see description of fcn. otherwise,
#         info is set as follows.
#
#         info = 0  improper input parameters.
#
#         info = 1  both actual and predicted relative reductions
#                   in the sum of squares are at most ftol.
#
#         info = 2  relative error between two consecutive iterates
#                   is at most xtol.
#
#         info = 3  conditions for info = 1 and info = 2 both hold.
#
#         info = 4  the cosine of the angle between fvec and any
#                   column of the jacobian is at most gtol in
#                   absolute value.
#
#         info = 5  number of calls to fcn with iflag = 1 has
#                   reached maxfev.
#
#         info = 6  ftol is too small. no further reduction in
#                   the sum of squares is possible.
#
#         info = 7  xtol is too small. no further improvement in
#                   the approximate solution x is possible.
#
#         info = 8  gtol is too small. fvec is orthogonal to the
#                   columns of the jacobian to machine precision.
#
#       nfev is an integer output variable set to the number of
#         calls to fcn with iflag = 1.
#
#       njev is an integer output variable set to the number of
#         calls to fcn with iflag = 2.
#
#       ipvt is an integer output array of length n. ipvt
#         defines a permutation matrix p such that jac*p = q*r,
#         where jac is the final calculated jacobian, q is
#         orthogonal (not stored), and r is upper triangular
#         with diagonal elements of nonincreasing magnitude.
#         column j of p is column ipvt(j) of the identity matrix.
#
#       qtf is an output array of length n which contains
#         the first n elements of the vector (q transpose)*fvec.
#
#       wa1, wa2, and wa3 are work arrays of length n.
#
#       wa4 is a work array of length m.
#
#     subprograms called
#
#       user-supplied ...... fcn
#
#       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
#
#       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
#
#     argonne national laboratory. minpack project. march 1980.
#     burton s. garbow, kenneth e. hillstrom, jorge j. more
#
#     **********
def lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
          maxfev,diag,mode,factor,nprint,info,nfev,njev,
          ipvt,qtf,wa1,wa2,wa3,wa4):


#
#     epsmch is the machine precision.
#
    epsmch = dpmpar[0]
#
    info = 0
    iflag = 0
    nfev = 0
    njev = 0
#
#     check the input parameters for errors.
#
    if (n < 0 or m < n or ldfjac < m
           or ftol < zero or xtol < zero or gtol < zero
           or maxfev <= 0 or factor <= zero):
        if (iflag < 0): info = iflag
        iflag = 0
        if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
        return

    if (mode == 2):
        for j in range(1,n+1):
            if (diag[j+1] <= zero):
                if (iflag < 0): info = iflag
                iflag = 0
                if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
                return
#
#     evaluate the function at the starting point
#     and calculate its norm.
#
    iflag = 1
    fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    nfev = 1
    if (iflag < 0):
        if (iflag < 0): info = iflag
        iflag = 0
        if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
        return

    fnorm = enorm(m,fvec)
#
#     initialize levenberg-marquardt parameter and iteration counter.
#
    par = zero
    fiter = 1
#
#     beginning of the outer loop.
#
    while True:
   
#
#        calculate the jacobian matrix.
#
        iflag = 2
        fcn(m,n,x,fvec,fjac,ldfjac,iflag)
        njev = njev + 1
        if (iflag < 0):
            if (iflag < 0): info = iflag
            iflag = 0
            if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            return
#
#        if requested, call fcn to enable printing of iterates.
#
        if (nprint > 0):
            iflag = 0
            if (mod(fiter-1,nprint) == 0):
                fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            if (iflag < 0):
                if (iflag < 0): info = iflag
                iflag = 0
                if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
                return
#
#        compute the qr factorization of the jacobian.
#
        qrfac(m,n,fjac,ldfjac,True,ipvt,n,wa1,wa2,wa3)
#
#        on the first iteration and if mode is 1, scale according
#        to the norms of the columns of the initial jacobian.
#
        if (fiter == 1): # 80
            if (mode != 2):
                for j in range(1,n+1):
                    diag[j-1] = wa2[j-1]
                if (wa2[j-1] == zero): diag[j-1] = one
#
#        on the first iteration, calculate the norm of the scaled x
#        and initialize the step bound delta.
#
            for j in range(1,n+1):
                wa3[j-1] = diag[j-1]*x[j-1]
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            if (delta == zero): delta = factor
#
#        form (q transpose)*fvec and store the first n components in
#        qtf.
#
        for i in range(1,m+1):
            wa4[i-1] = fvec[i-1]

        for j in range(1,n+1):
            if (fjac[j-1,j-1] != zero):
                sum = zero
                for i in range(j, m+1):
                    sum = sum + fjac[i-1,j-1]*wa4[i-1]
                temp = -sum/fjac[j-1,j-1]
                for i in range(j, m+1):
                    wa4[i-1] = wa4[i-1] + fjac[i-1,j-1]*temp
            fjac[j-1,j-1] = wa1[j-1]
            qtf[j-1] = wa4[j-1]

#
#        compute the norm of the scaled gradient.
#
        gnorm = zero
        if (fnorm != zero):
            for j in range(1,n+1):
                l = ipvt[j-1]
                if (wa2[l-1] != zero):
                    sum = zero
                    for i in range(1,j+1):
                        sum = sum + fjac[i-1,j-1]*(qtf[i-1]/fnorm)
                    gnorm = dmax1(gnorm,dabs(sum/wa2(l)))

#
#        test for convergence of the gradient norm.
#
        if (gnorm <= gtol): info = 4
        if (info != 0):
            if (iflag < 0): info = iflag
            iflag = 0
            if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            return
#
#        rescale if necessary.
#
        if (mode != 2):
            for j in range(1,n+1):
                diag[j-1] = dmax1(diag[j-1],wa2[j-1])
#
#        beginning of the inner loop.
#
        while True:
#
#           determine the levenberg-marquardt parameter.
#
            lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4)
#
#           store the direction p and x + p. calculate the norm of p.
#
            for j in range(1,n+1):
                wa1[j-1] = -wa1[j-1]
                wa2[j-1] = x[j-1] + wa1[j-1]
                wa3[j-1] = diag[j-1]*wa1[j-1]
            pnorm = enorm(n,wa3)
#
#           on the first iteration, adjust the initial step bound.
#
            if (fiter == 1):
                delta = dmin1(delta,pnorm)
#
#           evaluate the function at x + p and calculate its norm.
#
            iflag = 1
            fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
            nfev = nfev + 1
            if (iflag < 0):
                if (iflag < 0): info = iflag
                iflag = 0
                if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
                return
            fnorm1 = enorm(m,wa4)
#
#           compute the scaled actual reduction.
#
            actred = -one
            if (p1*fnorm1 < fnorm): actred = one - (fnorm1/fnorm)**2
#
#           compute the scaled predicted reduction and
#           the scaled directional derivative.
#
            for j in range(1,n+1):
                wa3[j-1] = zero
                l = ipvt[j-1]
                temp = wa1[l-1]
                for i in range(1,j+1):
                    wa3[i-1] = wa3[i-1] + fjac[i-1,j-1]*temp

            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
#
#           compute the ratio of the actual to the predicted
#           reduction.
#
            ratio = zero
            if (prered != zero): ratio = actred/prered
#
#           update the step bound.
#
            if (ratio <= p25):
                if (actred >= zero): temp = p5
                if (actred < zero):
                    temp = p5*dirder/(dirder + p5*actred)
                if (p1*fnorm1 >= fnorm or temp < p1): temp = p1
                delta = temp*dmin1(delta,pnorm/p1)
                par = par/temp
            elif par == zero or ratio >= p75:
                delta = pnorm/p5
                par = p5*par
#
#           test for successful iteration.
#
            if (ratio >= p0001):
#
#           successful iteration. update x, fvec, and their norms.
#
                for j in range(1,n+1):
                    x[j-1] = wa2[j-1]
                    wa2[j-1] = diag[j-1]*x[j-1]
                for i in range(1,m+1):
                    fvec[i-1] = wa4[i-1]
                xnorm = enorm(n,wa2)
                fnorm = fnorm1
                fiter = fiter + 1

#
#           tests for convergence.
#
            if (dabs(actred) <= ftol and prered <= ftol and p5*ratio <= one): info = 1
            if (delta <= xtol*xnorm): info = 2
            if (dabs(actred) <= ftol and prered <= ftol and p5*ratio <= one and info == 2): info = 3
            if (info != 0):
                if (iflag < 0): info = iflag
                iflag = 0
                if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
                return
#
#           tests for termination and stringent tolerances.
#
            if (nfev >= maxfev): info = 5
            if (dabs(actred) <= epsmch and prered <= epsmch and p5*ratio <= one): info = 6
            if (delta <= epsmch*xnorm): info = 7
            if (gnorm <= epsmch): info = 8
            if (info != 0):
                if (iflag < 0): info = iflag
                iflag = 0
                if (nprint > 0): fcn(m,n,x,fvec,fjac,ldfjac,iflag)
                return
#
#           end of the inner loop. repeat if iteration unsuccessful.
#
            if (ratio >= p0001): break
#
#        end of the outer loop.
#
#
#     last card of subroutine lmder.
#

#     **********
#
#     subroutine lmpar
#
#     given an m by n matrix a, an n by n nonsingular diagonal
#     matrix d, an m-vector b, and a positive number delta,
#     the problem is to determine a value for the parameter
#     par such that if x solves the system
#
#           a*x = b ,     sqrt(par)*d*x = 0 ,
#
#     in the least squares sense, and dxnorm is the euclidean
#     norm of d*x, then either par is zero and
#
#           (dxnorm-delta) .le. 0.1*delta ,
#
#     or par is positive and
#
#           abs(dxnorm-delta) .le. 0.1*delta .
#
#     this subroutine completes the solution of the problem
#     if it is provided with the necessary information from the
#     qr factorization, with column pivoting, of a. that is, if
#     a*p = q*r, where p is a permutation matrix, q has orthogonal
#     columns, and r is an upper triangular matrix with diagonal
#     elements of nonincreasing magnitude, then lmpar expects
#     the full upper triangle of r, the permutation matrix p,
#     and the first n components of (q transpose)*b. on output
#     lmpar also provides an upper triangular matrix s such that
#
#            t   t                   t
#           p *(a *a + par*d*d)*p = s *s .
#
#     s is employed within lmpar and may be of separate interest.
#
#     only a few iterations are generally needed for convergence
#     of the algorithm. if, however, the limit of 10 iterations
#     is reached, then the output par will contain the best
#     value obtained so far.
#
#     the subroutine statement is
#
#       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
#                        wa1,wa2)
#
#     where
#
#       n is a positive integer input variable set to the order of r.
#
#       r is an n by n array. on input the full upper triangle
#         must contain the full upper triangle of the matrix r.
#         on output the full upper triangle is unaltered, and the
#         strict lower triangle contains the strict upper triangle
#         (transposed) of the upper triangular matrix s.
#
#       ldr is a positive integer input variable not less than n
#         which specifies the leading dimension of the array r.
#
#       ipvt is an integer input array of length n which defines the
#         permutation matrix p such that a*p = q*r. column j of p
#         is column ipvt(j) of the identity matrix.
#
#       diag is an input array of length n which must contain the
#         diagonal elements of the matrix d.
#
#       qtb is an input array of length n which must contain the first
#         n elements of the vector (q transpose)*b.
#
#       delta is a positive input variable which specifies an upper
#         bound on the euclidean norm of d*x.
#
#       par is a nonnegative variable. on input par contains an
#         initial estimate of the levenberg-marquardt parameter.
#         on output par contains the final estimate.
#
#       x is an output array of length n which contains the least
#         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
#         for the output par.
#
#       sdiag is an output array of length n which contains the
#         diagonal elements of the upper triangular matrix s.
#
#       wa1 and wa2 are work arrays of length n.
#
#     subprograms called
#
#       minpack-supplied ... dpmpar,enorm,qrsolv
#
#       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
#
#     argonne national laboratory. minpack project. march 1980.
#     burton s. garbow, kenneth e. hillstrom, jorge j. more
#
#     **********

def lmpar (n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2):

#
#     dwarf is the smallest positive magnitude.
#
    dwarf = dpmpar[1]
#
#     compute and store in x the gauss-newton direction. if the
#     jacobian is rank-deficient, obtain a least squares solution.
#
    nsing = n
    for j in range(1,n+1):
        wa1[j-1] = qtb[j-1]
        if (r[j-1,j-1] == 0 and nsing == n): nsing = j - 1
        if (nsing < n): wa1[j-1] = 0

    if (nsing >= 1):
        for k in range(1, nsing+1):
            j = nsing - k + 1
            wa1[j-1] = wa1[j-1]/r[j-1,j-1]
            temp = wa1[j-1]
            jm1 = j - 1
            if (jm1 >= 1):
                for i in range(1,jm1+1):
                    wa1[i-1] = wa1[i-1] - r[i-1,j-1]*temp

    for j in range(1,n+1):
        l = ipvt[j-1]
        x[l-1] = wa1[j-1]
#
#     initialize the iteration counter.
#     evaluate the function at the origin, and test
#     for acceptance of the gauss-newton direction.
#
    fiter = 0

    for j in range(1,n+1):
        wa2[j-1] = diag[j-1]*x[j-1]

    dxnorm = enorm(n,wa2)
    fp = dxnorm - delta
    if (fp <= p1*delta):
        if (fiter == 0): par = 0
        return
#
#     if the jacobian is not rank deficient, the newton
#     step provides a lower bound, parl, for the zero of
#     the function. otherwise set this bound to zero.
#
    parl = 0
    if (nsing >= n):
        for j in range(1,n+1):
            l = ipvt[j-1]
            wa1[j-1] = diag[l-1]*(wa2[l-1]/dxnorm)

        for j in range(1,n+1): # 110
            sum = 0
            jm1 = j - 1
            if (jm1 >= 1):
                for i in range(1,jm+1):
                    sum = sum + r[i-1,j-1]*wa1[i-1]
            wa1[j-1] = (wa1[j-1] - sum)/r[j-1,j-1]

        temp = enorm(n,wa1)
        parl = ((fp/delta)/temp)/temp
#
#     calculate an upper bound, paru, for the zero of the function.
#
    for j in range(1,n+1): # 140
        sum = zero
        for i in range(1,j+1):
            sum = sum + r[i-1,j-1]*qtb[i-1]
        l = ipvt[j-1]
        wa1[j-1] = sum/diag[l-1]
    gnorm = enorm(n,wa1)
    paru = gnorm/delta
    if (paru == 0): paru = dwarf/dmin1(delta,p1)
#
#     if the input par lies outside of the interval (parl,paru),
#     set par to the closer endpoint.
#
    par = dmax1(par,parl)
    par = dmin1(par,paru)
    if (par == 0): par = gnorm/dxnorm
#
#     beginning of an iteration.
#

    while True: # 150
        fiter = fiter + 1
#
#        evaluate the function at the current value of par.
#
        if (par == zero): par = dmax1(dwarf,p001*paru)
        temp = dsqrt(par)
        for j in range(1,n+1):
            wa1[j-1] = temp*diag[j-1]
        qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
        for j in range(1,n+1):
            wa2[j-1] = diag[j-1]*x[j-1]
        dxnorm = enorm(n,wa2)
        temp = fp
        fp = dxnorm - delta
#
#        if the function is small enough, accept the current value
#        of par. also test for the exceptional cases where parl
#        is zero or the number of iterations has reached 10.
#
        if (dabs(fp) < p1*delta or parl == zero and fp <= temp and temp < 0 or fiter == 10):
            if (fiter == 0): par = zero
            return
#
#        compute the newton correction.
#
        for j in range(1,n+1):
            l = ipvt[j-1]
            wa1[j-1] = diag[l-1]*(wa2[l-1]/dxnorm)

        for j in range(1,n+1):
            wa1[j-1] = wa1[j-1]/sdiag[j-1]
            temp = wa1[j-1]
            jp1 = j + 1
            if (n >= jp1):
                for i in range(jp1,n+1):
                    wa1[i-1] = wa1[i-1] - r[i-1,j-1]*temp

        temp = enorm(n,wa1)
        parc = ((fp/delta)/temp)/temp
#
#        depending on the sign of the function, update parl or paru.
#
        if (fp > 0): parl = dmax1(parl,par)
        if (fp < 0): paru = dmin1(paru,par)
#
#        compute an improved estimate for par.
#
        par = dmax1(parl,par+parc)
#
#        end of an iteration.
#
