c Delivery-Date:  8 November 1991 16:31 mst
c Delivery-By:  Network_Server.Daemon (LISTSERV@HEARN.BITNET@UNCANVE)
c Date:  Friday, 8 November 1991 10:19 mst
c From:  LISTSERV at HEARN
c To:  BOGERT at UNCAMULT
c
c GCVSPL.FOR, 1986-05-12
c
c***********************************************************************
c
c SUBROUTINE GCVSPL (REAL*8)
c
c Purpose:
c *******
c
c       Natural B-spline data smoothing subroutine, using the Generali-
c       zed Cross-Validation and Mean-Squared Prediction Error Criteria
c       of Craven & Wahba (1979). Alternatively, the amount of smoothing
c       can be given explicitly, or it can be based on the effective
c       number of degrees of freedom in the smoothing process as defined
c       by Wahba (1980). The model assumes uncorrelated, additive noise
c       and essentially smooth, underlying functions. The noise may be
c       non-stationary, and the independent co-ordinates may be spaced
c       non-equidistantly. Multiple datasets, with common independent
c       variables and weight factors are accomodated.
c
c
c Calling convention:
c ******************
c
c       CALL GCVSPL ( X, Y, NY, WX, WY, M, N, K, MD, VAL, C, NC, WK, IER
c )
c
c Meaning of parameters:
c *********************
c
c       X(N)    ( I )   Independent variables: strictly increasing knot
c                       sequence, with X(I-1).lt.X(I), I=2,...,N.
c       Y(NY,K) ( I )   Input data to be smoothed (or interpolated).
c       NY      ( I )   First dimension of array Y(NY,K), with NY.ge.N.
c       WX(N)   ( I )   Weight factor array; WX(I) corresponds with
c                       the relative inverse variance of point Y(I,*).
c                       If no relative weighting information is
c                       available, the WX(I) should be set to ONE.
c                       All WX(I).gt.ZERO, I=1,...,N.
c       WY(K)   ( I )   Weight factor array; WY(J) corresponds with
c                       the relative inverse variance of point Y(*,J).
c                       If no relative weighting information is
c                       available, the WY(J) should be set to ONE.
c                       All WY(J).gt.ZERO, J=1,...,K.
c                       NB: The effective weight for point Y(I,J) is
c                       equal to WX(I)*WY(J).
c       M       ( I )   Half order of the required B-splines (spline
c                       degree 2*M-1), with M.gt.0. The values M =
c                       1,2,3,4 correspond to linear, cubic, quintic,
c                       and heptic splines, respectively.
c       N       ( I )   Number of observations per dataset, with N.ge.2*
cM.
c       K       ( I )   Number of datasets, with K.ge.1.
c       MD      ( I )   Optimization mode switch:
c                       |MD| = 1: Prior given value for p in VAL
c                                 (VAL.ge.ZERO). This is the fastest
c                                 use of GCVSPL, since no iteration
c                                 is performed in p.
c                       |MD| = 2: Generalized cross validation.
c                       |MD| = 3: True predicted mean-squared error,
c                                 with prior given variance in VAL.
c                       |MD| = 4: Prior given number of degrees of
c                                 freedom in VAL (ZERO.le.VAL.le.N-M).
c                        MD  < 0: It is assumed that the contents of
c                                 X, W, M, N, and WK have not been
c                                 modified since the previous invoca-
c                                 tion of GCVSPL. If MD < -1, WK(4)
c                                 is used as an initial estimate for
c                                 the smoothing parameter p.  At the
c                                 first call to GCVSPL, MD must be > 0.
c                       Other values for |MD|, and inappropriate values
c                       for VAL will result in an error condition, or
c                       cause a default value for VAL to be selected.
c                       After return from MD.ne.1, the same number of
c                       degrees of freedom can be obtained, for identica
cl
c                       weight factors and knot positions, by selecting
c                       |MD|=1, and by copying the value of p from WK(4)
c                       into VAL. In this way, no iterative optimization
c                       is required when processing other data in Y.
c       VAL     ( I )   Mode value, as described above under MD.
c       C(NC,K) ( O )   Spline coefficients, to be used in conjunction
c                       with function SPLDER. NB: the dimensions of C
c                       in GCVSPL and in SPLDER are different! In SPLDER
c,
c                       only a single column of C(N,K) is needed, and th
ce
c                       proper column C(1,J), with J=1...K should be use
cd
c                       when calling SPLDER.
c       NC       ( I )  First dimension of array C(NC,K), NC.ge.N.
c       WK(IWK) (I/W/O) Work vector, with length IWK.ge.6*(N*M+1)+N.
c                       On normal exit, the first 6 values of WK are
c                       assigned as follows:
c
c                       WK(1) = Generalized Cross Validation value
c                       WK(2) = Mean Squared Residual.
c                       WK(3) = Estimate of the number of degrees of
c                               freedom of the residual sum of squares
c                               per dataset, with 0.lt.WK(3).lt.N-M.
c                       WK(4) = Smoothing parameter p, multiplicative
c                               with the splines' derivative constraint.
c                       WK(5) = Estimate of the true mean squared error
c                               (different formula for |MD| = 3).
c                       WK(6) = Gauss-Markov error variance.
c
c                       If WK(4) -->  0 , WK(3) -->  0 , and an inter-
c                       polating spline is fitted to the data (p --> 0).
c                       A very small value > 0 is used for p, in order
c                       to avoid division by zero in the GCV function.
c
c                       If WK(4) --> inf, WK(3) --> N-M, and a least-
c                       squares polynomial of order M (degree M-1) is
c                       fitted to the data (p --> inf). For numerical
c                       reasons, a very high value is used for p.
c
c                       Upon return, the contents of WK can be used for
c                       covariance propagation in terms of the matrices
c                       B and WE: see the source listings. The variance
c                       estimate for dataset J follows as WK(6)/WY(J).
c
c       IER     ( O )   Error parameter:
c
c                       IER = 0:        Normal exit
c                       IER = 1:        M.le.0 .or. N.lt.2*M
c                       IER = 2:        Knot sequence is not strictly
c                                       increasing, or some weight
c                                       factor is not positive.
c                       IER = 3:        Wrong mode  parameter or value.
c
c Remarks:
c *******
c
c       (1) GCVSPL calculates a natural spline of order 2*M (degree
c       2*M-1) which smoothes or interpolates a given set of data
c       points, using statistical considerations to determine the
c       amount of smoothing required (Craven & Wahba, 1979). If the
c       error variance is a priori known, it should be supplied to
c       the routine in VAL, for |MD|=3. The degree of smoothing is
c       then determined to minimize an unbiased estimate of the true
c       mean squared error. On the other hand, if the error variance
c       is not known, one may select |MD|=2. The routine then deter-
c       mines the degree of smoothing to minimize the generalized
c       cross validation function. This is asymptotically the same
c       as minimizing the true predicted mean squared error (Craven &
c       Wahba, 1979). If the estimates from |MD|=2 or 3 do not appear
c       suitable to the user (as apparent from the smoothness of the
c       M-th derivative or from the effective number of degrees of
c       freedom returned in WK(3) ), the user may select an other
c       value for the noise variance if |MD|=3, or a reasonably large
c       number of degrees of freedom if |MD|=4. If |MD|=1, the proce-
c       dure is non-iterative, and returns a spline for the given
c       value of the smoothing parameter p as entered in VAL.
c
c       (2) The number of arithmetic operations and the amount of
c       storage required are both proportional to N, so very large
c       datasets may be accomodated. The data points do not have
c       to be equidistant in the independant variable X or uniformly
c       weighted in the dependant variable Y. However, the data
c       points in X must be strictly increasing. Multiple dataset
c       processing (K.gt.1) is numerically more efficient dan
c       separate processing of the individual datasets (K.eq.1).
c
c       (3) If |MD|=3 (a priori known noise variance), any value of
c       N.ge.2*M is acceptable. However, it is advisable for N-2*M
c       be rather large (at least 20) if |MD|=2 (GCV).
c
c       (4) For |MD| > 1, GCVSPL tries to iteratively minimize the
c       selected criterion function. This minimum is unique for |MD|
c       = 4, but not necessarily for |MD| = 2 or 3. Consequently,
c       local optima rather that the global optimum might be found,
c       and some actual findings suggest that local optima might
c       yield more meaningful results than the global optimum if N
c       is small. Therefore, the user has some control over the
c       search procedure. If MD > 1, the iterative search starts
c       from a value which yields a number of degrees of freedom
c       which is approximately equal to N/2, until the first (local)
c       minimum is found via a golden section search procedure
c       (Utreras, 1980). If MD < -1, the value for p contained in
c       WK(4) is used instead. Thus, if MD = 2 or 3 yield too noisy
c       an estimate, the user might try |MD| = 1 or 4, for suitably
c       selected values for p or for the number of degrees of
c       freedom, and then run GCVSPL with MD = -2 or -3. The con-
c       tents of N, M, K, X, WX, WY, and WK are assumed unchanged
c       since the last call to GCVSPL if MD < 0.
c
c       (5) GCVSPL calculates the spline coefficient array C(N,K);
c       this array can be used to calculate the spline function
c       value and any of its derivatives up to the degree 2*M-1
c       at any argument T within the knot range, using subrou-
c       tines SPLDER and SEARCH, and the knot array X(N). Since
c       the splines are constrained at their Mth derivative, only
c       the lower spline derivatives will tend to be reliable
c       estimates of the underlying, true signal derivatives.
c
c       (6) GCVSPL combines elements of subroutine CRVO5 by Utre-
c       ras (1980), subroutine SMOOTH by Lyche et al. (1983), and
c       subroutine CUBGCV by Hutchinson (1985). The trace of the
c       influence matrix is assessed in a similar way as described
c       by Hutchinson & de Hoog (1985). The major difference is
c       that the present approach utilizes non-symmetrical B-spline
c       design matrices as described by Lyche et al. (1983); there-
c       fore, the original algorithm by Erisman & Tinney (1975) has
c       been used, rather than the symmetrical version adopted by
c       Hutchinson & de Hoog.
c
c References:
c **********
c
c       P. Craven & G. Wahba (1979), Smoothing noisy data with
c       spline functions. Numerische Mathematik 31, 377-403.
c
c       A.M. Erisman & W.F. Tinney (1975), On computing certain
c       elements of the inverse of a sparse matrix. Communications
c       of the ACM 18(3), 177-179.
c
c       M.F. Hutchinson & F.R. de Hoog (1985), Smoothing noisy data
c       with spline functions. Numerische Mathematik 47(1), 99-106.
c
c       M.F. Hutchinson (1985), Subroutine CUBGCV. CSIRO Division of
c       Mathematics and Statistics, P.O. Box 1965, Canberra, ACT 2601,
c       Australia.
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori (1983), Fortran
c       subroutines for computing smoothing and interpolating natural
c       splines. Advances in Engineering Software 5(1), 2-5.
c
c       F. Utreras (1980), Un paquete de programas para ajustar curvas
c       mediante funciones spline. Informe Tecnico MA-80-B-209, Depar-
c       tamento de Matematicas, Faculdad de Ciencias Fisicas y Matema-
c       ticas, Universidad de Chile, Santiago.
c
c       Wahba, G. (1980). Numerical and statistical methods for mildly,
c       moderately and severely ill-posed problems with noisy data.
c       Technical report nr. 595 (February 1980). Department of Statis-
c       tics, University of Madison (WI), U.S.A.
c
c Subprograms required:
c ********************
c
c       BASIS, PREP, SPLC, BANDET, BANSOL, TRINV
c
c***********************************************************************
c
c
      subroutine gcvspl(x, y, ny, wx, wy, m, n, k, md, val, c, nc, wk, 
     &ier)
cf2py intent(out) :: C
cf2py intent(out) :: WK
cf2py intent(out) :: IER
      implicit double precision (o-z, a-h)
      parameter (ratio = 2d0, tau = 1.618033983d0, ibwe = 7, zero = 0d0
     &, half = 5d-1, one = 1d0, tol = 1d-6, eps = 1d-15, epsinv = one / 
     &eps)
      dimension x(n), y(ny, k), wx(n), wy(k), c(nc, k), wk(n + (6 * ((n
     & * m) + 1)))
      save el, nm1, m2
c
c***  Parameter check and work array initialization
c
c***  Check on mode parameter
      data m2 / 0 /
      data nm1 / 0 /
      data el / 0d0 /
      ier = 0
      if (((((iabs(md) .gt. 4) .or. (md .eq. 0)) .or. ((iabs(md) .eq. 1)
     & .and. (val .lt. zero))) .or. ((iabs(md) .eq. 3) .and. (val .lt. 
     &zero))) .or. ((iabs(md) .eq. 4) .and. ((val .lt. zero) .or. (val
     & .gt. (n - m))))) then
cWrong mode value                                 
      ier = 3
      return 
c***  Check on M and N
      end if
      if (md .gt. 0) then
      m2 = 2 * m
      nm1 = n - 1
      else
      if ((m2 .ne. (2 * m)) .or. (nm1 .ne. (n - 1))) then
cM or N modified since previous call           
      ier = 3
      return 
      end if
      end if
      if ((m .le. 0) .or. (n .lt. m2)) then
cM or N invalid                                   
      ier = 1
      return 
c***  Check on knot sequence and weights
      end if
      if (wx(1) .le. zero) ier = 2
      do 10 i = 2, n
      if ((wx(i) .le. zero) .or. (x(i - 1) .ge. x(i))) ier = 2
      if (ier .ne. 0) return 
   10 continue
      do 15 j = 1, k
      if (wy(j) .le. zero) ier = 2
      if (ier .ne. 0) return 
c
c***  Work array parameters (address information for covariance
c***  propagation by means of the matrices STAT, B, and WE). NB:
c***  BWE cannot be used since it is modified by function TRINV.
c
   15 continue
      nm2p1 = n * (m2 + 1)
c     ISTAT = 1            !Statistics array STAT(6)
c     IBWE  = ISTAT + 6      !Smoothing matrix BWE( -M:M  ,N)
      nm2m1 = n * (m2 - 1)
cDesign matrix    B  (1-M:M-1,N)       
      ib = ibwe + nm2p1
c     IWK   = IWE   + NM2P1      !Total work array length N + 6*(N*M+1)
c
c***  Compute the design matrices B and WE, the ratio
c***  of their L1-norms, and check for iterative mode.
c
cDesign matrix    WE ( -M:M  ,N)       
      iwe = ib + nm2m1
      if (md .gt. 0) then
      call basis(m, n, x, wk(ib), r1, wk(ibwe))
      call prep(m, n, x, wx, wk(iwe), el)
cL1-norms ratio (SAVEd upon RETURN)          
      el = el / r1
      end if
c***     Prior given value for p
      if (iabs(md) .ne. 1) goto 20
      r1 = val
c
c***  Iterate to minimize the GCV function (|MD|=2),
c***  the MSE function (|MD|=3), or to obtain the prior
c***  given number of degrees of freedom (|MD|=4).
c
      goto 100
   20 if (md .lt. (-1)) then
cUser-determined starting value                
      r1 = wk(4)
      else
cDefault (DOF ~ 0.5)                        
      r1 = one / el
      end if
      r2 = r1 * ratio
      gf2 = splc(m,n,k,y,ny,wx,wy,md,val,r2,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
   40 gf1 = splc(m,n,k,y,ny,wx,wy,md,val,r1,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      if (gf1 .gt. gf2) goto 50
cInterpolation         
      if (wk(4) .le. zero) goto 100
      r2 = r1
      gf2 = gf1
      r1 = r1 / ratio
      goto 40
   50 r3 = r2 * ratio
   60 gf3 = splc(m,n,k,y,ny,wx,wy,md,val,r3,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      if (gf3 .gt. gf2) goto 70
cLeast-squares polynomial  
      if (wk(4) .ge. epsinv) goto 100
      r2 = r3
      gf2 = gf3
      r3 = r3 * ratio
      goto 60
   70 r2 = r3
      gf2 = gf3
      alpha = (r2 - r1) / tau
      r4 = r1 + alpha
      r3 = r2 - alpha
      gf3 = splc(m,n,k,y,ny,wx,wy,md,val,r3,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      gf4 = splc(m,n,k,y,ny,wx,wy,md,val,r4,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
   80 if (gf3 .le. gf4) then
      r2 = r4
      gf2 = gf4
      err = (r2 - r1) / (r1 + r2)
      if ((((err * err) + one) .eq. one) .or. (err .le. tol)) goto 90
      r4 = r3
      gf4 = gf3
      alpha = alpha / tau
      r3 = r2 - alpha
      gf3 = splc(m,n,k,y,ny,wx,wy,md,val,r3,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      else
      r1 = r3
      gf1 = gf3
      err = (r2 - r1) / (r1 + r2)
      if ((((err * err) + one) .eq. one) .or. (err .le. tol)) goto 90
      r3 = r4
      gf3 = gf4
      alpha = alpha / tau
      r4 = r1 + alpha
      gf4 = splc(m,n,k,y,ny,wx,wy,md,val,r4,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      end if
      goto 80
c
c***  Calculate final spline coefficients
c
   90 r1 = half * (r1 + r2)
c
c***  Ready
c
  100 gf1 = splc(m,n,k,y,ny,wx,wy,md,val,r1,eps,c,nc,wk,wk(ib),wk(iwe),
     &el,wk(ibwe))
      return 
c BASIS.FOR, 1985-06-03
c
c***********************************************************************
c
c SUBROUTINE BASIS (REAL*8)
c
c Purpose:
c *******
c
c       Subroutine to assess a B-spline tableau, stored in vectorized
c       form.
c
c Calling convention:
c ******************
c
c       CALL BASIS ( M, N, X, B, BL, Q )
c
c Meaning of parameters:
c *********************
c
c       M               ( I )   Half order of the spline (degree 2*M-1),
c                               M > 0.
c       N               ( I )   Number of knots, N >= 2*M.
c       X(N)            ( I )   Knot sequence, X(I-1) < X(I), I=2,N.
c       B(1-M:M-1,N)    ( O )   Output tableau. Element B(J,I) of array
c                               B corresponds with element b(i,i+j) of
c                               the tableau matrix B.
c       BL              ( O )   L1-norm of B.
c       Q(1-M:M)        ( W )   Internal work array.
c
c Remark:
c ******
c
c       This subroutine is an adaptation of subroutine BASIS from the
c       paper by Lyche et al. (1983). No checking is performed on the
c       validity of M and N. If the knot sequence is not strictly in-
c       creasing, division by zero may occur.
c
c Reference:
c *********
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      subroutine basis(m, n, x, b, bl, q)
      implicit double precision (o-z, a-h)
      parameter (zero = 0d0, one = 1d0)
c
      dimension x(n), b(1 - m:m - 1, n), q(1 - m:m)
c***         Linear spline
      if (m .eq. 1) then
      do 3 i = 1, n
      b(0,i) = one
    3 continue
      bl = one
      return 
c
c***  General splines
c
      end if
      mm1 = m - 1
      mp1 = m + 1
      m2 = 2 * m
c***     1st row
      do 15 l = 1, n
      do 5 j = - mm1, m
      q(j) = zero
    5 continue
      q(mm1) = one
c***     Successive rows
      if ((l .ne. 1) .and. (l .ne. n)) q(mm1) = one / (x(l + 1) - x(l - 
     &1))
      arg = x(l)
      do 13 i = 3, m2
      ir = mp1 - i
      v = q(ir)
c***               Left-hand B-splines
      if (l .lt. i) then
      do 6 j = l + 1, i
      u = v
      v = q(ir + 1)
      q(ir) = u + ((x(j) - arg) * v)
      ir = ir + 1
    6 continue
      end if
      j1 = max0((l - i) + 1,1)
      j2 = min0(l - 1,n - i)
c***               Ordinary B-splines
      if (j1 .le. j2) then
      if (i .lt. m2) then
      do 8 j = j1, j2
      y = x(i + j)
      u = v
      v = q(ir + 1)
      q(ir) = u + (((v - u) * (y - arg)) / (y - x(j)))
      ir = ir + 1
    8 continue
      else
      do 10 j = j1, j2
      u = v
      v = q(ir + 1)
      q(ir) = ((arg - x(j)) * u) + ((x(i + j) - arg) * v)
      ir = ir + 1
   10 continue
      end if
      end if
      nmip1 = (n - i) + 1
c***           Right-hand B-splines
      if (nmip1 .lt. l) then
      do 12 j = nmip1, l - 1
      u = v
      v = q(ir + 1)
      q(ir) = ((arg - x(j)) * u) + v
      ir = ir + 1
   12 continue
      end if
   13 continue
      do 14 j = - mm1, mm1
      b(j,l) = q(j)
   14 continue
c
c***  Zero unused parts of B
c
   15 continue
      do 17 i = 1, mm1
      do 16 k = i, mm1
      b(- k,i) = zero
      b(k,(n + 1) - i) = zero
   16 continue
c
c***  Assess L1-norm of B
c
   17 continue
      bl = 0d0
      do 19 i = 1, n
      do 18 k = - mm1, mm1
      bl = bl + abs(b(k,i))
   18 continue
   19 continue
c
c***  Ready
c
      bl = bl / n
      return 
c PREP.FOR, 1985-07-04
c
c***********************************************************************
c
c SUBROUTINE PREP (REAL*8)
c
c Purpose:
c *******
c
c       To compute the matrix WE of weighted divided difference coeffi-
c       cients needed to set up a linear system of equations for sol-
c       ving B-spline smoothing problems, and its L1-norm EL. The matrix
c       WE is stored in vectorized form.
c
c Calling convention:
c ******************
c
c       CALL PREP ( M, N, X, W, WE, EL )
c
c Meaning of parameters:
c *********************
c
c       M               ( I )   Half order of the B-spline (degree
c                               2*M-1), with M > 0.
c       N               ( I )   Number of knots, with N >= 2*M.
c       X(N)            ( I )   Strictly increasing knot array, with
c                               X(I-1) < X(I), I=2,N.
c       W(N)            ( I )   Weight matrix (diagonal), with
c                               W(I).gt.0.0, I=1,N.
c       WE(-M:M,N)      ( O )   Array containing the weighted divided
c                               difference terms in vectorized format.
c                               Element WE(J,I) of array E corresponds
c                               with element e(i,i+j) of the matrix
c                               W**-1 * E.
c       EL              ( O )   L1-norm of WE.
c
c Remark:
c ******
c
c       This subroutine is an adaptation of subroutine PREP from the pap
cer
c       by Lyche et al. (1983). No checking is performed on the validity
c       of M and N. Division by zero may occur if the knot sequence is
c       not strictly increasing.
c
c Reference:
c *********
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      subroutine prep(m, n, x, w, we, el)
      implicit double precision (o-z, a-h)
      parameter (zero = 0d0, one = 1d0)
c
c***  Calculate the factor F1
c
cWE(-M:M,N)              
      dimension x(n), w(n), we(((2 * m) + 1) * n)
      m2 = 2 * m
      mp1 = m + 1
      m2m1 = m2 - 1
      m2p1 = m2 + 1
      nm = n - m
      f1 = - one
      if (m .ne. 1) then
      do 5 i = 2, m
      f1 = - (f1 * i)
    5 continue
      do 6 i = mp1, m2m1
      f1 = f1 * i
    6 continue
c
c***  Columnwise evaluation of the unweighted design matrix E
c
      end if
      i1 = 1
      i2 = m
      jm = mp1
      do 17 j = 1, n
      inc = m2p1
      if (j .gt. nm) then
      f1 = - f1
      f = f1
      else
      if (j .lt. mp1) then
      inc = 1
      f = f1
      else
      f = f1 * (x(j + m) - x(j - m))
      end if
      end if
      if (j .gt. mp1) i1 = i1 + 1
      if (i2 .lt. n) i2 = i2 + 1
c***     Loop for divided difference coefficients
      jj = jm
      ff = f
      y = x(i1)
      i1p1 = i1 + 1
      do 11 i = i1p1, i2
      ff = ff / (y - x(i))
   11 continue
      we(jj) = ff
      jj = jj + m2
      i2m1 = i2 - 1
      if (i1p1 .le. i2m1) then
      do 14 l = i1p1, i2m1
      ff = f
      y = x(l)
      do 12 i = i1, l - 1
      ff = ff / (y - x(i))
   12 continue
      do 13 i = l + 1, i2
      ff = ff / (y - x(i))
   13 continue
      we(jj) = ff
      jj = jj + m2
   14 continue
      end if
      ff = f
      y = x(i2)
      do 16 i = i1, i2m1
      ff = ff / (y - x(i))
   16 continue
      we(jj) = ff
      jj = jj + m2
      jm = jm + inc
c
c***  Zero the upper left and lower right corners of E
c
   17 continue
      kl = 1
      n2m = (m2p1 * n) + 1
      do 19 i = 1, m
      ku = (kl + m) - i
      do 18 k = kl, ku
      we(k) = zero
      we(n2m - k) = zero
   18 continue
      kl = kl + m2p1
c
c***  Weighted matrix WE = W**-1 * E and its L1-norm
c
   19 continue
   20 jj = 0
      el = 0d0
      do 22 i = 1, n
      wi = w(i)
      do 21 j = 1, m2p1
      jj = jj + 1
      we(jj) = we(jj) / wi
      el = el + abs(we(jj))
   21 continue
   22 continue
c
c***  Ready
c
      el = el / n
      return 
c SPLC.FOR, 1985-12-12
c
c Author: H.J. Woltring
c
c Organizations: University of Nijmegen, and
c                Philips Medical Systems, Eindhoven
c                (The Netherlands)
c
c***********************************************************************
c
c FUNCTION SPLC (REAL*8)
c
c Purpose:
c *******
c
c       To assess the coefficients of a B-spline and various statistical
c       parameters, for a given value of the regularization parameter p.
c
c Calling convention:
c ******************
c
c       FV = SPLC ( M, N, K, Y, NY, WX, WY, MODE, VAL, P, EPS, C, NC,
c       1           STAT, B, WE, EL, BWE)
c
c Meaning of parameters:
c *********************
c
c       SPLC            ( O )   GCV function value if |MODE|.eq.2,
c                               MSE value if |MODE|.eq.3, and absolute
c                               difference with the prior given number o
cf
c                               degrees of freedom if |MODE|.eq.4.
c       M               ( I )   Half order of the B-spline (degree 2*M-1
c),
c                               with M > 0.
c       N               ( I )   Number of observations, with N >= 2*M.
c       K               ( I )   Number of datasets, with K >= 1.
c       Y(NY,K)         ( I )   Observed measurements.
c       NY              ( I )   First dimension of Y(NY,K), with NY.ge.N
c.
c       WX(N)           ( I )   Weight factors, corresponding to the
c                               relative inverse variance of each measur
ce-
c                               ment, with WX(I) > 0.0.
c       WY(K)           ( I )   Weight factors, corresponding to the
c                               relative inverse variance of each datase
ct,
c                               with WY(J) > 0.0.
c       MODE            ( I )   Mode switch, as described in GCVSPL.
c       VAL             ( I )   Prior variance if |MODE|.eq.3, and
c                               prior number of degrees of freedom if
c                               |MODE|.eq.4. For other values of MODE,
c                               VAL is not used.
c       P               ( I )   Smoothing parameter, with P >= 0.0. If
c                               P.eq.0.0, an interpolating spline is
c                               calculated.
c       EPS             ( I )   Relative rounding tolerance*10.0. EPS is
c                               the smallest positive number such that
c                               EPS/10.0 + 1.0 .ne. 1.0.
c       C(NC,K)         ( O )   Calculated spline coefficient arrays. NB
c:
c                               the dimensions of in GCVSPL and in SPLDE
cR
c                               are different! In SPLDER, only a single
c                               column of C(N,K) is needed, and the prop
cer
c                               column C(1,J), with J=1...K, should be u
csed
c                               when calling SPLDER.
c       NC              ( I )   First dimension of C(NC,K), with NC.ge.N
c.
c       STAT(6)         ( O )   Statistics array. See the description in
c                               subroutine GCVSPL.
c       B (1-M:M-1,N)   ( I )   B-spline tableau as evaluated by subrout
cine
c                               BASIS.
c       WE( -M:M  ,N)   ( I )   Weighted B-spline tableau (W**-1 * E) as
c                               evaluated by subroutine PREP.
c       EL              ( I )   L1-norm of the matrix WE as evaluated by
c                               subroutine PREP.
c       BWE(-M:M,N)     ( O )   Central 2*M+1 bands of the inverted
c                               matrix ( B  +  p * W**-1 * E )**-1
c
c Remarks:
c *******
c
c       This subroutine combines elements of subroutine SPLC0 from the
c       paper by Lyche et al. (1983), and of subroutine SPFIT1 by
c       Hutchinson (1985).
c
c References:
c **********
c
c       M.F. Hutchinson (1985), Subroutine CUBGCV. CSIRO division of
c       Mathematics and Statistics, P.O. Box 1965, Canberra, ACT 2601,
c       Australia.
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      function splc(m, n, k, y, ny, wx, wy, mode, val, p, eps, c, nc, 
     &stat, b, we, el, bwe)
      implicit double precision (o-z, a-h)
      parameter (zero = 0d0, one = 1d0, two = 2d0)
c
c***  Check on p-value
c
      dimension y(ny, k), wx(n), wy(k), c(nc, k), stat(6), b(1 - m:m - 1
     &, n), we(- m:m, n), bwe(- m:m, n)
      dp = p
      stat(4) = p
c***  Pseudo-interpolation if p is too small
      pel = p * el
      if (pel .lt. eps) then
      dp = eps / el
      stat(4) = zero
c***  Pseudo least-squares polynomial if p is too large
      end if
      if ((pel * eps) .gt. one) then
      dp = one / (el * eps)
      stat(4) = dp
c
c***  Calculate  BWE  =  B  +  p * W**-1 * E
c
      end if
      do 40 i = 1, n
      km = - min0(m,i - 1)
      kp = min0(m,n - i)
      do 30 l = km, kp
      if (iabs(l) .eq. m) then
      bwe(l,i) = dp * we(l,i)
      else
      bwe(l,i) = b(l,i) + (dp * we(l,i))
      end if
   30 continue
c
c***  Solve BWE * C = Y, and assess TRACE [ B * BWE**-1 ]
c
   40 continue
      call bandet(bwe, m, n)
      call bansol(bwe, y, ny, c, nc, m, n, k)
ctrace * p = res. d.o.
      stat(3) = trinv(we,bwe,m,n) * dp
c
c***  Compute mean-squared weighted residual
c
      trn = stat(3) / n
      esn = zero
      do 70 j = 1, k
      do 60 i = 1, n
      dt = - y(i,j)
      km = - min0(m - 1,i - 1)
      kp = min0(m - 1,n - i)
      do 50 l = km, kp
      dt = dt + (b(l,i) * c(i + l,j))
   50 continue
      esn = esn + (((dt * dt) * wx(i)) * wy(j))
   60 continue
   70 continue
c
c***  Calculate statistics and function value
c
      esn = esn / (n * k)
cEstimated variance               
      stat(6) = esn / trn
cGCV function value               
      stat(1) = stat(6) / trn
c     STAT(3) = trace [p*B * BWE**-1] !Estimated residuals' d.o.f.
c     STAT(4) = P                     !Normalized smoothing factor
cMean Squared Residual            
      stat(2) = esn
c***     Unknown variance: GCV
      if (iabs(mode) .ne. 3) then
      stat(5) = stat(6) - esn
      if (iabs(mode) .eq. 1) splc = zero
      if (iabs(mode) .eq. 2) splc = stat(1)
      if (iabs(mode) .eq. 4) splc = dabs(stat(3) - val)
c***     Known variance: estimated mean squared error
      else
      stat(5) = esn - (val * ((two * trn) - one))
      splc = stat(5)
c
      end if
      return 
c BANDET.FOR, 1985-06-03
c
c***********************************************************************
c
c SUBROUTINE BANDET (REAL*8)
c
c Purpose:
c *******
c
c       This subroutine computes the LU decomposition of an N*N matrix
c       E. It is assumed that E has M bands above and M bands below the
c       diagonal. The decomposition is returned in E. It is assumed that
c       E can be decomposed without pivoting. The matrix E is stored in
c       vectorized form in the array E(-M:M,N), where element E(J,I) of
c       the array E corresponds with element e(i,i+j) of the matrix E.
c
c Calling convention:
c ******************
c
c       CALL BANDET ( E, M, N )
c
c Meaning of parameters:
c *********************
c
c       E(-M:M,N)       (I/O)   Matrix to be decomposed.
c       M, N            ( I )   Matrix dimensioning parameters,
c                               M >= 0, N >= 2*M.
c
c Remark:
c ******
c
c       No checking on the validity of the input data is performed.
c       If (M.le.0), no action is taken.
c
c***********************************************************************
c
      end
c
      subroutine bandet(e, m, n)
      implicit double precision (o-z, a-h)
c
      dimension e(- m:m, n)
      if (m .le. 0) return 
      do 40 i = 1, n
      di = e(0,i)
      mi = min0(m,i - 1)
      if (mi .ge. 1) then
      do 10 k = 1, mi
      di = di - (e(- k,i) * e(k,i - k))
   10 continue
      e(0,i) = di
      end if
      lm = min0(m,n - i)
      if (lm .ge. 1) then
      do 30 l = 1, lm
      dl = e(- l,i + l)
      km = min0(m - l,i - 1)
      if (km .ge. 1) then
      du = e(l,i)
      do 20 k = 1, km
      du = du - (e(- k,i) * e(l + k,i - k))
      dl = dl - (e((- l) - k,l + i) * e(k,i - k))
   20 continue
      e(l,i) = du
      end if
      e(- l,i + l) = dl / di
   30 continue
      end if
c
c***  Ready
c
   40 continue
      return 
c BANSOL.FOR, 1985-12-12
c
c***********************************************************************
c
c SUBROUTINE BANSOL (REAL*8)
c
c Purpose:
c *******
c
c       This subroutine solves systems of linear equations given an LU
c       decomposition of the design matrix. Such a decomposition is pro-
c       vided by subroutine BANDET, in vectorized form. It is assumed
c       that the design matrix is not singular.
c
c Calling convention:
c ******************
c
c       CALL BANSOL ( E, Y, NY, C, NC, M, N, K )
c
c Meaning of parameters:
c *********************
c
c       E(-M:M,N)       ( I )   Input design matrix, in LU-decomposed,
c                               vectorized form. Element E(J,I) of the
c                               array E corresponds with element
c                               e(i,i+j) of the N*N design matrix E.
c       Y(NY,K)         ( I )   Right hand side vectors.
c       C(NC,K)         ( O )   Solution vectors.
c       NY, NC, M, N, K ( I )   Dimensioning parameters, with M >= 0,
c                               N > 2*M, and K >= 1.
c
c Remark:
c ******
c
c       This subroutine is an adaptation of subroutine BANSOL from the
c       paper by Lyche et al. (1983). No checking is performed on the
c       validity of the input parameters and data. Division by zero may
c       occur if the system is singular.
c
c Reference:
c *********
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      subroutine bansol(e, y, ny, c, nc, m, n, k)
      implicit double precision (o-z, a-h)
c
c***  Check on special cases: M=0, M=1, M>1
c
      dimension e(- m:m, n), y(ny, k), c(nc, k)
      nm1 = n - 1
c
c***  M = 0: Diagonal system
c
      if (m - 1) 10, 40, 80
   10 do 30 i = 1, n
      do 20 j = 1, k
      c(i,j) = y(i,j) / e(0,i)
   20 continue
   30 continue
c
c***  M = 1: Tridiagonal system
c
      return 
   40 do 70 j = 1, k
      c(1,j) = y(1,j)
cForward sweep                          
      do 50 i = 2, n
      c(i,j) = y(i,j) - (e(-1,i) * c(i - 1,j))
   50 continue
      c(n,j) = c(n,j) / e(0,n)
cBackward sweep                          
      do 60 i = nm1, 1, -1
      c(i,j) = (c(i,j) - (e(1,i) * c(i + 1,j))) / e(0,i)
   60 continue
   70 continue
c
c***  M > 1: General system
c
      return 
   80 do 130 j = 1, k
      c(1,j) = y(1,j)
cForward sweep                         
      do 100 i = 2, n
      mi = min0(m,i - 1)
      d = y(i,j)
      do 90 l = 1, mi
      d = d - (e(- l,i) * c(i - l,j))
   90 continue
      c(i,j) = d
  100 continue
      c(n,j) = c(n,j) / e(0,n)
cBackward sweep                         
      do 120 i = nm1, 1, -1
      mi = min0(m,n - i)
      d = c(i,j)
      do 110 l = 1, mi
      d = d - (e(l,i) * c(i + l,j))
  110 continue
      c(i,j) = d / e(0,i)
  120 continue
  130 continue
c
      return 
c TRINV.FOR, 1985-06-03
c
c***********************************************************************
c
c FUNCTION TRINV (REAL*8)
c
c Purpose:
c *******
c
c       To calculate TRACE [ B * E**-1 ], where B and E are N * N
c       matrices with bandwidth 2*M+1, and where E is a regular matrix
c       in LU-decomposed form. B and E are stored in vectorized form,
c       compatible with subroutines BANDET and BANSOL.
c
c Calling convention:
c ******************
c
c       TRACE = TRINV ( B, E, M, N )
c
c Meaning of parameters:
c *********************
c
c       B(-M:M,N)       ( I ) Input array for matrix B. Element B(J,I)
c                             corresponds with element b(i,i+j) of the
c                             matrix B.
c       E(-M:M,N)       (I/O) Input array for matrix E. Element E(J,I)
c                             corresponds with element e(i,i+j) of the
c                             matrix E. This matrix is stored in LU-
c                             decomposed form, with L unit lower tri-
c                             angular, and U upper triangular. The unit
c                             diagonal of L is not stored. Upon return,
c                             the array E holds the central 2*M+1 bands
c                             of the inverse E**-1, in similar ordering.
c       M, N            ( I ) Array and matrix dimensioning parameters
c                             (M.gt.0, N.ge.2*M+1).
c       TRINV           ( O ) Output function value TRACE [ B * E**-1 ]
c
c Reference:
c *********
c
c       A.M. Erisman & W.F. Tinney, On computing certain elements of the
c       inverse of a sparse matrix. Communications of the ACM 18(1975),
c       nr. 3, pp. 177-179.
c
c***********************************************************************
c
      end
c
      double precision function trinv(b, e, m, n)
      implicit double precision (o-z, a-h)
      parameter (zero = 0d0, one = 1d0)
c
c***  Assess central 2*M+1 bands of E**-1 and store in array E
c
      dimension b(- m:m, n), e(- m:m, n)
cNth pivot                             
      e(0,n) = one / e(0,n)
      do 40 i = n - 1, 1, -1
      mi = min0(m,n - i)
c***     Save Ith column of L and Ith row of U, and normalize U row
cIth pivot                             
      dd = one / e(0,i)
      do 10 k = 1, mi
cIth row of U (normalized)    
      e(k,n) = e(k,i) * dd
cIth column of L                   
      e(- k,1) = e(- k,k + i)
   10 continue
c***     Invert around Ith pivot
      dd = dd + dd
      do 30 j = mi, 1, -1
      du = zero
      dl = zero
      do 20 k = 1, mi
      du = du - (e(k,n) * e(j - k,i + k))
      dl = dl - (e(- k,1) * e(k - j,i + j))
   20 continue
      e(j,i) = du
      e(- j,j + i) = dl
      dd = dd - ((e(j,n) * dl) + (e(- j,1) * du))
   30 continue
      e(0,i) = 5d-1 * dd
c
c***  Assess TRACE [ B * E**-1 ] and clear working storage
c
   40 continue
      dd = zero
      do 60 i = 1, n
      mn = - min0(m,i - 1)
      mp = min0(m,n - i)
      do 50 k = mn, mp
      dd = dd + (b(k,i) * e(- k,k + i))
   50 continue
   60 continue
      trinv = dd
      do 70 k = 1, m
      e(k,n) = zero
      e(- k,1) = zero
c
c***  Ready
c
   70 continue
      return 
c SPLDER.FOR, 1985-06-11
c
c***********************************************************************
c
c FUNCTION SPLDER (REAL*8)
c
c Purpose:
c *******
c
c       To produce the value of the function (IDER.eq.0) or of the
c       IDERth derivative (IDER.gt.0) of a 2M-th order B-spline at
c       the point T. The spline is described in terms of the half
c       order M, the knot sequence X(N), N.ge.2*M, and the spline
c       coefficients C(N).
c
c Calling convention:
c ******************
c
c       SVIDER = SPLDER ( IDER, M, N, T, X, C, L, Q )
c
c Meaning of parameters:
c *********************
c
c       SPLDER  ( O )   Function or derivative value.
c       IDER    ( I )   Derivative order required, with 0.le.IDER
c                       and IDER.le.2*M. If IDER.eq.0, the function
c                       value is returned; otherwise, the IDER-th
c                       derivative of the spline is returned.
c       M       ( I )   Half order of the spline, with M.gt.0.
c       N       ( I )   Number of knots and spline coefficients,
c                       with N.ge.2*M.
c       T       ( I )   Argument at which the spline or its deri-
c                       vative is to be evaluated, with X(1).le.T
c                       and T.le.X(N).
c       X(N)    ( I )   Strictly increasing knot sequence array,
c                       X(I-1).lt.X(I), I=2,...,N.
c       C(N)    ( I )   Spline coefficients, as evaluated by
c                       subroutine GVCSPL.
c       L       (I/O)   L contains an integer such that:
c                       X(L).le.T and T.lt.X(L+1) if T is within
c                       the range X(1).le.T and T.lt.X(N). If
c                       T.lt.X(1), L is set to 0, and if T.ge.X(N),
c                       L is set to N. The search for L is facili-
c                       tated if L has approximately the right
c                       value on entry.
c       Q(2*M)  ( W )   Internal work array.
c
c Remark:
c ******
c
c       This subroutine is an adaptation of subroutine SPLDER of
c       the paper by Lyche et al. (1983). No checking is performed
c       on the validity of the input parameters.
c
c Reference:
c *********
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      double precision function splder(ider, m, n, t, x, c, l, q)
      implicit double precision (o-z, a-h)
      parameter (zero = 0d0, one = 1d0)
cf2py intent(out) :: splder
cf2py intent(hide) :: q
c
c***  Derivatives of IDER.ge.2*M are alway zero
c
      dimension x(n), c(n), q(2 * m)
      m2 = 2 * m
      k = m2 - ider
      if (k .lt. 1) then
      splder = zero
      return 
c
c***  Search for the interval value L
c
      end if
c
c***  Initialize parameters and the 1st row of the B-spline
c***  coefficients tableau
c
      call search(n, x, t, l)
      tt = t
      mp1 = m + 1
      npm = n + m
      m2m1 = m2 - 1
      k1 = k - 1
      nk = n - k
      lk = l - k
      lk1 = lk + 1
      lm = l - m
      jl = l + 1
      ju = l + m2
      ii = n - m2
      ml = - l
      do 2 j = jl, ju
      if ((j .ge. mp1) .and. (j .le. npm)) then
      q(j + ml) = c(j - m)
      else
      q(j + ml) = zero
      end if
c
c***  The following loop computes differences of the B-spline
c***  coefficients. If the value of the spline is required,
c***  differencing is not necessary.
c
    2 continue
      if (ider .gt. 0) then
      jl = jl - m2
      ml = ml + m2
      do 6 i = 1, ider
      jl = jl + 1
      ii = ii + 1
      j1 = max0(1,jl)
      j2 = min0(l,ii)
      mi = m2 - i
      j = j2 + 1
      if (j1 .le. j2) then
      do 3 jin = j1, j2
      j = j - 1
      jm = ml + j
      q(jm) = (q(jm) - q(jm - 1)) / (x(j + mi) - x(j))
    3 continue
      end if
      if (jl .ge. 1) goto 6
      i1 = i + 1
      j = ml + 1
      if (i1 .le. ml) then
      do 5 jin = i1, ml
      j = j - 1
      q(j) = - q(j - 1)
    5 continue
      end if
    6 continue
      do 7 j = 1, k
      q(j) = q(j + ider)
    7 continue
c
c***  Compute lower half of the evaluation tableau
c
      end if
cTableau ready if IDER.eq.2*M-1            
      if (k1 .ge. 1) then
      do 14 i = 1, k1
      nki = nk + i
      ir = k
      jj = l
      ki = k - i
c***        Right-hand B-splines
      nki1 = nki + 1
      if (l .ge. nki1) then
      do 9 j = nki1, l
      q(ir) = q(ir - 1) + ((tt - x(jj)) * q(ir))
      jj = jj - 1
      ir = ir - 1
    9 continue
c***        Middle B-splines
      end if
      lk1i = lk1 + i
      j1 = max0(1,lk1i)
      j2 = min0(l,nki)
      if (j1 .le. j2) then
      do 11 j = j1, j2
      xjki = x(jj + ki)
      z = q(ir)
      q(ir) = z + (((xjki - tt) * (q(ir - 1) - z)) / (xjki - x(jj)))
      ir = ir - 1
      jj = jj - 1
   11 continue
c***        Left-hand B-splines
      end if
      if (lk1i .le. 0) then
      jj = ki
      lk1i1 = 1 - lk1i
      do 13 j = 1, lk1i1
      q(ir) = q(ir) + ((x(jj) - tt) * q(ir - 1))
      jj = jj - 1
      ir = ir - 1
   13 continue
      end if
   14 continue
c
c***  Compute the return value
c
      end if
c***  Multiply with factorial if IDER.gt.0
      z = q(k)
      if (ider .gt. 0) then
      do 16 j = k, m2m1
      z = z * j
   16 continue
      end if
c
c***  Ready
c
      splder = z
      return 
c SEARCH.FOR, 1985-06-03
c
c***********************************************************************
c
c SUBROUTINE SEARCH (REAL*8)
c
c Purpose:
c *******
c
c       Given a strictly increasing knot sequence X(1) < ... < X(N),
c       where N >= 1, and a real number T, this subroutine finds the
c       value L such that X(L) <= T < X(L+1).  If T < X(1), L = 0;
c       if X(N) <= T, L = N.
c
c Calling convention:
c ******************
c
c       CALL SEARCH ( N, X, T, L )
c
c Meaning of parameters:
c *********************
c
c       N       ( I )   Knot array dimensioning parameter.
c       X(N)    ( I )   Stricly increasing knot array.
c       T       ( I )   Input argument whose knot interval is to
c                       be found.
c       L       (I/O)   Knot interval parameter. The search procedure
c                       is facilitated if L has approximately the
c                       right value on entry.
c
c Remark:
c ******
c
c       This subroutine is an adaptation of subroutine SEARCH from
c       the paper by Lyche et al. (1983). No checking is performed
c       on the input parameters and data; the algorithm may fail if
c       the input sequence is not strictly increasing.
c
c Reference:
c *********
c
c       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       Advances in Engineering Software 5(1983)1, pp. 2-5.
c
c***********************************************************************
c
      end
c
      subroutine search(n, x, t, l)
      implicit double precision (o-z, a-h)
c
      dimension x(n)
c***     Out of range to the left
      if (t .lt. x(1)) then
      l = 0
      return 
      end if
c***     Out of range to the right
      if (t .ge. x(n)) then
      l = n
      return 
c***  Validate input value of L
      end if
      l = max0(l,1)
c
c***  Often L will be in an interval adjoining the interval found
c***  in a previous call to search
c
      if (l .ge. n) l = n - 1
      if (t .ge. x(l)) goto 5
      l = l - 1
c
c***  Perform bisection
c
      if (t .ge. x(l)) return 
      il = 1
    3 iu = l
    4 l = (il + iu) / 2
      if ((iu - il) .le. 1) return 
      if (t .lt. x(l)) goto 3
      il = l
      goto 4
    5 if (t .lt. x(l + 1)) return 
      l = l + 1
      if (t .lt. x(l + 1)) return 
      il = l + 1
      iu = n
c
      goto 4
      end
