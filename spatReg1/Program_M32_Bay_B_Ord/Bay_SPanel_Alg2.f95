Module Precision
IMPLICIT none
 Integer, Parameter :: sp = Kind(0.0E0)
 Integer, Parameter :: dp = Kind(0.0D0)
 Integer, Parameter :: wp = dp
End Module Precision


MODULE DataParameters
 USE Precision
 USE F95_LAPACK, ONLY:                    LA_POSV, LA_SYSV, LA_POTRF
 USE F95_LAPACK, ONLY:                    LA_SYEV
 USE F77_LAPACK, ONLY:                    LA_GETRF
IMPLICIT none
 INTEGER, PARAMETER::                     nobs = 10766 
 INTEGER, PARAMETER::                     N = 2408 
 INTEGER, PARAMETER::                     TYear = 58 
 INTEGER, PARAMETER::                     KVars = 25 
 INTEGER, DIMENSION(nobs)::               Year
 INTEGER, DIMENSION(nobs) ::               CountryLocationID
 INTEGER, DIMENSION(nobs) ::               UrbanLocationID
 INTEGER, DIMENSION(nobs)::               CityID
 REAL(wp), DIMENSION(nobs)::              Y
 REAL(wp), DIMENSION(nobs,KVars)::        X
 INTEGER, DIMENSION(N) :: CountryID
 REAL(wp), DIMENSION(2, N)::              coordinates

 INTEGER, DIMENSION(TYear)::              YStart, YEnd
 REAL(wp), DIMENSION(N, N)::              W    ! weight matrix
 INTEGER, ALLOCATABLE, DIMENSION(:)::     CTObserved
 REAL(wp), ALLOCATABLE, DIMENSION(:,:)::  CYear, WYear, XYear, WXYear
 REAL(wp), ALLOCATABLE, DIMENSION(:)::    YYear, WYYear, RowSum
 REAL(wp), ALLOCATABLE, DIMENSION(:)::    eigsYear 
 REAL(wp), DIMENSION(nobs)::              eigs 
 REAL(wp), DIMENSION(nobs, KVars)::       WX
 REAL(wp), DIMENSION(nobs)::              WY

 INTEGER, PARAMETER::                     ndraw = 1000  ! number of draw
 INTEGER, PARAMETER::                     nburn = 0  ! number of burn-in
 INTEGER::                                iter
 REAL(wp), DIMENSION(KVars)::             beta
 REAL(wp)::                               sigma_sq
 REAL(wp)::                               rho
 REAL(wp), DIMENSION(N)::                 MU
 REAL(wp)::                               sigma_mu_sq

 INTEGER, PARAMETER::                     npar = KVars + 3 + N 
 REAL(wp), DIMENSION(ndraw-nburn, npar):: AllDraws
 REAL(wp), DIMENSION(KVars)::             XTY
 REAL(wp), DIMENSION(KVars, KVars)::      XTX
 REAL(wp), DIMENSION(nobs)::              e       ! OLS residuals
 REAL(wp), DIMENSION(KVars)::             beta0
 REAL(wp), DIMENSION(KVars, KVars)::      M0
 REAL(wp)::                               s0, v0
 REAL(wp)::                               h0, p0
 REAL(wp)::                               rho_min, rho_max
 REAL(wp)::                               cc
 INTEGER::                                acc

 REAL(wp), DIMENSION(nobs)::              YTilde
 REAL(wp), DIMENSION(nobs,KVars)::        XTilde
 REAL(wp), DIMENSION(KVars, KVars)::      identity_M1
 REAL(wp)::                               rhotemp, logdet
 REAL(wp)::                               lnconditionalrho

 REAL(wp), DIMENSION(npar)::              PostMeans, PostVars
 REAL(wp), DIMENSION(npar)::              NSE, RNE, CD
 REAL(wp)::                               frac1, frac2


END MODULE DataParameters


MODULE random
! This code a modified version of Alan Miller's random.f90.
! This module contains the random_normal and random_gamma
! distributions. I includes a code for choosing precision.

! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma

! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

!     Author: Alan Miller
!     e-mail: amiller @ bigpond.net.au

USE Precision
IMPLICIT NONE

REAL(wp), PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = TINY(1.0), vlarge = HUGE(1.0)


CONTAINS


FUNCTION random_normal() RESULT(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp) :: fn_val

!     Local variables
REAL(wp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal



FUNCTION random_gamma(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp), INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL(wp)                :: fn_val

IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
END IF

RETURN
END FUNCTION random_gamma



FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp), INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL(wp)                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s - one/3.
  c = one/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = random_normal()
    v = (one + c*x)**3
    IF (v > zero) EXIT
  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < one - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp), INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL(wp)                :: fn_val

!     Local variables
REAL       :: r, x, w
REAL, SAVE :: a, p, c, uf, vr, d

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2


FUNCTION random_exponential() RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp)  :: fn_val

!     Local variable
REAL(wp)  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

fn_val = -LOG(r)
RETURN

END FUNCTION random_exponential


END MODULE random



PROGRAM Bay_SPanel_Alg2 
USE Precision
USE random 
USE DataParameters
IMPLICIT none
CHARACTER(100)::              DataPath = "./../InData/"
CHARACTER(100)::              input1, input2 
INTEGER::                     i, j, t, id, ios
REAL::                        start, finish
!****************************************************************
!****************************************************************
!* This is Bayesian estimation of unbalanced panel data/spatial *
!* model using Metropolis-Within-Gibbs algorithm.               *
!* Ord's method is used for log-determinant.                    *
!****************************************************************
!****************************************************************

!*****************************************************************
! read data                                                      *
!*****************************************************************
!(1) import data
 input1 = adjustl(trim(DataPath))//"GrowthSpatialData_M32.raw"
 OPEN (Unit = 3, File = input1, ACTION="READ")
 DO i = 1, nobs
   READ (Unit = 3, Fmt = *, IOSTAT=ios) &
         Year(i), CountryLocationID(i), UrbanLocationID(i), CityID(i), &
                            Y(i), X(i,2:KVars)
 END DO
 CLOSE (Unit = 3)
 print *, "IO status for Y, X data is ",ios

! Data include constant terms
 X(:,1) = 1.0_wp

! (2) import latitude and logitude data 
 input2 = adjustl(trim(datapath))//"LatLong.raw"
 OPEN (Unit = 3, File = input2, ACTION="READ")
 DO i = 1, N
   READ (Unit = 3, Fmt = *, IOSTAT=ios) CountryID(i), coordinates(:,i)
 END DO
 CLOSE (Unit = 3)
 print *, "IO status for weight matrix is ",ios

open(unit=10, file='./TimeCheck_Alg2.raw', action="write")
!*************************************************************
! Find the row range for each year and the total             *
! number of years                                            *
!*************************************************************
YStart = 0
YEnd = 0

i = 1
t = 0
DO WHILE ( i <= nobs)
   id = Year(i)
   t = t + 1
   YStart(t) = i
   DO WHILE (Year(i) == id .and. i <= nobs)
    i = i + 1
    if (i==nobs+1) EXIT
   END DO
   YEnd(t) = i-1
END DO

!print *, "The total number of Year is", t


!***************************************************************
! Specify spatial weight matrix                                *
!***************************************************************
 CALL SWeightMatrix
 CALL WXWY
 CALL Prep_LogDet_Ord


!***************************************************************
! HyperParameters                                              *
!***************************************************************
do i =1, KVars
beta0(i) = 0.0_wp
end do

M0 = 0.0_wp
do i=1, KVars
 M0(i,i) = 0.0001_wp
end do

! diffuse priors for sigma_mu_sq and sigma_sq
v0 = 0.0_wp ; S0 = 0.0_wp
h0 = 0.0_wp ; p0 = 0.0_wp

rho_min = 0.0_wp
rho_max = 1.0_wp


!******************************************************************
! set initial values                                              *
!******************************************************************
! inital values
  CALL  OLS
  rho = 0.5_wp
  sigma_mu_sq = 100.0_wp
  do i = 1, N
   MU(i) = random_normal()*sqrt(sigma_mu_sq)
  end do

!*****************************************************************
! Preparations for Sampling                                      *
!*****************************************************************
identity_M1 = 0.0_wp
do i = 1, KVars
 identity_M1(i,i) = 1.0_wp
enddo

!Used in MH algorithm
cc = 1.0_wp
acc = 0



print *, " "
print *, "   ************************************************   "
print *, "   * The number of draws is ", ndraw
print *, "   * The average runtime is 108 min. per 100 draws"
print *, "   ************************************************   "
print *, " "

CALL CPU_TIME(start)

iter = 1
CALL RANDOM_SEED
DO WHILE (iter <= ndraw)

  CALL GibbsSampler

  ! Check time
   if ((iter==100) .or. (iter==200) .or. (iter==300) .or. &
      (iter==400) .or. (iter==500) .or. (iter==600) .or. &
      (iter==700) .or. (iter==800) .or. (iter==900) .or. &
      (iter==1000) ) then

      CALL CPU_TIME(finish)
      write(10,*) iter, " ", finish - start
   end if

 !save samples
   if (iter > nburn) then
      AllDraws(iter-nburn, 1:KVars) = beta
      AllDraws(iter-nburn, KVars+1) = sqrt(sigma_mu_sq)
      AllDraws(iter-nburn, KVars+2) = sqrt(sigma_sq)
      AllDraws(iter-nburn, KVars+3) = rho
      AllDraws(iter-nburn, KVars+4:N+KVars+3) = MU
   end if

iter = iter + 1
END DO

CALL CPU_TIME(finish)
  print *, "end of sampling: time", finish-start, "sec."

!******************************************************************
! save draws                                                      *
!******************************************************************
open(unit=16, file='./Draws_part1_bay_SPanel_alg2.raw', action="write")
open(unit=7, file='./Draws_part2_bay_SPanel_alg2.raw', action="write")
do i = 1, ndraw-nburn
  write(16,*) AllDraws(i, 1:KVars), AllDraws(i, KVars+1:KVars+3)
  write(7,*) AllDraws(i, KVars+4:KVars+3+N)
end do

close(unit=16); close(unit=7)


!*****************************************************************
! Geweke (1991)'s convergence dignostics                         *
!***************************************************************** 
CALL MC_Statistics(AllDraws, ndraw-nburn, npar, PostMeans, &
                   PostVars, NSE, RNE)

frac1 = 0.1_wp
frac2 = 0.5_wp

CALL Geweke_Diagnostic(AllDraws, ndraw-nburn, npar, frac1, frac2, CD)
  write(*,*) "Geweke_Diag. done"

do i = 1, KVars+3
write(*,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
               sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do


open(unit=3, file='./GewekeResults_bay_spanel_alg2', action="write")
write(3,*) "The number of draws", ndraw
write(3,*) "The number of burns", nburn
write(3,*) "Post.Mean, Post.S.Error, NSE, RNE, CD"
write(3,*) "beta"
do i = 1, KVars 
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
               sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do
write(3,*) " "
write(3,*) "sigma_u, sigma_e and rho"
do i = KVars+1, KVars+3
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
                 sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do
write(3,*) " "
write(3,*) "Random effects u"
do i = Kvars+4, KVars + 3 + N
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
               sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do

close (unit = 3)

END PROGRAM Bay_SPanel_Alg2



SUBROUTINE GibbsSampler
USE Precision
USE DataParameters 
USE Random
IMPLICIT none
REAL(wp), DIMENSION(nobs)::          YYTilde, ZTilde
REAL(wp), DIMENSION(nobs,KVars)::    XXTilde
REAL(wp), DIMENSION(N)::             YTildeBar
REAL(wp), DIMENSION(N, KVars)::      XTildeBar
INTEGER, DIMENSION(N) ::             NObsCity
REAL(wp)::                           theta 
REAL(wp), DIMENSION(KVars, KVars)::  M1
REAL(wp), DIMENSION(KVars, KVars)::  inverseM1
REAL(wp), DIMENSION(KVars)::         Term1
REAL(wp), DIMENSION(KVars)::         str
REAL(wp), DIMENSION(KVars)::         beta1
CHARACTER(5)::                       UPLO
REAL(wp)::                           p1, h1
LOGICAL ::                           First = .true.
REAL(wp)::                           v1, s1
REAL(WP), DIMENSION(nobs)::          BQ
REAL(wp), DIMENSION(KVars)::         Mb
REAL(wp)::                           inversePsi, DTBQ, MU1
REAL(wp)::                           rhostar
INTEGER::                            accept
REAL(wp)::                           lnp, lnpstar, u
REAL(wp)::                           rnd
REAL(wp), DIMENSION(ndraw)::         acc_rate
REAL(wp)::                           prob
INTEGER::                            i, j
 

!****************************************
! sampling beta
! Normal-gamma prior for beta and sigma^2
!****************************************
! first transformation
 XTilde = X - rho*WX
 YTilde = Y - rho*WY

! second transformation
Do i = 1, N
 NObsCity(i) = COUNT(CityID == i)
 YTildeBar(i) = SUM( PACK(YTilde, CityID == i) )/REAL(NObsCity(i))
 do j = 1, KVars
  XTildeBar(i,j) = SUM( PACK(XTilde(:,j), CityID == i), 1)/ &
                    REAL(NObsCity(i))
 end do
END DO

 DO i = 1, nobs
  theta = 1.0_wp - sqrt(sigma_sq/ &
         ( (real(NObsCity(CityID(i))))*sigma_mu_sq + sigma_sq ))
  YYTilde(i) = YTilde(i) - theta * YTildeBar(CityID(i))
  XXTilde(i,:) = XTilde(i,:) - theta *  XTildeBar(CityID(i),:)
 END DO


! used in sigma_sq 
 DO i = 1, nobs 
  ZTilde(i) = YTilde(i) - MU(CityID(i))
 END DO

DO i=1, KVars
  DO j=1, KVars
     XTX(i,j)=dot_product(XXTilde(:,i), XXTilde(:,j))
  END DO 
  XTY(i)=dot_product(XXTilde(:,i), YYTilde)
END DO

  M1 = M0 + XTX 
  Term1 = matmul(M0, beta0) + XTY

 !invert M1
  inverseM1 = identity_M1
  call la_posv(M1, inverseM1)

  beta1 = matmul(inverseM1, Term1)


! generate random sample from N(0,1)
  do i = 1, KVars
     str(i) = random_normal()
  end do

 ! cholesky decomposiiton of inverseM1 using la_potrf
  ! make upper part zero
    call la_potrf(inverseM1, UPLO = 'L')
    do i = 1, KVars 
      do j = i+1, KVars 
        inverseM1(i,j) = 0.0_wp
      end do
    end do

    beta =  sqrt(sigma_sq)*matmul(inverseM1, str) +  beta1

 
 !***************************************************
 ! sampling u
 !***************************************************
  BQ = YTilde - matmul(XTilde, beta)
  DO i = 1, N
   DTBQ = SUM(PACK(BQ, CityID == i))
  
   inversePsi = (sigma_mu_sq*sigma_sq)/(sigma_mu_sq*NObsCity(i) + sigma_sq)
   MU1 = (sigma_mu_sq*DTBQ)/(sigma_mu_sq*NObsCity(i) + sigma_sq)

   MU(i) = sqrt(inversePsi)*random_normal() + MU1

  END DO


 !*****************************************
 ! sampling sigma_e2
 !*****************************************
  v1 = real(nobs) + v0
  e = ZTilde - matmul(XTilde, beta)
  Mb = matmul(M0, beta-beta0)
  s1 = dot_product(e, e) + &
       dot_product(beta-beta0, Mb) + s0 

  sigma_sq = s1/(2.0_wp * random_gamma1(v1/2.0_wp, First))


 !**************************************************
 ! sampling sigma_u2
 !**************************************************
  h1 = h0 + real(N)
  p1 = dot_product(MU, MU) + p0

  sigma_mu_sq = p1/(2.0_wp * random_gamma1(h1/2.0_wp, First))


 !***************************************************************** 
 ! Sampling rho with Metropolis-Hastings Algorithm                *
 !*****************************************************************
  !(1) Take a cadidate from the candidate density, N(rho, cc^2) 
   rhostar = rho + cc * random_normal()
     accept = 0
     do while (accept == 0)
       if (( rhostar > rho_min).and.(rhostar < rho_max)) then
          accept = 1
       else
          rhostar = rho + cc * random_normal()
       end if
     end do

  !(2) Calculate ln p(rhostar) and ln p(rho)
       rhotemp = rho
       call WEIGHT
       lnp = lnconditionalrho

       rhotemp = rhostar
       call WEIGHT
       lnpstar = lnconditionalrho

  !(3) Accept the candidate with probability p(rhostar)/p(rho)
      call random_number(u)

      if ( (lnpstar - lnp) < 0.0_wp ) then
         prob = exp(lnpstar - lnp)
      else
         prob = 1.0_wp
      end if

     if (u < prob) then
        rho = rhostar
        acc = acc + 1
     end if

     acc_rate(iter) = real(acc)/real(iter)

 !(4) update cc depending on the acceptance rate 
    if (acc_rate(iter) < 0.4_wp) then
        cc = cc/1.1_wp
    else if (acc_rate(iter) > 0.6_wp) then
        cc = cc*1.1_wp
    end if

END SUBROUTINE GibbsSampler



SUBROUTINE WEIGHT 
USE Precision
USE DataParameters
IMPLICIT none
INTEGER::                              i, j
!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************
  CALL LogDet_Ord  
  

!*****************************************************************
! value of log of the full conditional distribution of rho 
! under noninformative prior
! It is used to calculate acceptance probability
! double-check the use of this function 
!*****************************************************************
XTilde = X - rhotemp*WX
YTilde = Y - rhotemp*WY
DO i = 1, nobs
 YTilde(i) = YTilde(i) - MU(CityID(i))
END DO
e = YTilde - matmul(XTilde, beta)


!*****************************************************************
! log of conditional distribution of rho 
!*****************************************************************
lnconditionalrho = logdet - dot_product(e, e)/(2.0_wp*sigma_sq)
             

RETURN
END SUBROUTINE WEIGHT



SUBROUTINE LogDet_Ord
use Precision
USE DataParameters
IMPLICIT none
INTEGER ::           i, j
REAL(wp), DIMENSION(nobs) :: detp
!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************
  detp = 1.0_wp - rhotemp*eigs

  logdet = sum( log(detp) )

RETURN
END SUBROUTINE LogDet_Ord



SUBROUTINE OLS
USE Precision
USE DataParameters
IMPLICIT none
INTEGER::          i, j
!*****************************************************************
! Solve for OLS estimator of beta, using the LU decomposition of
! X'X and right-hand side vector X'Y. Make use of the columns of
! the X matrix and dot-products for efficiency.
!*****************************************************************
  do i = 1, KVars
    do j = 1, KVars
      XTX(i,j) = dot_product(X(:,i),X(:,j))
    enddo
    XTY(i) = dot_product(X(:,i),Y)
  enddo

  call LA_POSV(XTX,XTY)
  beta = XTY
  print *, "OLS estimators for beta is", beta

  e = Y - matmul(X, beta)
  sigma_sq = dot_product(e, e)/real(nobs)
  print *, "OLS sigma squared is ", sigma_sq

RETURN
END SUBROUTINE OLS



SUBROUTINE SWeightMatrix
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
REAL(wp) ::                 toradians, r
!*************************************************************
!* This specifies spatial weight matrix by using geographic  *
!* latitude-logitude coordinates of spatial units.           *
!*************************************************************

!*************************************************************
! Calculate the distance between two spatial units           *
!*************************************************************
! Transform latitude and longitude from degrees to radians
toradians = 57.29577951_wp
coordinates = coordinates/toradians

! Use the value of r that yields kilometers (at the latitude
! of the first admin area, since the right value for the
! radius varies with latitude)
! One reference suggested using
 r = 6378d0 - 21d0 * sin(coordinates(1,1))
! The following radius is what is used in the Stata add-on
! program globdist, and gives results that agree with the
! Stata program up to the second or third decial
! (for Malawi)

!r = 6365.0_wp

! Implement the Haversine great circle distance formula
!  Note: Can improve code by going down columns of dist
DO i = 1, N
 DO j = 1, N
  W(i,j) = r * 2.0_wp *   &
  asin( min(1.0_wp, &
  sqrt( sin((coordinates(1,j)- coordinates(1,i))/2.0_wp)**2 &
   + cos(coordinates(1,i)) * cos(coordinates(1,j)) *    &
  sin( (coordinates(2,j)-coordinates(2,i))/2.0_wp )**2  &
  ) ))

  W(i,i) = 0.0_wp
 END DO
END DO

!print *, "max distance is", maxval(W(1,:))


!*************************************************************
! Specification: distance weight matrix                      *
!*************************************************************
DO i = 1, N
  where (W(:,i) > 0.0_wp)
     W(:,i) = 1.0_wp/W(:,i)
  end where
END DO

END SUBROUTINE SWeightMatrix



SUBROUTINE WXWY
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
!*************************************************************
! WX and WY                                                  *
!*************************************************************

DO t = 1, TYear

 ALLOCATE( CTObserved(YEnd(t)-YStart(t)+1) )
 ALLOCATE( CYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( WYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( XYear(size(CTObserved), KVars) )
 ALLOCATE( WXYear(size(CTObserved), KVars) )
 ALLOCATE( YYear(size(CTObserved)) )
 ALLOCATE( WYYear(size(CTObserved)) )
 ALLOCATE( RowSum(size(CTObserved)) )

 CTObserved = CityID ( YStart(t):YEnd(t) )
 XYear = X( YStart(t):YEnd(t), :)
 YYear = Y( YStart(t):YEnd(t))
 CYear = W( CTObserved, CTObserved )

!*************************************************************
! row standardization for each year.                         *
!*************************************************************
! Inverse of row sum calculated by column sum since CYear
! is symmetric.
DO i = 1, size(CTObserved)
  RowSum(i) = 1.0_wp/SUM(CYear(:,i))
END DO

! element-by-element multiplication 
 DO i = 1, size(CTObserved)
    WYear(:,i) = CYear(:,i)*RowSum
 END DO

DO i = 1, size(CTObserved)
  DO j = 1, KVars
     WXYear(i,j) = dot_product( WYear(i,:), XYear(:,j) )
  END DO
  WYYear(i) = dot_product( WYear(i,:), YYear)
END DO

WX(YStart(t):YEnd(t), :) = WXYear
WY(YStart(t):YEnd(t)) = WYYear


DEALLOCATE (CTObserved, CYear, WYear, XYear, WXYear, YYear, WYYear, RowSum)
END DO

END SUBROUTINE WXWY


SUBROUTINE Prep_LogDet_Ord
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
!*************************************************************
! Eigenvalues of row-standardized weight matrix              *
!*************************************************************

DO t = 1, TYear

 ALLOCATE( CTObserved(YEnd(t)-YStart(t)+1) )
 ALLOCATE( CYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( WYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( eigsYear(size(CTObserved)) )
 ALLOCATE( RowSum(size(CTObserved)) )

 CTObserved = CityID ( YStart(t):YEnd(t) )
 CYear = W( CTObserved, CTObserved )

!*************************************************************
! row standardization for each year.                         *
!*************************************************************
! Inverse of row sum calculated by column sum since CYear
! is symmetric.
DO i = 1, size(CTObserved)
  RowSum(i) = 1.0_wp/SUM(CYear(:,i))
END DO

! element-by-element multiplication 
 DO i = 1, size(CTObserved)
    WYear(:,i) = CYear(:,i)*RowSum
 END DO

!*************************************************************
! Create the symmetric eigenvalue equvalent matrix to the    *
! row-stardardized matrix.                                   *
!*************************************************************
 RowSum = SQRT(RowSum)

! Element-by-element multiplication.
! Here, WYear is no longer the row-stardardized matrix.
Do i = 1, size(CTObserved)
  WYear(:,i) = CYear(:,i)*RowSum
END DO

! Here, WYear is the symmetrix eigenvalue equivalent matrix.
DO i = 1, size(CTObserved)
  WYear(:,i) = WYear(:,i)*RowSum(i)
END DO

!LA_SYEV is used for a symmetric matrix.
!WYear is corrupted.
CALL LA_SYEV(WYear, eigsYear)
eigs(YStart(t):YEnd(t)) = eigsYear

DEALLOCATE (CTObserved, CYear, WYear, RowSum, eigsYear)
END DO

print *, "max eigenvalue of W is ", maxval(eigs)
print *, "min eigenvalue of W is ", minval(eigs)

END SUBROUTINE Prep_LogDet_Ord



SUBROUTINE Geweke_Diagnostic(Draws, ndraws, npars, frac1, frac2, CD )
USE Precision
IMPLICIT none
INTEGER, INTENT(IN)::               ndraws, npars
REAL(wp), INTENT(IN)::              frac1, frac2
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN):: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT):: CD
REAL(wp), DIMENSION(npars)::        PostMeans_a, PostMeans_b
REAL(wp), DIMENSION(npars)::        PostVars_a, PostVars_b
REAL(wp), DIMENSION(npars)::        NSE_a, NSE_b
REAL(wp), DIMENSION(npars)::        RNE_a, RNE_b
INTEGER::                           ndraws_a, ndraws_b
INTEGER::                           start_a, start_b
INTEGER::                           i

!********************************************************************
! Using MCMC draws, it returns Geweke (1992)'s convergence          * 
! diagonstic, CD.                                                   *
! Arguments                                                         *
!  Draws: (ndraws * npars) matrix of MCMC draws                     *
!  ndraws: the number of draws                                      *
!  npars: the number of parameters                                  *
!  frac1: first fraction of draws used for CD                       *
!  frac2: second fraction of draws used for CD                      * 
!  CD: convergence diagonstic statistic                             *
!********************************************************************

ndraws_b = floor(frac2*real(ndraws))
start_b =  ndraws - ndraws_b + 1
CALL MC_Statistics(Draws(start_b:ndraws,:), ndraws_b, npars, &
                   PostMeans_b, PostVars_b, NSE_b, RNE_b)

ndraws_a = floor(frac1*real(ndraws))
start_a = 1 
CALL MC_Statistics(Draws(start_a:ndraws_a,:), ndraws_a, npars, &
                   PostMeans_a, PostVars_a, NSE_a, RNE_a)

CD = (PostMeans_a - PostMeans_b)/(NSE_a + NSE_b)

END SUBROUTINE


SUBROUTINE MC_Statistics(Draws, ndraws, npars, PostMeans, PostVars, NSE, RNE)
USE Precision
IMPLICIT none
INTEGER, INTENT(IN) ::                            ndraws, npars
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN) :: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT) ::        PostMeans, PostVars
REAL(wp), DIMENSION(npars), INTENT(OUT)::         NSE, RNE
REAL(wp), ALLOCATABLE, DIMENSION(:,:)::           Covars
REAL(wp), DIMENSION(npars) ::                     SG 
INTEGER ::                                        lags
INTEGER ::                                        i, j

!***************************************************************************
! Using MCMC draws, it returns posterior mean, variance, NSE (numerical    *
! Standard error, and RNE (relative numerical error).                      *
! Arguments                                                                *
!  Draws: (ndraws * npars) matrix of draws                                 *
!  ndraws: the number of draws                                             *
!  npars: the number of parameters                                         *
!  PostMeans: posterior means of the parameters                            *
!  PostVars: posterior variances of the parameters                         *
!  NSE: numerical standard error                                           *
!  RNE: relative numerical error                                           *
!***************************************************************************

! posterior means 
PostMeans = SUM(Draws,1)/real(ndraws)

! posterior variances
do i = 1, npars
 PostVars(i) = SUM( (Draws(:,i)**2) )/real(ndraws) - PostMeans(i)**2
end do
! variance of posterior mean is SG/ndraws
! numerical standard error (nse) is dsqrt(SG/ndraws)
! relative numerical error (rne) is PostVar/SG
lags = floor( real(ndraws)**(0.25_wp) ) + 1
Allocate( Covars(lags, npars) )

 do i = 1, npars
  do j = 1, lags 
    Covars(j, i) = &
        dot_product( Draws(1:ndraws-j, i) - PostMeans(i), &
                     Draws(j+1:ndraws, i) - PostMeans(i) )/ &
        real(ndraws)
  end do
 end do  

 do i = 1, npars
  SG(i) = PostVars(i) + ( 2.0_wp*dot_product( real((/(i,i=lags,1,-1)/)), &
                       Covars(:,i)) )/ real(lags+1) 
 end do

 RNE = PostVars/SG
 NSE = sqrt( SG/real(ndraws) ) 

Deallocate ( Covars )

END SUBROUTINE 

SUBROUTINE GR_Diagnostic(draws, ndraws, npars, nchains, B, W, V, R)
USE Precision
IMPLICIT none
INTEGER, INTENT(IN):: ndraws, npars, nchains
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN):: draws
REAL(wp), DIMENSION(npars), INTENT(OUT):: B, W, V, R
REAL(wp), DIMENSION(npars) :: TotalMeans, ChainMeans
INTEGER :: N, start, end
REAL(wp), DIMENSION(ndraws, npars):: Deviation_sq
INTEGER:: i, j
!*****************************************************************
!Gelman-Rubin (1992)'s convergence diagnostic                    *
!Input arguments:                                                *
! draws: (nchain*ndraws by npars) matrix                         *
! ndraws: number of draws per chain                              *
! npars: number of parameter                                     *
! nchains: number of chains                                      *
!Output arguments:                                               *
! B: Between chain variance                                      *
! W: Within chain variance                                       *
! V: Overall estimate of variance                                *
! R: GR statistic                                                *
!***************************************************************** 

! posterior means
TotalMeans = SUM(draws,1)/real(ndraws)

! N: number of draws per chain
N = int(real(ndraws)/real(nchains))
start = 1
end = N
B = 0.0_wp

do i = 1, nchains
  ChainMeans = SUM(draws(start:end,:),1)/real(N)
  do j = 1, npars
    Deviation_sq(start:end, j) = (draws(start:end, j) - ChainMeans(j))**2
  end do
  B = B + (ChainMeans - TotalMeans)**2

  start = start + N
  end = end + N
end do

B = (B*real(N))/real(nchains-1)

W = SUM( Deviation_sq, 1)
W = W/(real(nchains)*real(N-1))

do i = 1, ndraws
end do
V = (1.0_wp - 1.0_wp/real(N))*W + B/real(N)

R = sqrt(V/W)

END SUBROUTINE   

