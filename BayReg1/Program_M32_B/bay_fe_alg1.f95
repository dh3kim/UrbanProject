MODULE GlobalVars 
USE LA_PRECISION, ONLY:                 wp => dp 
USE F95_LAPACK, ONLY:                   LA_POSV, LA_POTRF
USE F77_LAPACK, ONLY:                   LA_GETRF
IMPLICIT none
INTEGER, PARAMETER::                    nobs = 10766
INTEGER, PARAMETER::                    KVars = 17 
INTEGER, PARAMETER::                    NGroups = 2408 
INTEGER, DIMENSION(nobs)::              UrbanID        
INTEGER, DIMENSION(nobs)::              Year           
REAL(wp), DIMENSION(nobs)::             Y 
REAL(wp), DIMENSION(nobs, KVars)::      X
REAL(wp), ALLOCATABLE, DIMENSION(:)::   YN
REAL(wp), ALLOCATABLE, DIMENSION(:,:):: XN
INTEGER, DIMENSION(NGroups)::           T
REAL(wp), DIMENSION(NGroups)::          Ybar
REAL(wp), DIMENSION(NGroups, KVars)::   Xbar

INTEGER, DIMENSION(nobs) ::             pid 
INTEGER, DIMENSION(NGroups)::           CTStart, CTEnd

REAL(wp), DIMENSION(KVars) ::           beta, beta0,  beta1
REAL(wp), DIMENSION(NGroups) ::         u
LOGICAL ::                              First = .true.
REAL(wp) ::                             sigma_u_sq, sigma_e_sq
REAL(wp), DIMENSION(KVars, KVars) ::    M0
REAL(wp) ::                             v0, s0
REAL(wp), DIMENSION(KVars, KVars) ::    invM1, idenM1

INTEGER, PARAMETER::                    ndraw = 50000 
INTEGER, PARAMETER::                    nburn = 40000 
INTEGER, PARAMETER::                    npar = KVars + 1 + NGroups
REAL(wp), DIMENSION(ndraw-nburn, npar):: AllDraws

REAL(wp), DIMENSION(npar)::             PostMeans, PostVars
REAL(wp), DIMENSION(npar)::             NSE, RNE, CD
REAL(wp)::                              frac1, frac2

END MODULE GlobalVars 


MODULE Random
CONTAINS

FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
REAL  :: fn_val

! Local variables

REAL(8)            :: u, sum
REAL(8), SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL(8), PARAMETER :: one = 1.0, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
  fn_val = u*sln
END IF

RETURN
END FUNCTION rnorm

FUNCTION random_gamma1(s, first) RESULT(fn_val)
!************************************************************
! I change random Normal generator
! x = rnorm() instead of x = random_normal()
!***********************************************************

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

REAL(8), INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL(8)                :: fn_val

! Local variables
REAL, SAVE  :: c, d
REAL        :: u, v, x

IF (first) THEN
  d = s -1.0/3.0
  c = 1.0/SQRT(9.0*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = rnorm()
    v = (1.0 + c*x)**3.0
    IF (v > 0.0) EXIT

  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1.0 - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(1.0 - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1

END MODULE Random


PROGRAM bayes_fe
USE GlobalVars
USE Random
IMPLICIT none
CHARACTER(100) ::                  DataPath = "../InData/"
CHARACTER(100) ::                  CityData
REAL ::                            start, finish
INTEGER ::                         i, k, id, ios
INTEGER ::                         iter

!*****************************************************************
! read data                                                      *
!*****************************************************************
  CityData = adjustl(trim(DataPath))//"GrowthData_M32_FE.raw"
  open(unit=3, file=CityData, action="read")

  i=1
  DO WHILE (i <= nobs)
    read(3,*, IOSTAT=ios) UrbanID(i), Year(i), Y(i), X(i, 1:KVars)
    if (ios /= 0) EXIT
    i=i+1
  END DO
  print *, "ios is ", ios
  close(unit=3)

! intercept (excluded in this fixed-effects model)
! X(:,1) = 1.0_wp


!******************************************************************
! Find the total number of cities and the row range for each city *
!******************************************************************
i = 1
k = 0
DO WHILE ( i <= nobs)
   id = UrbanID(i)
   k = k + 1
   CTStart(k) = i
   DO WHILE (UrbanID(i) == id .and. i <= nobs)
    pid(i) = k
    i = i + 1
    if (i == nobs+1) EXIT
   END DO
   CTEnd(k) = i-1
END DO

!******************************************************************
! set hyperparameters in priors                                   *
!******************************************************************
! beta ~ N(beta0, MO)
do i = 1, KVars 
 beta0(i) = 0.0_wp
end do

M0 = 0.0_wp
do i=1, KVars 
 M0(i,i) = 0.000001_wp
end do

! u ~ N(0, sigma_u_sq) where sigma_u_sq is a hyperparameter.
sigma_u_sq = 100000.0_wp

! sigma ~ iG(v0/2, s0/2)
v0 = 0.0_wp ; s0 = 0.0_wp

!******************************************************************
! set initial values                                              *
!******************************************************************
CALL OLS

!sigma_u_sq = 1.25_wp

do i = 1, NGroups
 u(i) = rnorm()*5.0_wp
end do

! Preps for Sampling
 ! identity matrix
 idenM1 = 0.0_wp
 do i = 1, KVars
  idenM1(i,i) = 1.0_wp
 end do

 DO i = 1, NGroups
    ALLOCATE( YN(CTEnd(i)-CTStart(i)+1) )
    ALLOCATE( XN(size(YN), KVars) )

    YN = Y(CTStart(i):CTEnd(i))
    XN = X(CTStart(i):CTEnd(i),:)

    T(i) = SIZE(YN)
    Ybar(i) = SUM(YN)/DBLE(SIZE(YN))
    Xbar(i,:) = SUM(XN,1)/DBLE(SIZE(YN))

    DEALLOCATE( YN, XN )
 END DO

First = .true.

print *, " "
print *, "   ************************************************   "
print *, "   * The number of draws is ", ndraw
print *, "   * The average runtime is 116.31 sec. per 100 draws"
print *, "   ************************************************   "
print *, "Sampling starts ..."

CALL CPU_TIME(start)

iter = 1
CALL RANDOM_SEED
DO WHILE (iter <= ndraw)

CALL GibbsSampler

 ! save samples
  if (iter > nburn) then
     AllDraws(iter-nburn, 1:KVars) = beta
     !AllDraws(iter-nburn, KVars+1) = sqrt(sigma_u_sq)
     AllDraws(iter-nburn, KVars+1) = sqrt(sigma_e_sq)
     AllDraws(iter-nburn, KVars+2:NGroups+KVars+1) = u
  end if

iter = iter + 1
END DO

CALL CPU_TIME(finish)
  print *, "end of sampling: time", finish-start, "sec."

!******************************************************************
! save draws of some important parameters for graphs              *
!******************************************************************
open(unit=16, file='./Draws_part1_bay_fe_alg1.raw', action="write")
open(unit=7, file='./Draws_part2_bay_fe_alg1.raw', action="write")
do i = 1, ndraw-nburn
  write(16,*) AllDraws(i, 1:KVars), AllDraws(i, KVars+1:KVars+1)
  write(7,*) AllDraws(i, KVars+2:npar)
end do
close(unit=16)
close(unit=7)



CALL MC_Statistics(AllDraws, ndraw-nburn, NPar, PostMeans, PostVars, NSE, RNE)

frac1 = 0.1_wp
frac2 = 0.5_wp

CALL Geweke_Diagnostic(AllDraws, ndraw-nburn, NPar, frac1, frac2, CD)

open(unit=3, file='./GewekeResults_bay_fe_alg1', action="write")
write(3,*) "The number of draws", ndraw
write(3,*) "The number of burns", nburn
write(3,*) "Post.Mean,t-stat, Post.S.Error, NSE, RNE, CD"
write(3,*) "beta"
do i = 1, KVars
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
            sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do
write(3,*) " "
do i = KVars +1, KVars+2
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
            sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do
write(3,*) " "
do i = KVars + 3, npar
write(3,*) PostMeans(i), PostMeans(i)/sqrt(PostVars(i)), &
            sqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do
close (unit=3)


END PROGRAM bayes_fe




SUBROUTINE OLS
USE GlobalVars 
IMPLICIT none
INTEGER::          i, j
REAL(wp), DIMENSION(KVars, KVars) :: XTX
REAL(wp), DIMENSION(KVars) :: XTY
REAL(wp), DIMENSION(nobs) :: e
!*****************************************************************
! Solve for OLS estimator of beta, using the LU decomposition of
! X'X and right-hand side vector X'Y. Make use of the columns of
! the X matrix and dot-products for efficiency.
!*****************************************************************
  do i = 1, KVars
   do j = 1, KVars
    XTX(i,j) = dot_product(X(:,i), X(:,j))
   end do
   XTY(i) = dot_product(X(:,i), Y)
  enddo

  call LA_POSV(XTX,XTY)
  beta1 = XTY
  print *, ""
  print *, "OLS estimators for beta is"
  do i = 1, KVars
  print *,  beta1(i)
  end do

  e = Y - matmul(X, beta1)
  sigma_e_sq = dot_product(e, e)/real(nobs-KVars)
  print *, "OLS sigma is ", sqrt(sigma_e_sq)
  print *, ""

RETURN
END SUBROUTINE OLS




SUBROUTINE GibbsSampler
USE GlobalVars
USE Random
IMPLICIT none
INTEGER :: i, j, iter
REAL(wp), DIMENSION(nobs) :: W
REAL(wp), DIMENSION(KVars) :: XTW
REAL(wp), DIMENSION(KVars, KVars) :: XTX
REAL(wp), DIMENSION(KVars, KVars) :: M1
REAL(wp), DIMENSION(KVars) :: Term1
REAL(wp), DIMENSION(KVars) :: str
CHARACTER(5) :: UPLO
!REAL(wp) :: p1, h1
REAL(wp):: v1, s1
REAL(wp), DIMENSION(nobs):: epsilon
!REAL(wp), DIMENSION(KVars):: Mb
REAL(wp):: VarU, MeanU
REAL(wp):: theta
REAL(wp), DIMENSION(nobs):: YTilde
REAL(wp), DIMENSION(nobs, KVars):: XTilde


 !****************************************
 ! sampling beta
 ! Norma-gamma prior for beta and sigma^2
 !****************************************
   do i = 1, nobs
     W(i) = Y(i) - u(pid(i))
   end do 

  do i = 1, KVars 
    do j = 1, KVars 
       XTX(i,j) = dot_product(X(:,i), X(:,j))
    end do
    XTW(i) = dot_product(X(:,i), W)
  end do
 
  
  M1 = M0 + XTX 
  Term1 = matmul(M0, beta0) + XTW

 !invert M1
  invM1 = idenM1
  call la_posv(M1, invM1)

  beta1 = matmul(invM1, Term1)


! generate random sample from N(0,1)
  do i = 1, KVars
     str(i) = rnorm()
  end do

 ! cholesky decomposiiton of inverseM1 using la_potrf
  ! make upper part zero
    call la_potrf(invM1, UPLO = 'L')
    do i = 1, KVars 
      do j = i+1, KVars 
        invM1(i,j) = 0.0_wp
      end do
    end do

    beta =  dsqrt(sigma_e_sq)*matmul(invM1, str) +  beta1

   !*****************************************
   ! sampling sigma_e2
   !*****************************************
      v1 = nobs + v0
      epsilon = (W - matmul(X, beta))
      Mb = matmul(M0, beta-beta0)
      s1 = dot_product(epsilon, epsilon) + &
           dot_product(beta-beta0, Mb) + s0

     sigma_e_sq = s1/(2.0 * random_gamma1(v1/2.0, First))

 !***************************************************
 ! sampling u
 !***************************************************
  do i = 1, NGroups
    ALLOCATE( YN(CTEnd(i)-CTStart(i)+1) )
    ALLOCATE( XN(size(YN), KVars) )

    YN = Y(CTStart(i):CTEnd(i))
    XN = X(CTStart(i):CTEnd(i),:)

    VarU = (sigma_e_sq*sigma_u_sq)/ &
         (sigma_e_sq + real(size(YN))*sigma_u_sq)
    MeanU = (sigma_u_sq * sum(YN-matmul(XN,beta)))/&
         (sigma_e_sq + real(size(YN))*sigma_u_sq)
    u(i) = dsqrt(VarU)*rnorm() + MeanU

    DEALLOCATE (YN, XN) 
   end do

END SUBROUTINE GibbsSampler



SUBROUTINE Geweke_Diagnostic(Draws, ndraws, npars, frac1, frac2, CD )
USE LA_PRECISION, ONLY:    wp => dp
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
USE LA_PRECISION, ONLY:    wp => dp
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
! numerical standard error (nse) is sqrt(SG/ndraws)
! relative numerical error (rne) is PostVar/SG
lags = floor( ndraws**(0.25_wp) ) + 1

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
  SG(i) = PostVars(i) + ( 2.0_wp*dot_product( (/(i,i=lags,1,-1)/), &
                       Covars(:,i)) )/ real(lags+1) 
 end do

 RNE = PostVars/SG 
 NSE = sqrt( SG/real(ndraws) ) 

Deallocate ( Covars )

END SUBROUTINE



SUBROUTINE GR_Diagnostic(draws, ndraws, npars, nchains, B, W, V, R)
USE LA_PRECISION, ONLY:    wp => dp
IMPLICIT none
INTEGER, INTENT(IN):: ndraws, npars, nchains
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN):: draws
REAL(wp), DIMENSION(npars), INTENT(OUT):: B, W, V, R
REAL(wp), DIMENSION(npars) :: TotalMeans, ChainMeans
INTEGER :: N, start, end
REAL(wp), DIMENSION(ndraws, npars):: Deviation_sq
INTEGER:: i, j


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
