MODULE DataPaths
IMPLICIT none
!**************************************************************
!* This declares the data path, and input and output          * 
!* variables.                                                 * 
!************************************************************** 
CHARACTER(100)::      datapath = "/home/dokim/Sep09/spatReg1/InData/" 
CHARACTER(200)::      input
CHARACTER(200)::      latlong_input
END MODULE DataPaths


MODULE DataParameters
!**************************************************************
!* This declares (1) function(s) from LAPACK95 and            *
!* (2) variables in the program. Note that this program       *
!* uses LAPACK(Linear Algebra PACKage).                       *
!**************************************************************
 USE LA_PRECISION, ONLY:                  WP => DP
 USE F95_LAPACK, ONLY:                    LA_POSV, LA_SYSV
IMPLICIT none
 INTEGER, PARAMETER ::                     nobs = 10766
 INTEGER, PARAMETER ::                     N = 2408 
 INTEGER, PARAMETER ::                     TYear = 58 
 INTEGER, PARAMETER ::                     KVars = 25   
 INTEGER, DIMENSION(nobs) ::               CountryLocationID
 INTEGER, DIMENSION(nobs) ::               UrbanLocationID
 INTEGER, DIMENSION(nobs) ::               Year
 INTEGER, DIMENSION(nobs) ::               CityID
 REAL(wp), DIMENSION(nobs) ::              Y
 REAL(wp), DIMENSION(nobs,KVars) ::        X

 INTEGER, DIMENSION(N) :: CountryID 
 REAL(wp), DIMENSION(2, N) ::              coordinates


 ! Used in OLS
 REAL(wp), DIMENSION(KVars) ::             XTY
 REAL(wp), DIMENSION(KVars, KVars) ::      XTX
 REAL(wp), DIMENSION(KVars) ::             beta
 REAL(wp), DIMENSION(KVars, KVars) ::      var_beta
 REAL(wp), DIMENSION(nobs) ::              e       ! OLS residuals
 REAL(wp) ::                               sigma_sq

 ! Used for indices of unbalanced panel data 
 INTEGER, DIMENSION(TYear) ::              YStart, YEnd

  
 ! Used in GM 
 REAL(wp), DIMENSION(N, N)::                   W
 INTEGER, ALLOCATABLE, DIMENSION(:)::          CTObserved
 INTEGER ::                                    NYear
 REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::      XYear,WYear,WXYear
 REAL(wp), ALLOCATABLE, DIMENSION(:) ::        YYear, WYYear
 REAL(wp), ALLOCATABLE, DIMENSION(:) ::        eYear, WeYear, WWeYear
 REAL(wp), DIMENSION(nobs) ::                  YTilde
 REAL(wp), DIMENSION(nobs, KVars) ::           XTilde
 REAL(wp), DIMENSION(nobs) ::                  We, WWe

 REAL(wp), ALLOCATABLE, DIMENSION(:,:)::       WWYear
 REAL(wp), DIMENSION(N,N) ::                   DDYear
 REAL(wp), DIMENSION(N) ::                     DWWD
 INTEGER, DIMENSION(N) ::                      NObsCity 
 
 REAL(wp)::                                    trWTW
 REAL(wp)::                                    trWPW, trWQW
 
 REAL(wp) ::                                   trDWQWD, trDWWD, trDWPWD
 REAL(wp), DIMENSION(N) ::                     DWDDWD
 REAL(wp), DIMENSION(N,N) ::                   DWD


 ! moments
 REAL(wp), DIMENSION(nobs) ::                 Qe, QWe, QWWe
 REAL(wp), DIMENSION(nobs) ::                 Pe, PWe
 
 REAL(wp) ::                                  eQWe, eWQWe, eQe 
 REAL(wp) ::                                  eWWQWe, eWWQWWe
 REAL(wp) ::                                  eQWWe
 REAL(wp) ::                                  ePe, ePWe, eWPWe
 
 INTEGER, PARAMETER ::                        NumMoments = 4  
 INTEGER, PARAMETER ::                        NumVar = 4 
 INTEGER, PARAMETER ::                        NumParameters = 3 
 REAL(wp), DIMENSION(NumMoments,NumVar) ::    LHG
 REAL(wp), DIMENSION(NumMoments) ::           RHg
 REAL(wp), DIMENSION(NumParameters) ::        delta, delta0, ddelta
 REAL(wp),  DIMENSION(NumMoments) ::          Error, SE
 REAL(wp), DIMENSION(NumParameters) ::        g
 REAL(wp), DIMENSION(NumParameters,NumParameters) ::     H
 REAL(wp), DIMENSION(NumParameters,NumParameters) ::     Hinverse
 REAL(wp), DIMENSION(NumParameters,NumParameters) ::     identity_H
 REAL(wp), DIMENSION(KVars,KVars) ::          identity_beta


END MODULE DataParameters



PROGRAM GM_UNBALANCED_SPANEL_PARTONE 
USE DataPaths
USE DataParameters
IMPLICIT none
INTEGER ::                                 i, ios, j, t
INTEGER ::                                 iter
REAL(wp) ::                                SSE, newSSE
REAL ::                                    start, finish
REAL(wp)::                                 RowSum
!**************************************************************
!**************************************************************
!* This is a program for generalized moments (GM)             *
!* estimation of spatial error model in the case of           *
!* unbalanced and irregular-spaced panel data.                *
!**************************************************************
!**************************************************************

 input = adjustl(trim(datapath))//"GrowthSpatialData_M32.raw"
 latlong_input = adjustl(trim(datapath))//"LatLong.raw"

! import data
 OPEN (Unit = 8, File = input, ACTION="READ")
 DO i = 1, nobs
   READ (Unit = 8, Fmt = *, IOSTAT=ios) &
      Year(i), CountryLocationID(i), UrbanLocationID(i), CityID(i), &
                            Y(i), X(i,2:KVars)
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for Y, X data is ",ios

! import spatial Weights matrix
 OPEN (Unit = 8, File = latlong_input, ACTION="READ")
 DO i = 1, N
   READ (Unit = 8, Fmt = *, IOSTAT=ios) CountryID(i), coordinates(:,i)
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for weight matrix is ",ios
 print *, "number of city (latlong) is ", i-1

! Open a file recording estimation results
  OPEN (Unit = 10, File="./EstResults12", ACTION="WRITE")


!***************************************************************
! Establish indices for unbalanced panel                       *
!***************************************************************
 CALL UNBALANCEDPANEL_INDEX


!***************************************************************
! Run OLS to get OLS residuals                                 *
!***************************************************************
! Data include constant terms
 X(:,1) = 1.0_wp

 identity_beta = 0.0_wp
 do i = 1, KVars
   identity_beta(i,i) = 1.0_wp
 enddo


 CALL OLS

CALL CPU_TIME(start)

!***************************************************************
! Specify spatial weight matrix                                *
!***************************************************************
 CALL SWeightMatrix

!***************************************************************
! WY, WX, We, WWe, and traces of matrices                      *
!***************************************************************
 CALL PanelGM
 print *, "PanelGM done"

!***************************************************************
! Calculate sample moments                                     *
!***************************************************************
 CALL MOMENTS
 print *, "MOMENTS done"

!***************************************************************
! Find inital GM estimators of rho and sigma_sq                *
!***************************************************************
! initial guesses
  delta(1) = 0.1_wp       ! initial value of rho
  delta(2) = 1.0_wp 
  delta(3) = 1.0_wp 


! Create an identity matrix to be used in
  identity_H = 0.0_wp
    do i = 1, NumParameters
    identity_H(i,i) = 1.0_wp
  end do

! Newton optimization algorithm

  
  delta0 = delta

  call SSEfunc(SSE)


 iter = 30 

 DO i = 1, iter
   print *, 'iteration', i
   call direction

   delta0 = delta + ddelta
   call SSEfunc(newSSE)
   if (newSSE.lt.SSE) then
     delta = delta + ddelta
     print *, 'newSSE = ', newSSE
     SSE = newSSE
   else
     print *, 'adjusting step length'
     call step(SSE, newSSE)
     delta = delta + ddelta
     SSE = newSSE
   endif

 END DO

  print *, 'GM estimator of rho', delta(1)
  print *, 'GM estimator of sigma_sq', delta(2)
  print *, 'GM estimator of sigma_u_sq', delta(3)


!***************************************************************
!***************************************************************
! GLS estimation                                               *
!***************************************************************
 CALL UnbalancedGLS

 CALL CPU_TIME(finish)
 print *, "time", finish-start, "sec."


! Write estimation results
do i = 2, KVars
  write(10, "(2F12.3)" ) beta(i), beta(i)/sqrt(var_beta(i,i)) 
enddo
  write(10, "(2F12.3)" ) beta(1), beta(1)/sqrt(var_beta(1,1))
  write(10, "(3F12.3)" )   delta(1), sqrt(delta(3)), sqrt(sigma_sq)





END PROGRAM GM_UNBALANCED_SPANEL_PARTONE


SUBROUTINE OLS
USE DataParameters
IMPLICIT none
INTEGER::          i, j
!*************************************************************
! This solves for OLS estimator of beta, using the LU decomp.
! X'X and right-hand side vector X'Y. Make use of the columnf
! the X matrix and dot-products for efficiency.
!*************************************************************
  do i = 1, KVars
    do j = 1, KVars
      XTX(i,j) = dot_product(X(:,i),X(:,j))
    enddo
    XTY(i) = dot_product(X(:,i),Y)
  enddo

  call LA_POSV(XTX,XTY)
  beta = XTY
  print *, "   ------------------------------------"
  print *, "    OLS estimators for beta is" 
  print *, "   ------------------------------------"
  do i = 1, KVars
  print *, "   ", beta(i)
  end do

  e = Y - matmul(X, beta)
  sigma_sq = dot_product(e, e)/real(nobs)
  print *, "     OLS sigma squared is ", sigma_sq

RETURN
END SUBROUTINE OLS



SUBROUTINE UNBALANCEDPANEL_INDEX
 USE DataParameters
 IMPLICIT none
 INTEGER :: i, t
 INTEGER :: id
!*************************************************************
! This finds the row range for each year and the total       *
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
    if (i == nobs+1) EXIT
   END DO
   YEnd(t) = i-1
END DO

!print *, "The total number of Year is", t

END SUBROUTINE UNBALANCEDPANEL_INDEX


SUBROUTINE SWeightMatrix
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k
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
 r = 6378.0_wp - 21.0_wp * sin(coordinates(1,1))
! The following radius is what is used in the Stata add-on
! program globdist, and gives results that agree with the
! Stata program up to the second or third decial
! (for Malawi)

!r = 6365.0_wp

! Implement the Haversine great circle distance formula
!  Note: Can improve code by going down columns of dist
DO i = 1, N
 DO j = 1, N

 if ((CountryID(i) == CountryID(j)) .and. (i /= j) ) then
  W(i,j) = r * 2.0_wp *   &
  asin( min(1.0_wp, &
  sqrt( sin((coordinates(1,j)- coordinates(1,i))/2.0_wp)**2 &
   + cos(coordinates(1,i)) * cos(coordinates(1,j)) *    &
  sin( (coordinates(2,j)-coordinates(2,i))/2.0_wp )**2  &
  ) ))
 else
  W(i,j) = 0.0_wp
 end if

 END DO
END DO

print *, "max distance is", maxval(W(1,:))



!*************************************************************
! Specification: distance weight matrix                      *
!*************************************************************
DO i = 1, N
  DO j = 1, N
    if (W(i,j) /= 0.0_wp) then
      W(i,j) = sqrt(1.0_wp/W(i,j))
    end if
  END DO
END DO


END SUBROUTINE SWeightMatrix




SUBROUTINE PanelGM
USE DataParameters
IMPLICIT none
INTEGER :: i, j, t
REAL(wp):: RowSum
!***************************************************************
! WY, WX, We, WWe, and traces of several matrices              *
!***************************************************************

 trWTW = 0.0; DWWD = 0.0; trDWWD = 0.0; DWD = 0.0

YEAR_ROUTINE: DO t = 1, TYear
 print *, "Time", t 

 NYear = YEnd(t)-YStart(t)+1
 ALLOCATE( CTObserved(NYear) )
 ALLOCATE( WYear(NYear, NYear) )

 CTObserved = CityID ( YStart(t):YEnd(t) )
 WYear = W( CTObserved, CTObserved )

! row standardization for each year.
 DO i = 1, NYear
  RowSum = 0.0_wp
  RowSum = SUM( WYear(i,:) )
  DO j = 1, NYear
   if (WYear(i,j) /= 0) then
    WYear(i,j) = WYear(i,j)/RowSum
   end if
  END DO
 END DO


!*************************************************************
! Here, YTilde means WY and XTilde means WX                  *
!*************************************************************
! print *, "NYear is", NYear
 ALLOCATE( XYear(NYear, KVars) )
 ALLOCATE( YYear(NYear) )
 ALLOCATE( eYear(NYear) )
 ALLOCATE( WXYear(NYear, KVars) )
 ALLOCATE( WYYear(NYear) )
 ALLOCATE( WeYear(NYear) )
 ALLOCATE( WWeYear(NYear) )
 ALLOCATE( WWYear(NYear, NYear) )


 XYear = X( YStart(t):YEnd(t), :)
 YYear = Y( YStart(t):YEnd(t) )
 eYear = e( YStart(t):YEnd(t) )

 WXYear = MATMUL(WYear, XYear)
 WYYear = MATMUL(WYear, YYear)
 WeYear = MATMUL(WYear, eYear)
 WWeYear = MATMUL(WYear, WeYear)

 XTilde(YStart(t):YEnd(t), :) = WXYear
 YTilde(YStart(t):YEnd(t)) = WYYear
 We(YStart(t):YEnd(t)) = WeYear
 WWe(YStart(t):YEnd(t)) = WWeYear

 DEALLOCATE( YYear, XYear, eYear, WXYear, WYYear )
 DEALLOCATE( WeYear, WWeYear )

!**************************************************************
!**************************************************************
! Calculate the following traces                             **
! (1) Trace(W'PW) = trWPW = tr(DWWD * diag(1/Ti))            **
! (2) Trace(W'QW) = trWQW =  trWTW - trWPW                   **
! (3) Trace(W'PWDD') = trDWPWD = tr( DWDDWD*diag(1/Ti) )     **
! (4) Trace(W'QWDD') = trWWDD - trWPWDD = trDWWD - trDWPWD   **
!**************************************************************
!**************************************************************

!**************************************************************
! trWTW used in (2)                                           *
!**************************************************************
DO i = 1, NYear
 DO j = 1, NYear
   trWTW = trWTW + WYear(i,j)*WYear(i,j)
 END DO
END DO


!**************************************************************
! trDWWD means trD'W'WD = trWTW used in (4)                   *
! it is same as trWTW                                         *
!**************************************************************
 WWYear = 0.0_wp
 ! here, WWYear = WTW for a given year
 DO i = 1, NYear
  DO j = 1, NYear
   WWYear(i,j) = dot_product( WYear(:,i), WYear(:,j) )
  END DO
 END DO


!**************************************************************
! DWWD is D'WW'D used in (1)                                  *
!**************************************************************
! here, DDYear is DWWDYear
 DDYear = 0.0_wp
 DO i = 1, NYear
  DO j = 1, NYear
   DDYear(CTObserved(i), CTObserved(j)) = WWYear(i,j)
  END DO
 END DO

! Only Diagonal elements are necessary
DO i = 1, N
 DWWD(i) = DWWD(i) + DDYear(i,i)
END DO

!**************************************************************
! DWD                                                         *
!**************************************************************
! here, DDYear = DWDYear
DDYear = 0.0_wp
DO i = 1, NYear
 DO j = 1, NYear
  DDYear(CTObserved(i), CTObserved(j)) = WYear(i,j)
 END DO
END DO

 DWD = DWD + DDYear

DEALLOCATE( CTObserved, WYear, WWYear )

END DO YEAR_ROUTINE


trDWWD = trWTW

!**************************************************************
! (1) trWPW and (2) trWQW                                     *
!**************************************************************
! Total Number of Observations per city
DO i = 1, N
 NObsCity(i) = count( CityID == i )
END DO

trWPW = 0.0_wp
DO i = 1, N
 trWPW = trWPW + DWWD(i)*( 1.0_wp/REAL(NObsCity(i)) )
END DO
trWQW = trWTW - trWPW

!**************************************************************
! Diagonal elements of DWDDWD                                 *
!**************************************************************
DO i = 1, N
 DWDDWD(i) = dot_product( DWD(i,:), DWD(:,i) )
END DO

!*************************************************************
! (3) trDWPWD and (4) trDWQWD                                *
!*************************************************************
trDWPWD = 0.0_wp
DO i = 1, N
 trDWPWD = trDWPWD + DWDDWD(i) * ( 1.0_wp/REAL(NObsCity(i)) )
END DO
trDWQWD = trDWWD - trDWPWD

END SUBROUTINE PanelGM



SUBROUTINE MOMENTS
USE DataParameters
IMPLICIT none
INTEGER :: i, j
REAL(wp), ALLOCATABLE, DIMENSION(:):: temp1, temp2, temp3
REAL(wp), DIMENSION(N):: Meane, MeanWe, MeanWWe
INTEGER:: iwt
!**************************************************************
! This calculates elements of sample moment conditions by     *
! using the OLS residuals, e, and spatially weighted          *
! versions, We and WWe.                                       *
!**************************************************************
! Pe is average of e_it over time for each spatial unit.      *
! Qe is deviation from average of e_it over time.             *
!**************************************************************


 Do i = 1, N
  iwt = count(CityID == i)
  allocate ( temp1(iwt) )
  allocate ( temp2(iwt) )
  allocate ( temp3(iwt) )

  temp1 = PACK(e, CityID == i)
  temp2 = PACK(We, CityID == i)
  temp3 = PACK(WWe, CityID == i)

  Meane(i) = SUM(temp1)/REAL(NObsCity(i))
  MeanWe(i) = SUM(temp2)/REAL(NObsCity(i))
  MeanWWe(i) = SUM(temp3)/REAL(NObsCity(i))

 deallocate (temp1, temp2, temp3)

 END DO

 DO i = 1, nobs
  Pe(i) = Meane(CityID(i))
  PWe(i) = MeanWe(CityID(i))

  Qe(i) = e(i) - Meane(CityID(i))
  QWe(i) = We(i) - MeanWe(CityID(i))
  QWWe(i) = WWe(i) - MeanWWe(CityID(i))
 END DO

 eQWe = dot_product(e, QWe)
 eWQWe = dot_product(We, QWe)
 eQe = dot_product(e, Qe)

 eWWQWe = dot_product(WWe, QWe)
 eWWQWWe = dot_product(WWe, QWWe)
 eQWWe = dot_product(e, QWWe)

 ePe = dot_product(e, Pe)
 ePWe = dot_product(e, PWe)
 eWPWe = dot_product(We, PWe)

! Elements of samples moments
 LHG = 0.0_wp
 LHG(1,1) = 2.0_wp*eQWe
 LHG(1,2) = -eWQWe
 LHG(1,3) = nobs - N 
 
 LHG(2,1) = 2.0_wp*eWWQWe
 LHG(2,2) = -eWWQWWe
 LHG(2,3) = trWQW              !tr(W'QW) 
 LHG(2,4) = trDWQWD            !tr(W'QWDD') 
 
 LHG(3,1) = eQWWe + eWQWe
 LHG(3,2) = -eWWQWe
 
 LHG(4,1) = 2.0_wp*ePWe
 LHG(4,2) = -eWPWe
 LHG(4,3) = N
 LHG(4,4) = nobs

 RHg = 0.0_wp
 RHg(1) = eQe
 RHg(2) = eWQWe
 RHg(3) = eQWe
 RHg(4) = ePe

RETURN
END SUBROUTINE MOMENTS



SUBROUTINE SSEfunc(SSE)
USE DataParameters
IMPLICIT none
INTEGER ::  i, iter
REAL(wp), INTENT(OUT) ::     SSE


! Calculate SSE (Sum of Squared Error)
   do i = 1,  NumMoments
     Error(i) =LHG(i,1)*delta0(1) + LHG(i,2)*delta0(1)**2 +&
               LHG(i,3)*delta0(2) + LHG(i,4)*delta0(3)-RHg(i)
     SE(i) = Error(i)*Error(i)
   end do
    

    SSE = sum(SE)

RETURN
END SUBROUTINE SSEfunc



SUBROUTINE direction
USE DataParameters
IMPLICIT none
INTEGER ::  i

 ! gradient g
   g = 0.0_wp
   do i = 1, NumMoments
     g(1) = g(1) + 2.0_wp*Error(i)*(LHG(i,1)+2.0_wp*delta0(1)*LHG(i,2))
     g(2) = g(2) + 2.0_wp*Error(i)*LHG(i,3)
     g(3) = g(3) + 2.0_wp*Error(i)*LHG(i,4)
   end do

 ! hessian H
   H = 0.0_wp
   do i = 1, NumMoments
     H(1,1) = H(1,1) + 2.0_wp*(LHG(i,1)+2.0_wp*delta0(1)*LHG(i,2))**2.0 &
                 + 4.0_wp*Error(i)*LHG(i,2)
     H(1,2) = H(1,2) + 2.0_wp*(LHG(i,1)+2.0_wp*delta0(1)*LHG(i,2))*LHG(i,3)
     H(1,3) = H(1,3) + 2.0_wp*(LHG(i,1)+2.0_wp*delta0(1)*LHG(i,2))*LHG(i,4)
     H(2,2) = H(2,2) + 2.0_wp*LHG(i,3)*LHG(i,3)
     H(2,3) = H(2,3) + 2.0_wp*LHG(i,4)*LHG(i,3)
     H(3,3) = H(3,3) + 2.0_wp*LHG(i,4)*LHG(i,4)
   end do
     H(2,1) = H(1,2)
     H(3,1) = H(1,3)
     H(3,2) = H(2,3)
  Hinverse = identity_H
  call la_sysv(H, Hinverse)

! Search direction
  ddelta =-matmul(Hinverse, g)


RETURN
END SUBROUTINE direction



SUBROUTINE step(SSE, newSSE)
USE  DataParameters
IMPLICIT none
REAL(wp), INTENT(IN)    :: SSE
REAL(wp), INTENT(INOUT) :: newSSE

DO WHILE(newSSE < SSE-1d-8)
   print *, 'SSE and last newSSE', SSE, newSSE

   ddelta = 0.8_wp*ddelta
   delta0 = delta + ddelta
   call SSEfunc(newSSE)
END DO

RETURN
END SUBROUTINE step


SUBROUTINE UnbalancedGLS
USE DataParameters
IMPLICIT none
INTEGER :: i, j, t
REAL(wp), ALLOCATABLE, DIMENSION(:) ::   YTemp
REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: XTemp 
REAL(wp), DIMENSION(N) ::                Theta
REAL(wp), DIMENSION(N) ::                MeanY
REAL(wp), DIMENSION(N, KVars) ::         MeanX
REAL(wp), DIMENSION(KVars,KVars)::       XX
!********************************************************************
! Generalized Least Square (GLS) estimation
!********************************************************************

! first transformation
! XTilde in the right-hand size is WX.
  XTilde = X - delta(1) * XTilde
  YTilde = Y - delta(1) * YTilde

! second transformation
 DO i = 1, N
  Theta(i) = 1.0_wp - SQRT(delta(2)/(COUNT(CityID==i)*delta(3) + delta(2)))

  ALLOCATE( YTemp(COUNT(CityID==i)) )
  ALLOCATE( XTemp(COUNT(CityID==i), KVars) ) 

   YTemp = PACK(YTilde, CityID == i)
   MeanY(i) = SUM(YTemp)/REAL(COUNT(CityID == i))
 
   DO j = 1, KVars
    XTemp(:,j) = PACK(XTilde(:,j), CityID == i) 
   END DO
   MeanX(i,:) = SUM(XTemp,1)/REAL(COUNT(CityID == i))

  DEALLOCATE( YTemp, XTemp )
 END DO

 DO i = 1, nobs
  XTilde(i,:) = XTilde(i,:) - theta(CityID(i))*MeanX(CityID(i),:)
  YTilde(i) = YTilde(i) - theta(CityID(i))*MeanY(CityID(i))
 END DO

 DO i = 1, KVars
   DO j = 1, KVars
     XTX(i,j) = dot_product(XTilde(:,i),XTilde(:,j))
   END DO 
   XTY(i) = dot_product(XTilde(:,i),YTilde)
 END DO

  XX = XTX

  call la_posv(XTX,XTY)
  beta = XTY
  print *, "Final beta is ", beta

  e = YTilde - matmul(XTilde,beta)
  sigma_sq = dot_product(e,e)/real(nobs)
  print *, "Final sigma squared is ", sigma_sq

  ! Invert XTX
  call la_posv(XX,identity_beta)

  ! Multiply the inverse by sigma-squared to obtain
  !  the variance matrix for beta
  var_beta = sigma_sq*identity_beta

  do i = 1, KVars
   print *, beta(i), sqrt(var_beta(i,i)), beta(i)/sqrt(var_beta(i,i))
  enddo


RETURN
END SUBROUTINE UnbalancedGLS

