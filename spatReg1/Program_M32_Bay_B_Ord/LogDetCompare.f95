MODULE Precision
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
 INTEGER, PARAMETER ::                     nobs = 10766 
 INTEGER, PARAMETER ::                     N = 2408 
 INTEGER, PARAMETER ::                     TYear = 58 
 INTEGER, PARAMETER ::                     KVars = 25 
 INTEGER, DIMENSION(nobs) ::               Year
 INTEGER, DIMENSION(nobs) ::               CountryLocationID
 INTEGER, DIMENSION(nobs) ::               UrbanLocationID
 INTEGER, DIMENSION(nobs) ::               CityID
 REAL(wp), DIMENSION(nobs) ::              Y
 REAL(wp), DIMENSION(nobs,KVars) ::        X
 INTEGER, DIMENSION(nobs)::                CountryID
 REAL(wp), DIMENSION(2, N) ::              coordinates

 INTEGER, DIMENSION(TYear) ::              YStart, YEnd
 REAL(wp), DIMENSION(N, N)::               W    ! weight matrix
 INTEGER, ALLOCATABLE, DIMENSION(:)::      CTObserved
 REAL(wp), ALLOCATABLE, DIMENSION(:)::      RowSum
 REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  CYear, WYear

!Used in Ord's method for log-determinant
 REAL(wp), ALLOCATABLE, DIMENSION(:) ::    eigsYear
 REAL(wp), DIMENSION(nobs) ::              eigs

! Used in Chebyshev approx. methof for log_determinant
 REAL(wp), PARAMETER::                     pi = 3.1415927_wp
 INTEGER, PARAMETER::                      nposs = 3
 REAL(wp)::                                trWW
 REAL(wp), DIMENSION(nposs)::              tvec
 REAL(wp), DIMENSION(nposs)::              xk
 REAL(wp), DIMENSION(nposs, nposs)::       vcon

 REAL(wp)::                                rhotemp, logdet
 INTEGER, PARAMETER::                      nseq = 2001
 REAL(wp), DIMENSION(nseq,5)::             logdetseq

! Used in LU decomposition method
REAL(wp), DIMENSION(nobs, nobs)::          AllW, IRWM 

END MODULE DataParameters



PROGRAM Logdet_test 
USE Precision
USE DataParameters
IMPLICIT none
CHARACTER(100)::              DataPath = "./../InData/"
CHARACTER(100)::              CityData, LatLonData
INTEGER ::                    i, j, t, id, ios
REAL ::                       start, finish
CHARACTER(5)::                skip
!**************************************************************
! This code calculates the log-determinants for the values
! of rho ranging from -1 to 1 with the following methods.
!  (1) Ord's method
!  (2) Chebyshev approximation method
!  (3) LU decomposition method
!**************************************************************


!*****************************************************************
! read data                                                      *
!*****************************************************************
! import data
 CityData = adjustl(trim(DataPath))//"GrowthSpatialData_M32.raw"

 OPEN (Unit = 3, File = CityData, ACTION="READ")
 DO i = 1, nobs
   READ (Unit = 3, Fmt = *, IOSTAT=ios)  Year(i), CountryLocationID(i), &
                            UrbanLocationID(i), CityID(i), &
                            Y(i), X(i,2:KVars)
 END DO
 CLOSE (Unit = 3)
 print *, "IO status for Y, X data is ",ios

! Here, only CityID is used.

! import latitude and logitude data
 LatLonData = adjustl(trim(datapath))//"LatLong.raw"

 OPEN (Unit = 3, File = LatLonData, ACTION="READ")
 DO i = 1, N
   READ (Unit = 3, Fmt = *, IOSTAT=ios) CountryID(i), coordinates(:,i)
 END DO
 CLOSE (Unit = 3)
 print *, "IO status for weight matrix is ",ios

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

! ranges of rho to caclulate log-determinants.
  logdetseq(1,1) = -1.0_wp
 do i = 2, nseq
  logdetseq(i,1) = logdetseq(i-1,1) + 0.001_wp
 end do

!***************************************************************
! Ord's method                                                 *
!***************************************************************
CALL CPU_TIME(start)

 CALL Prep_LogDet_Ord
 print *, "max eigenvalues of W is ", maxval(eigs)
 print *, "min eigenvalues of W is ", minval(eigs)

 do i = 1, nseq
   rhotemp = logdetseq(i,1)
   CALL LogDet_Ord
   logdetseq(i,2) = logdet
 end do

CALL CPU_TIME(finish)
  print *, "time for Ord's method ", finish-start, "sec."

!***************************************************************
! Chebyshev approximation method                               *
!***************************************************************
CALL CPU_TIME(start)

 CALL Prep_LogDet_Chebyshev
 
 do i = 1, nseq
   rhotemp = logdetseq(i,1)
   CALL LogDet_Chebyshev
   logdetseq(i,3) = logdet
 end do

CALL CPU_TIME(finish)
  print *, "time for Chebyshev ", finish-start, "sec."


!***************************************************************
! LU decomposition method                                      *
!***************************************************************
skip = "no"
if (skip == "no") then
CALL CPU_TIME(start)

 CALL  Prep_LogDet_LU
  
 do i = 1, nseq
   rhotemp = logdetseq(i,1) 
   CALL LogDet_LU
   logdetseq(i,4) = logdet
 end do

CALL CPU_TIME(finish)
  print *, "time for LU ", finish-start, "sec."
end if

!***************************************************************
! LU decomposition method                                      *
!***************************************************************
if (skip == "no") then
CALL CPU_TIME(start)

 do i = 1, nseq
   rhotemp = logdetseq(i,1)
   CALL LogDet_LU_alt
   logdetseq(i,5) = logdet
 end do

CALL CPU_TIME(finish)
  print *, "time for LU ", finish-start, "sec."
end if


! results
open(3, file = './LogDet_ByRho.raw', action = 'write')
 do i = 1, nseq
  write(3,*) logdetseq(i,:)
 end do

END PROGRAM Logdet_test



SUBROUTINE Prep_LogDet_Ord
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
!**************************************************************
! Ord's method uses eigenvalues of the weight matrix.         *
! As a preparation, calculate the eigenvalues.                *
!**************************************************************

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

!************************************************************
!Notes:
! (1) LA_SYEV is used for a symmetric matrix.
! (2) Input matrix WYear is corrupted after LA_SYEV.
!************************************************************
CALL LA_SYEV(WYear, eigsYear)
eigs(YStart(t):YEnd(t)) = eigsYear

DEALLOCATE (CTObserved, CYear, WYear, RowSum, eigsYear)
END DO


END SUBROUTINE Prep_LogDet_Ord


SUBROUTINE LogDet_Ord
use Precision
USE DataParameters
IMPLICIT none
INTEGER ::           i, j
REAL(wp), DIMENSION(nobs) :: detp

!*****************************************************************
! Get the log of the Jacobian for this value of rho, which is
! sum of the (1 - rho*eigenvalue) for all eigenvalues.
!*****************************************************************
  detp = 1.0_wp - rhotemp*eigs

  logdet = sum (log(detp) )

RETURN
END SUBROUTINE LogDet_Ord



SUBROUTINE Prep_LogDet_Chebyshev
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
INTEGER, DIMENSION(nposs):: seq
!****************************************************************
! Chebyshev approximation method uses Chebyshev coefficents and *
! polynomials:                                                  *
! c1*tr(I) + c2*tr(W) + c3*(2*tr(WW)-tr(I)) - 0.2*c1*tr(I)      *
!****************************************************************

!****************************************************************
! Chebyshev polynomials                                         *
!****************************************************************
! calculate trace(WW)
trWW = 0.0_wp

DO t = 1, TYear

 ALLOCATE( CTObserved(YEnd(t)-YStart(t)+1) )
 ALLOCATE( CYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( WYear(size(CTObserved), size(CTObserved)) )
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


! trWW
DO i = 1, size(CTObserved)
 DO j = 1, size(CTObserved)
  trWW = trWW + WYear(i,j)**2.0_wp
 END DO
END DO

DEALLOCATE (CTObserved, CYear, WYear, RowSum)
END DO

tvec(1) = 0.5_wp*real(nobs)
tvec(2) = 0.0_wp
tvec(3) = 2.0_wp*trWW - real(nobs)

!****************************************************************
! Chebyshev coefficients: calculate the terms in c1, c2 and c3  *
! which does not involve rho, as a preparation.                 *
!****************************************************************
seq = (/1,2,3/)

DO j = 1, 3
 vcon(:,j) = COS ( pi*(j-1)*(real(seq)-0.5_wp)/ real(nposs) )
END DO

xk = COS( pi*(seq - 0.5_wp)/ real(nposs) )

RETURN
END SUBROUTINE Prep_LogDet_Chebyshev



SUBROUTINE LogDet_Chebyshev
use Precision
USE DataParameters
IMPLICIT none
INTEGER ::           i, j
REAL(wp), DIMENSION(nposs):: cposs
!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************

do i = 1, nposs
  cposs(i) =(2.0_wp/real(nposs))* &
                       SUM ( log(1.0_wp - rhotemp*xk) * vcon(:,i) )
end do

logdet = dot_product(cposs, tvec)

END SUBROUTINE LogDet_Chebyshev



SUBROUTINE Prep_LogDet_LU
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
!*************************************************************
! LU decomposition method uses the LU decomposition of       *
! (I - rho*W) matrix. Log-determinant is the sum of the      *
! diagonal elements of upper triangular matrix U.            *
!*************************************************************
! As a preparation, it creates the weight matrix for all 
! years, that is, diag(W_1, ..., W_T) whose dimension is
! n by n.
!*************************************************************

DO t = 1, TYear

 ALLOCATE( CTObserved(YEnd(t)-YStart(t)+1) )
 ALLOCATE( CYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( WYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( RowSum(size(CTObserved)) )

 CTObserved = CityID ( YStart(t):YEnd(t) )
 CYear = W( CTObserved, CTObserved )

!*************************************************************
! row standardization for each year.                         *
!*************************************************************
! Invserse of row sum calculated by column sum since CYear
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


AllW(YStart(t):YEnd(t), YStart(t):YEnd(t)) = WYear

DEALLOCATE (CTObserved, CYear, WYear, RowSum)
END DO


END SUBROUTINE Prep_LogDet_LU



SUBROUTINE LogDet_LU
use Precision
USE DataParameters
IMPLICIT none
INTEGER ::           i, j
INTEGER, DIMENSION(nobs) ::  ipivot
INTEGER ::                   infor
!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************
  do i = 1, nobs
    do j = 1, nobs
     if (i == j) then
      IRWM(i,j) = 1.0_wp - rhotemp*AllW(i,j)
     else
      IRWM(i,j) = - rhotemp*AllW(i,j)
     end if
    end do
  end do

  call LA_GETRF(nobs, nobs, IRWM , nobs, ipivot, infor)
   logdet = 0.0
   do i = 1, nobs
     logdet = logdet + log(abs(IRWM(i,i)))
   enddo

END SUBROUTINE LogDet_LU



SUBROUTINE LogDet_LU_alt
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j, k, t
INTEGER, DIMENSION(nobs) ::  ipivot
INTEGER ::                   infor
!*************************************************************
! LU decomposition method uses the LU decomposition of       *
! (I - rho*W) matrix. Log-determinant is the sum of the      *
! diagonal elements of upper triangular matrix U.            *
!*************************************************************
! Instead of creating diag(W_1, ..., W_T), the whole weight
! matrix, it uses the property that determinant of a block
! diagonal matrix is equal to the product of determinants
! of the blocks.
!*************************************************************

logdet = 0.0_wp

DO t = 1, TYear

 ALLOCATE( CTObserved(YEnd(t)-YStart(t)+1) )
 ALLOCATE( CYear(size(CTObserved), size(CTObserved)) )
 ALLOCATE( WYear(size(CTObserved), size(CTObserved)) )
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

!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************
  do i = 1, size(CTObserved) 
    do j = 1, size(CTObserved) 
     if (i == j) then
      WYear(i,j) = 1.0_wp - rhotemp*WYear(i,j)
     else
      WYear(i,j) = - rhotemp*WYear(i,j)
     end if
    end do
  end do

  call LA_GETRF(size(CTObserved), size(CTObserved), WYear, &
                size(CTObserved), ipivot, infor)
 
  do i = 1, size(CTObserved) 
     logdet = logdet + log(abs(WYear(i,i)))
  enddo

DEALLOCATE (CTObserved, CYear, WYear, RowSum)
END DO


END SUBROUTINE LogDet_LU_alt




SUBROUTINE SWeightMatrix
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::                  i, j
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

RETURN
END SUBROUTINE SWeightMatrix

