PROGRAM GetTexForm
! This program puts together estimation results
! from several files.
INTEGER, PARAMETER:: KVars=25
INTEGER, PARAMETER:: OtherVars=3
INTEGER, PARAMETER:: Length=2*KVars+OtherVars
REAL, DIMENSION(KVars):: Results00_beta, Results00_t 
REAL, DIMENSION(KVars):: Results01_beta, Results01_t 
REAL, DIMENSION(KVars):: Results02_beta, Results02_t 
REAL, DIMENSION(KVars):: Results10_beta, Results10_t 
REAL, DIMENSION(KVars):: Results11_beta, Results11_t 
REAL, DIMENSION(KVars):: Results12_beta, Results12_t
CHARACTER(25), DIMENSION(KVars)::VarName= (/'Start-of-period Urban TFR', &
           'Start-of-period Urban Q5 ', &
           'Inland Water             ', &
           'LECZ                     ', &
           'Dry subhumid             ', &
           'Semiarid                 ', &
           'Arid                     ', &
           'LECZ * Dry subhumid      ', &
           'LECZ * (Semiarid or arid)', &
           '100$<=$City Size$<$500   ', &
           '500$<=$CIty SIze$<$1,000 ', &
           'City Size $>=$1,000      ', &
           'Unknown-Unknown          ', &
           'Unknown-Proper           ', &
           'Unknown-Agglomeration    ', &
           'Unknown-Metro.Area       ', &
           'Proper-Unknown           ', &
           'Proper-Proper            ', &
           'Proper-Agglomeration     ', &
           'Agglomeration-Unknown    ', &
           'Agglomeration-Proper     ', &
           'Agglomeration-Metro.Area ', &
           'Metro.Area-Metro.Area    ' ,&
           'Others-Others            ', &
           'Constant                 ' /)

REAL:: Results00_rho, Results00_sigmau, Results00_sigma
REAL:: Results01_rho, Results01_sigmau, Results01_sigma
REAL:: Results02_rho, Results02_sigmau, Results02_sigma
REAL:: Results10_rho, Results10_sigmau, Results10_sigma
REAL:: Results11_rho, Results11_sigmau, Results11_sigma
REAL:: Results12_rho, Results12_sigmau, Results12_sigma

INTEGER:: i, ios
 
 
! import files with estimation results
open(unit=3, file='./EstResults00', ACTION='READ')
open(unit=4, file='./EstResults01', ACTION='READ')
open(unit=5, file='./EstResults02', ACTION='READ')
open(unit=16, file='./EstResults10', ACTION='READ')
open(unit=7, file='./EstResults11', ACTION='READ')
open(unit=8, file='./EstResults12', ACTION='READ')

! output 
open(unit=10, file='./AllEstResults_SpatReg1A.tex', ACTION='WRITE')
open(unit=11, file='./AllEstResults_SpatReg1A_101112.tex', ACTION='WRITE')
open(unit=12, file='./AllEstResults_SpatReg1A_000102.tex', ACTION='WRITE')

do i = 1, KVars 
 read(unit=3, Fmt='(2F12.3)', IOSTAT=ios) Results00_beta(i), Results00_t(i)
 read(unit=4, Fmt='(2F12.3)', IOSTAT=ios) Results01_beta(i), Results01_t(i) 
 read(unit=5, Fmt='(2F12.3)', IOSTAT=ios) Results02_beta(i), Results02_t(i)
 read(unit=16, Fmt='(2F12.3)', IOSTAT=ios) Results10_beta(i), Results10_t(i)
 read(unit=7, Fmt='(2F12.3)', IOSTAT=ios) Results11_beta(i), Results11_t(i)
 read(unit=8, Fmt='(2F12.3)', IOSTAT=ios) Results12_beta(i), Results12_t(i)
end do 
 read(unit=3, Fmt='(3F12.3)', IOSTAT=ios) Results00_rho, Results00_sigmau, &
                                 Results00_sigma
 read(unit=4, Fmt='(3F12.3)', IOSTAT=ios) Results01_rho, Results01_sigmau, &
                                 Results01_sigma
 read(unit=5, Fmt='(3F12.3)', IOSTAT=ios) Results02_rho, Results02_sigmau, &
                                 Results02_sigma
 read(unit=6, Fmt='(3F12.3)', IOSTAT=ios) Results10_rho, Results10_sigmau, &
                                 Results10_sigma
 read(unit=7, Fmt='(3F12.3)', IOSTAT=ios) Results11_rho, Results11_sigmau, &
                                 Results11_sigma
 read(unit=8, Fmt='(3F12.3)', IOSTAT=ios) Results12_rho, Results12_sigmau, &
                                 Results12_sigma

print *, "IO status is ",ios


! Write results in LaTex form (long table)
! not used. tables below are used.
write(10,*) '\begin{center}'
write(10,*) '\begin{longtable}{}'
write(10,*) '\caption{\label{ } }\\\hline'
write(10,*) '&\multicolumn{}{c}{}&\multicolumn{}{c}{} \\'
write(10,*) '%\cmidrule{2-3}\cmidrule{5-6}'
write(10,*) '& & & &\\\hline'
write(10,*) '\endfirsthead'
write(10,*) '\multicolumn{}{l}{\emph{... table \thetable{} continued}} \\\hline'
write(10,*) '&\multicolumn{}{c}{}&\multicolumn{}{c}{}'
write(10,*) '\\%\cmidrule{2-3}\cmidrule{4-5}'
write(10,*) '& & & &\\\hline'
write(10,*) '\endhead '
write(10,*) '\hline'
write(10,*) '\multicolumn{}{r}{\emph{Continued on next page...}}\\'
write(10,*) '\endfoot'
write(10,*) '\endlastfoot'
do i = 1,KVars 
 write(unit=10,fmt='(A, A, 6(F7.3, A))' ) VarName(i), '&',  &
       Results10_beta(i), '&', Results11_beta(i), '&',  & 
       Results12_beta(i), '&', Results00_beta(i), '&',  &
       Results01_beta(i), '&', Results02_beta(i), '\\' 
 write(unit=10,fmt='(A, 6(F7.2, A, A))' ) &
       "&(",Results10_t(i),')&', '(',Results11_t(i), ')&', &
       "(",Results12_t(i),')&', "(",Results00_t(i), ')&',  &
       "(",Results01_t(i),')&', "(",Results02_t(i), ')', '\\'
end do 
 write(unit=10, fmt='(A, 6(F7.3, A))') &
          '$/rho$&', Results10_rho, '&', Results11_rho, '&', &
                    Results12_rho, '&', Results00_rho, '&', &
                    Results01_rho, '&', Results02_rho, '\\'
 write(unit=10,fmt='(A, 6(F7.3, A))') &
          '$/sigma_u$&', Results10_sigmau, '&', Results11_sigmau, '&', &
                         Results12_sigmau, '&', Results00_sigmau, '&', &
                         Results01_sigmau, '&', Results02_sigmau, '\\'
 write(unit=10,fmt='(A, 6(F7.3, A))') &
          '$/sigma$&', Results10_sigma, '&', Results11_sigma, '&', &
                       Results12_sigma, '&', Results00_sigma, '&', &
                       Results01_sigma, '&', Results02_sigma, '\\'
write(10,*) '\hline'
write(10,*) '\end{longtable}'
write(10,*) '\end{center}'

! Write results in LaTex form (long table)
! not used. See below
write(11,*) '\begin{center}'
write(11,*) '\begin{longtable}{}'
write(11,*) '\caption{\label{ } }\\\hline'
write(11,*) '&\multicolumn{}{c}{}&\multicolumn{}{c}{} \\'
write(11,*) '%\cmidrule{2-3}\cmidrule{5-6}'
write(11,*) '& & & &\\\hline'
write(11,*) '\endfirsthead'
write(11,*) '\multicolumn{}{l}{\emph{... table \thetable{} continued}} \\\hline'
write(11,*) '&\multicolumn{}{c}{}&\multicolumn{}{c}{}'
write(11,*) '\\%\cmidrule{2-3}\cmidrule{4-5}'
write(11,*) '& & & &\\\hline'
write(11,*) '\endhead '
write(11,*) '\hline'
write(11,*) '\multicolumn{}{r}{\emph{Continued on next page...}}\\'
write(11,*) '\endfoot'
write(11,*) '\endlastfoot'
do i = 1,KVars
 write(unit=11,fmt='(A, A, 3(F7.3, A))' ) VarName(i), '&', &
       Results10_beta(i), '&', Results11_beta(i), '&', &
       Results12_beta(i), '\\'
 write(unit=11,fmt='(A, 3(F7.2, A, A))' ) &
       "&(",Results10_t(i),')&', "(",Results11_t(i), ')&', &
       "(",Results12_t(i), ')', '\\'
end do
 write(unit=11, fmt='(A, 3(F7.3, A))') &
          '$/rho$&', Results10_rho, '&', Results11_rho, '&', &
                    Results12_rho, '\\'
 write(unit=11,fmt='(A, 3(F7.3, A))') &
          '$/sigma_u$&', Results10_sigmau, '&', Results11_sigmau, '&', &
                         Results12_sigmau, '\\'
 write(unit=11,fmt='(A, 3(F7.3, A))') &
          '$/sigma$&', Results10_sigma, '&', Results11_sigma, '&', &
                       Results12_sigma, '\\'
write(11,*) '\hline'
write(11,*) '\end{longtable}'
write(11,*) '\end{center}'


! Write results in LaTex form (sidewaystable, longtable)
write(11,*) '\begin{center}'
write(11,*) '\begin{sidewaystable}'
write(11,*) '\caption{Panel data city growth regression models'
write(11,*) ' with spatially correlated errors, Distance-based'
write(11,*) ' spatial weights ($d_{ij}$ denotes the distance between'
write(11,*) ' cities $i$ and $j$), Spatial weights are standardized.'
write(11,*) ' Assume that spatial interaction among cities occurs'
write(11,*) ' only within country'
write(11,*) ' \label{} }'
write(11,*) '\begin{longtable}{lcccccccc}\hline'
write(11,*) '&\multicolumn{2}{c}{Moderate spatial decay:}&&'
write(11,*) '\multicolumn{2}{c}{Fast spatial decay:}&&'
write(11,*) '\multicolumn{2}{c}{Slow spatial decay:} \\'
write(11,*) '&\multicolumn{2}{c}{Spatial weights are}&&'
write(11,*) '\multicolumn{2}{c}{Spatial weights are}&&'
write(11,*) '\multicolumn{2}{c}{Spatial weights are} \\'
write(11,*) '&\multicolumn{2}{c}{based on $(1/d_{ij})$}&&'
write(11,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^2$}&&'
write(11,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^{0.5}$} \\'
write(11,*) '\cmidrule{2-3}\cmidrule{5-6}\cmidrule{8-9}'
write(11,*) '&Coefficient &(Z-statistic) &&'
write(11,*) 'Coefficient &(Z-statistic) &&'
write(11,*) 'Coefficient &(Z-statistic) \\\hline'
write(11,*) '\endfirsthead'
write(11,*) '\multicolumn{}{l}{\emph{... table \thetable{} continued}} \\\hline'
write(11,*) '&\multicolumn{2}{c}{Moderate spatial decay:}&&'
write(11,*) '\multicolumn{2}{c}{Fast spatial decay:}&&'
write(11,*) '\multicolumn{2}{c}{Slow spatial decay:} \\'
write(11,*) '&\multicolumn{2}{c}{Spatial weights are}&&'
write(11,*) '\multicolumn{2}{c}{Spatial weights are}&&'
write(11,*) '\multicolumn{2}{c}{Spatial weights are} \\'
write(11,*) '&\multicolumn{2}{c}{based on $(1/d_{ij})$}&&'
write(11,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^2$}&&'
write(11,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^{0.5}$} \\'
write(11,*) '\cmidrule{2-3}\cmidrule{5-6}\cmidrule{8-9}'
write(11,*) '&Coefficient &(Z-statistic) &&'
write(11,*) 'Coefficient &(Z-statistic) &&'
write(11,*) 'Coefficient &(Z-statistic) \\\hline'
write(11,*) '\endhead '
write(11,*) '\hline'
write(11,*) '\multicolumn{}{r}{\emph{Continued on next page...}}\\'
write(11,*) '\endfoot'
write(11,*) '\endlastfoot'

do i = 1,KVars
 write(unit=11,fmt='(A, A, 3(F7.3, A, A, F7.2, A))' ) VarName(i), '&', &
       Results10_beta(i), '&',  "(",Results10_t(i),')&&',  &
       Results11_beta(i), '&',  "(",Results11_t(i), ')&&', &
       Results12_beta(i), '&',  "(",Results12_t(i), ')\\'
end do
 write(unit=11, fmt='(A, 3(F7.3, A))') &
          '$\rho$&', Results10_rho, '&&&', Results11_rho, '&&&', &
                    Results12_rho, '&\\'
 write(unit=11,fmt='(A, 3(F7.3, A))') &
          '$\sigma_u$&', Results10_sigmau, '&&&', Results11_sigmau, '&&&', &
                         Results12_sigmau, '&\\'
 write(unit=11,fmt='(A, 3(F7.3, A))') &
          '$\sigma$&', Results10_sigma, '&&&', Results11_sigma, '&&&', &
                       Results12_sigma, '&\\'
write(11,*) '\hline'
write(11,*) '\end{longtable}'
write(11,*) '\end{sidewaystable}'
write(11,*) '\end{center}'



! Write results in LaTex form (long table)
write(12,*) '\begin{center}'
write(12,*) '\begin{sidewaystable}'
write(12,*) '\caption{Panel data city growth regression models'
write(12,*) ' with spatially correlated errors, Distance-based'
write(12,*) ' spatial weights ($d_{ij}$ denotes the distance between'
write(12,*) ' cities $i$ and $j$), Spatial weights are standardized.'
write(12,*) ' Assume that spatial interaction among cities occurs'
write(12,*) ' worldwide ($d_{ij}$ is calculated for all cities)'
write(12,*) ' \label{} }'

write(12,*) '\begin{longtable}{lcccccccc}\hline'
write(12,*) '&\multicolumn{2}{c}{Moderate spatial decay:}&&'
write(12,*) '\multicolumn{2}{c}{Fast spatial decay:}&&'
write(12,*) '\multicolumn{2}{c}{Slow spatial decay:} \\'
write(12,*) '&\multicolumn{2}{c}{Spatial weights are}&&'
write(12,*) '\multicolumn{2}{c}{Spatial weights are}&&'
write(12,*) '\multicolumn{2}{c}{Spatial weights are} \\'
write(12,*) '&\multicolumn{2}{c}{based on $(1/d_{ij})$}&&'
write(12,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^2$}&&'
write(12,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^{0.5}$} \\'
write(12,*) '\cmidrule{2-3}\cmidrule{5-6}\cmidrule{8-9}'
write(12,*) '&Coefficient &(Z-statistic) &&'
write(12,*) 'Coefficient &(Z-statistic) &&'
write(12,*) 'Coefficient &(Z-statistic) \\\hline'
write(12,*) '\endfirsthead'
write(12,*) '\multicolumn{}{l}{\emph{... table \thetable{} continued}} \\\hline'
write(12,*) '&\multicolumn{2}{c}{Moderate spatial decay:}&&'
write(12,*) '\multicolumn{2}{c}{Fast spatial decay:}&&'
write(12,*) '\multicolumn{2}{c}{Slow spatial decay:} \\'
write(12,*) '&\multicolumn{2}{c}{Spatial weights are}&&'
write(12,*) '\multicolumn{2}{c}{Spatial weights are}&&'
write(12,*) '\multicolumn{2}{c}{Spatial weights are} \\'
write(12,*) '&\multicolumn{2}{c}{based on $(1/d_{ij})$}&&'
write(12,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^2$}&&'
write(12,*) '\multicolumn{2}{c}{based on $(1/d_{ij})^{0.5}$} \\'
write(12,*) '\cmidrule{2-3}\cmidrule{5-6}\cmidrule{8-9}'
write(12,*) '&&Coefficient &(Z-statistic) &&'
write(12,*) 'Coefficient &(Z-statistic) &&'
write(12,*) 'Coefficient &(Z-statistic) \\\hline'
write(12,*) '\endhead '
write(12,*) '\hline'
write(12,*) '\multicolumn{}{r}{\emph{Continued on next page...}}\\'
write(12,*) '\endfoot'
write(12,*) '\endlastfoot'
do i = 1,KVars
 write(unit=12,fmt='(A, A, 3(F7.3, A, A, F7.2, A))' ) VarName(i), '&', &
       Results00_beta(i), '&',  "(",Results00_t(i),')&&',  &
       Results01_beta(i), '&',  "(",Results01_t(i), ')&&', &
       Results02_beta(i), '&',  "(",Results02_t(i), ')\\'
end do
 write(unit=12, fmt='(A, 3(F7.3, A))') &
          '$\rho$&', Results00_rho, '&&&', Results01_rho, '&&&', &
                    Results02_rho, '&\\'
 write(unit=12,fmt='(A, 3(F7.3, A))') &
          '$\sigma_u$&', Results00_sigmau, '&&&', Results01_sigmau, '&&&', &
                         Results02_sigmau, '&\\'
 write(unit=12,fmt='(A, 3(F7.3, A))') &
          '$\sigma$&', Results00_sigma, '&&&', Results01_sigma, '&&&', &
                       Results02_sigma, '&\\'
write(12,*) '\hline'
write(12,*) '\end{longtable}'
write(12,*) '\end{sidewaystable}'
write(12,*) '\end{center}'



 
END PROGRAM GetTexForm
