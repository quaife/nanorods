PROGRAM FMM_TEST

use omp_lib

IMPLICIT NONE

INTEGER :: nsource, i
REAL *8, DIMENSION(100000) :: xs
REAL *8, DIMENSION(100000) :: ys
REAL *16 :: ostart, oend
COMPLEX *16, DIMENSION(100000) :: dip1
COMPLEX *16, DIMENSION(100000) :: dip2
COMPLEX *16, DIMENSION(100000) :: vel

nsource = 100000
xs = (/(i*0.1, i = 1,nsource)/)
ys = (/(i*0.1, i = 1,nsource)/)
dip1 = (/((1,1), i = 1,nsource)/)
dip2 = (/((1,1), i = 1,nsource)/)

!$ostart = omp_get_wtime()

call stokesDLPnew(nsource,xs,ys,dip1,dip2,vel)

!$oend = omp_get_wtime()

write(*,*) 'Walltime ellapsed', oend-ostart
END PROGRAM
