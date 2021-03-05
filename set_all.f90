!*****************************************************************************!
  SUBROUTINE set_all
!*****************************************************************************!
USE parameter, ONLY:nk,nx,ny,dx,dy,lxq,lyq,planf,plani,fftp,ffts &
                  ,laplacei,piq,alpha,gamma,zeta,reiq,dtq
use omp_lib
IMPLICIT NONE
INCLUDE './FFTW/include/fftw3.f'
INTEGER :: ierr,mcore=4,FFTW_METHOD = FFTW_PATIENT !fftw
REAL(16) :: alphaq(3)                               !coefficients of RK method
COMPLEX(16) :: laplace(0:nx/2,0:ny-1,0:3)           !for poisson and Helmholtz eq
COMPLEX(16) :: laplace_i(0:nx/2,0:ny-1,0:3)         !for poisson and Helmholtz eq
INTEGER :: l,m                                      !wave number
INTEGER :: n                                        !step of RK
!=============================================================================!
!==coefficients of RK method==!
alphaq(1)= 4.Q0/15.Q0 ; alphaq(2)=   1.Q0/15.Q0 ; alphaq(3)=  1.Q0/6.Q0
alpha(1) = 4.Q0/15.Q0 ; alpha(2) =   1.Q0/15.Q0 ; alpha(3) =  1.Q0/6.Q0
gamma(1) = 8.Q0/15.Q0 ; gamma(2) =   5.Q0/12.Q0 ; gamma(3) =  3.Q0/4.Q0
zeta(1)  = 0.Q0       ; zeta(2)  = -17.Q0/60.Q0 ; zeta(3)  = -5.Q0/12.Q0

!==FFTW==!
IF(nk==8) THEN
!  CALL dfftw_init_threads(ierr)
!  CALL dfftw_plan_with_nthreads(omp_get_max_threads()) 
  CALL dfftw_plan_many_dft_r2c( planf,2,(/nx,ny/),1,fftp,(/nx,ny/) &
       ,1,1,ffts,(/(nx/2+1),ny/),1,1,FFTW_METHOD )
  CALL dfftw_plan_many_dft_c2r( plani,2,(/nx,ny/),1,ffts,(/(nx/2+1),ny/) &
       ,1,1,fftp,(/nx,ny/),1,1,FFTW_METHOD )
ENDIF
IF(nk==4) THEN
!  CALL sfftw_init_threads(ierr)
!  CALL sfftw_plan_with_nthreads(omp_get_max_threads())
  CALL sfftw_plan_many_dft_r2c( planf,2,(/nx,ny/),1,fftp,(/nx,ny/) &
       ,1,1,ffts,(/(nx/2+1),ny/),1,1,FFTW_METHOD )
  CALL sfftw_plan_many_dft_c2r( plani,2,(/nx,ny/),1,ffts,(/(nx/2+1),ny/) &
       ,1,1,fftp,(/nx,ny/),1,1,FFTW_METHOD )
ENDIF

!==for poisson and Helmholtz Eq==!
!for poisson
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 0,ny-1
DO l = 0,nx/2
  laplace(l,m,0) = &
           (2*COS( 2*piq*l/nx )-2)*(nx/lxq)**2 &
          +(2*COS( 2*piq*m/ny )-2)*(ny/lyq)**2 
ENDDO
ENDDO
!$OMP END PARALLEL DO
!for Helmholtz
!$OMP PARALLEL DO COLLAPSE(3)
DO n = 1,3
DO m = 0,ny-1
DO l = 0,nx/2
  laplace(l,m,n) = &
           (2*COS( 2*piq*l/nx )-2)*(nx/lxq)**2 &
          +(2*COS( 2*piq*m/ny )-2)*(ny/lyq)**2 &
          -1/(alphaq(n)*reiq*dtq)
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(3)
DO n = 0,3
DO m = 0,ny-1
DO l = 0,nx/2
  IF(l==0.AND.m==0.AND.n==0) THEN !only for poisson
    laplace_i(l,m,n)=0
  ELSE
    laplace_i(l,m,n)= 1/laplace(l,m,n) 
  ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(3)
DO n = 0,3
DO m = 0,ny-1
DO l = 0,nx/2
  laplacei(l,m,n)= laplace_i(l,m,n) 
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE set_all
