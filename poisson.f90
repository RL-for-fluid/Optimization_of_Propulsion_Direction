!*****************************************************************************!
  SUBROUTINE poisson(phi,ua,va,n)
!*****************************************************************************!
#ifdef _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,planf,plani,nxnyi,dxi,dyi,laplacei,alpha,dt
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: phi(nx,ny)         !擬似圧力
REAL(nk),INTENT(IN) :: ua(nx,ny),va(nx,ny) !速度^ast
INTEGER,INTENT(IN) :: n                    !ルンゲクッタのstep

REAL(nk) :: ci                             !係数
REAL(nk) :: rh(nx,ny)                      !右辺＊正規化
COMPLEX(nk) :: rhf(0:nx/2,0:ny-1)          !右辺＊正規化のフーリエ変換
COMPLEX(nk) :: phif(0:nx/2,0:ny-1)         !擬似圧力のフーリエ変換
INTEGER :: i,j,im,jm                       !格子
INTEGER :: l,m                             !波数

!==check用==!
REAL(nk) :: error(nx,ny)
!=============================================================================!

!===poissonの右辺===!
!#ifdef _TIMER_
!  CALL timer_start(561)
!#endif

ci = 1/(2*alpha(n)*dt)

!$OMP PARALLEL DO PRIVATE(im,jm) COLLAPSE(2)
DO j = 1, ny 
DO i = 1, nx
  im = i-1 ; IF(i==1) im = nx
  jm = j-1 ; IF(j==1) jm = ny
  rh(i,j) = (  (ua(i,j)-ua(im,j))*dxi &
              +(va(i,j)-va(i,jm))*dyi ) *nxnyi*ci
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(561)
!#endif

!===右辺＊正規化のFourier変換===!
!#ifdef _TIMER_
!  CALL timer_start(562)
!#endif
  IF(nk==8) CALL dfftw_execute_dft_r2c(planf,rh,rhf)
  IF(nk==4) CALL sfftw_execute_dft_r2c(planf,rh,rhf)
!#ifdef _TIMER_
!  CALL timer_end(562)
!#endif

!===漸化式===!
!#ifdef _TIMER_
!  CALL timer_start(563)
!#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 0,ny-1
DO l = 0,nx/2
  phif(l,m) = laplacei(l,m,0)*rhf(l,m)
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(563)
!#endif

!===擬似圧力の逆Fourier変換===!
!#ifdef _TIMER_
!  CALL timer_start(564)
!#endif
  IF(nk==8) CALL dfftw_execute_dft_c2r(plani,phif,phi)
  IF(nk==4) CALL sfftw_execute_dft_c2r(plani,phif,phi)
!#ifdef _TIMER_
!  CALL timer_end(564)
!#endif

!===check===!
!#ifdef _TIMER_
!  CALL timer_start(565)
!#endif
!  CALL check(error,phi,rh,0)
!  WRITE(*,*) n, 'step, max error(poisson)',maxval(abs(error))
!#ifdef _TIMER_
!  CALL timer_end(565)
!#endif


END SUBROUTINE poisson
