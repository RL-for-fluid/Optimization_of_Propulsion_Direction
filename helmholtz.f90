!*****************************************************************************!
  SUBROUTINE helmholtz(ua,ut,fu,u,n)
!*****************************************************************************!
#ifdef _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,planf,plani,nxnyi,dx2i,dy2i,laplacei,alpha,dti,rei
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: ua(nx,ny) !速度^ast
REAL(nk),INTENT(IN) :: ut(nx,ny)  !速度^tilde
REAL(nk),INTENT(IN) :: fu(nx,ny)  !ibm用の外力
REAL(nk),INTENT(IN) :: u(nx,ny)   !速度
INTEGER,INTENT(IN) :: n           !ルンゲクッタのstep

REAL(nk) :: ci                              !係数
REAL(nk) :: rh(nx,ny)                       !右辺*正規化
COMPLEX(nk) :: rhf(0:nx/2,0:ny-1)           !右辺*正規化のフーリエ変換
COMPLEX(nk) :: uaf(0:nx/2,0:ny-1)           !速度^astのフーリエ変換
INTEGER :: i,j,im,jm,ip,jp                  !格子
INTEGER :: l,m                              !波数

!==check用==!
REAL(nk) :: error(nx,ny)
!=============================================================================!

!===右辺*正規化===!
!#ifdef _TIMER_
!  CALL timer_start(551)
!#endif

ci = 1/(rei*alpha(n))

!$OMP PARALLEL DO PRIVATE(im,ip,jm,jp) COLLAPSE(2)
DO j = 1, ny 
DO i = 1, nx
  im = MOD(i-2+nx,nx)+1 ; ip = MOD(i,nx)+1
  jm = MOD(j-2+ny,ny)+1 ; jp = MOD(j,ny)+1
  rh(i,j) = -(ut(i,j)*dti+fu(i,j))*ci &
            +(u(ip,j)+u(im,j)-2*u(i,j))*dx2i+(u(i,jp)+u(i,jm)-2*u(i,j))*dy2i
  rh(i,j)=rh(i,j)*nxnyi
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(551)
!#endif

!===右辺*正規化のFourier変換===!
!#ifdef _TIMER_
!  CALL timer_start(552)
!#endif
  IF(nk==8) CALL dfftw_execute_dft_r2c(planf,rh,rhf)
  IF(nk==4) CALL sfftw_execute_dft_r2c(planf,rh,rhf)
!#ifdef _TIMER_
!  CALL timer_end(552)
!#endif

!===漸化式===!
!#ifdef _TIMER_
!  CALL timer_start(553)
!#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 0,ny-1
DO l = 0,nx/2
  uaf(l,m) = laplacei(l,m,n)*rhf(l,m)
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(553)
!#endif

!===速度^astの逆Fourier変換===!
!#ifdef _TIMER_
!  CALL timer_start(554)
!#endif
  IF(nk==8) CALL dfftw_execute_dft_c2r(plani,uaf,ua)
  IF(nk==4) CALL sfftw_execute_dft_c2r(plani,uaf,ua)
!#ifdef _TIMER_
!  CALL timer_end(554)
!#endif

!===check===!
!#ifdef _TIMER_
!  CALL timer_start(555)
!#endif
!  CALL check(error,ua,rh,n)
!  WRITE(*,*) n, 'step, max error(helmholtz)',maxval(abs(error))
!#ifdef _TIMER_
!  CALL timer_end(555)
!#endif


END SUBROUTINE helmholtz
