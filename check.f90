!*****************************************************************************!
  SUBROUTINE check(error,p,rh,n)
!*****************************************************************************!
USE parameter, ONLY:nk,nx,ny,dx2i,dy2i,rei,dt,alpha
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: error(nx,ny)      !局所誤差
REAL(nk),INTENT(IN) :: p(nx,ny),rh(nx,ny) !解と右辺
INTEGER,INTENT(IN) :: n                   !ルンゲクッタのstep
REAL(nk) :: ci                            !係数
INTEGER :: i,j,im,ip,jm,jp                !格子
!=============================================================================!
!==poisson方程式用==!
IF(n==0)THEN 
!$OMP PARALLEL DO PRIVATE(im,ip,jm,jp) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  im = MOD(i-2+nx,nx)+1 ; ip = MOD(i,nx)+1
  jm = MOD(j-2+ny,ny)+1 ; jp = MOD(j,ny)+1
  error(i,j) = &
   (p(im,j)-2*p(i,j)+p(ip,j))* dx2i &
  +(p(i,jm)-2*p(i,j)+p(i,jp))* dy2i &
  -rh(i,j)*nx*ny
ENDDO
ENDDO
!$OMP END PARALLEL DO
ENDIF

!==Helmholtz方程式用==!
IF(n>0)THEN 
ci = 1/(alpha(n)*rei*dt)
!$OMP PARALLEL DO PRIVATE(im,ip,jm,jp) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  im = MOD(i-2+nx,nx)+1 ; ip = MOD(i,nx)+1
  jm = MOD(j-2+ny,ny)+1 ; jp = MOD(j,ny)+1
  error(i,j) = &
   (p(im,j)-2*p(i,j)+p(ip,j))* dx2i &
  +(p(i,jm)-2*p(i,j)+p(i,jp))* dy2i &
  -p(i,j)*ci                        &
  -rh(i,j)*nx*ny
ENDDO
ENDDO
!$OMP END PARALLEL DO
ENDIF

END SUBROUTINE check
