!*****************************************************************************!
 SUBROUTINE ibm(delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,pos,angle,vel,avel)
!*****************************************************************************!
#ifdef  _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,nl,np,dx,dy,dtheta,h,hi,radius,pi,dt,rei,dxi,dyi,xo,yo
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: delta_x(nx,ny,nl,np),delta_y(nx,ny,nl,np) !デルタ関数
REAL(nk),INTENT(OUT) :: Udl(nl,np),Vdl(nl,np)                     !Lagrange点の物体速度
REAL(nk),INTENT(OUT) :: Xl(nl,np),Yl(nl,np)                       !Lagrange点の座標
INTEGER,INTENT(OUT) :: ij(nl,np,2)                                !物体付近の格子
REAL(nk),INTENT(IN) :: pos(2,np),angle(np),vel(2,np),avel(np)     !各物体の物理量

REAL(nk) :: absr                  !Euler点とLagrange点の距離的なもの
REAL(nk) :: phi                   !デルタ関数の一部
REAL(nk) :: x,y                   !Euler座標
INTEGER :: i,j                    !格子
INTEGER :: l                      !表面要素
INTEGER :: m                      !物体の数
!=============================================================================! 


!==Lagrange点の物体位置と速度==!
!#ifdef _TIMER_
!  CALL timer_start(410)
!#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
  Xl(l,m) = xo+radius*cos(dtheta*l)!pos(1,m)+radius*cos(dtheta*l)
  Yl(l,m) = yo+radius*sin(dtheta*l)!pos(2,m)+radius*sin(dtheta*l)
  Udl(l,m) = 0!vel(1,m)-avel(m)*radius*sin(dtheta*l)
  Vdl(l,m) = 0!vel(2,m)+avel(m)*radius*cos(dtheta*l)
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(410)
!#endif


!==物体周りの格子==!
!#ifdef _TIMER_
!  CALL timer_start(420)
!#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
  ij(l,m,1)=INT(Xl(l,m)*dxi)-5
  ij(l,m,2)=INT(Yl(l,m)*dyi)-5
  IF(ij(l,m,1)   <1 ) ij(l,m,1)=1
  IF(ij(l,m,2)   <1 ) ij(l,m,2)=1
  IF(ij(l,m,1)+10>nx) ij(l,m,1)=nx
  IF(ij(l,m,2)+10>ny) ij(l,m,2)=ny
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(420)
!#endif

!==デルタ関数==!
!#ifdef _TIMER_
!  CALL timer_start(430)
!#endif
!$OMP PARALLEL DO PRIVATE(x,y,absr,phi) COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
DO j = ij(l,m,2), ij(l,m,2)+10
DO i = ij(l,m,1), ij(l,m,1)+10
  !=x方向速度のためのデルタ関数=!
  x=i*dx ; y=(j-0.5_nk)*dy                      !座標

  absr = ABS(x-Xl(l,m))*hi                     
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_x(i,j,l,m) = phi*hi
  
  absr = ABS(y-Yl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_x(i,j,l,m) = delta_x(i,j,l,m)*phi*hi
  
  !=y方向速度のためのデルタ関数=!
  x=(i-0.5_nk)*dx ; y=j*dy                      !座標

  absr = ABS(x-Xl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_y(i,j,l,m) = phi*hi
  
  absr = ABS(y-Yl(l,m))*hi
  phi = 0
  IF(0.5 .le. absr .and. absr .le. 1.5)THEN
    phi = (5-3*absr-sqrt(-3*(1-absr)**2+1))/6
  ENDIF
  IF(absr .le. 0.5)THEN
    phi = (1+sqrt(-3*absr**2+1))/3
  ENDIF
  
  delta_y(i,j,l,m) = delta_y(i,j,l,m)*phi*hi
  
ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO
!#ifdef _TIMER_
!  CALL timer_end(430)
!#endif


END SUBROUTINE ibm
