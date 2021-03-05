!*****************************************************************************!
  SUBROUTINE rg(un,vn,p,u,v,up,vp,pos,angle,vel,avel,theta,drag,delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,k)
!*****************************************************************************!
#ifdef _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,nl,np,dxi,dyi,hhdxi,hhdyi,dx2rei,dy2rei,dx2i,dy2i, &
                    alpha,gamma,zeta,dt,dti,h,dv,rei,rho,inertia,area,dx,dy
IMPLICIT NONE
REAL(nk),INTENT(OUT) :: un(nx,ny),vn(nx,ny)                      !次ルンゲクッタ時刻の速度
REAL(nk),INTENT(INOUT) :: p(nx,ny)                               !圧力
REAL(nk),INTENT(IN) :: u(nx,ny),v(nx,ny)                         !現ルンゲクッタ時刻の速度
REAL(nk),INTENT(IN) :: up(nx,ny),vp(nx,ny)                       !前ルンゲクッタ時刻の速度
REAL(nk),INTENT(INOUT) :: pos(2,np),angle(np),vel(2,np),avel(np)    !各物体の物理量
REAL(nk),INTENT(IN) :: theta
REAL(nk),INTENT(INOUT) :: drag(2)
REAL(nk),INTENT(IN) :: delta_x(nx,ny,nl,np),delta_y(nx,ny,nl,np) !デルタ関数
REAL(nk),INTENT(IN) :: Udl(nl,np),Vdl(nl,np)                     !Lagrange点の物体速度
REAL(nk),INTENT(IN) :: Xl(nl,np),Yl(nl,np)   !Lagrange点の座標
INTEGER,INTENT(IN):: ij(nl,np,2)
INTEGER,INTENT(IN) :: k                                          !ルンゲクッタのstep

REAL(nk) :: ut(nx,ny),vt(nx,ny)     !速度^tilde
REAL(nk) :: Utl(nl,np),Vtl(nl,np)   !Lagrange点の速度^tilde
REAL(nk) :: Ful(nl,np),Fvl(nl,np)   !Lagrange点のibm用の力
REAL(nk) :: fu(nx,ny),fv(nx,ny)     !ibm用の力
REAL(nk) :: ua(nx,ny),va(nx,ny)     !速度^ast
REAL(nk) :: phi(nx,ny)              !擬似圧力
INTEGER :: i,j,im,ip,jm,jp          !格子
INTEGER :: l                        !表面要素
INTEGER :: m                        !物体の個数

REAL(nk) :: sumFul(np),sumFvl(np),sumFl(np)
REAL(nk) :: x,y
!=============================================================================!


!===(12 a)===!
#ifdef _TIMER_
  CALL timer_start(510)
#endif
!$OMP PARALLEL DO PRIVATE(im,ip,jm,jp) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  im = MOD(i-2+nx,nx)+1 ; ip = MOD(i,nx)+1
  jm = MOD(j-2+ny,ny)+1 ; jp = MOD(j,ny)+1
  ut(i,j) =  &
    ( (u(im,j)-2*u(i,j)+u(ip,j))*dx2rei + (u(i,jm)-2*u(i,j)+u(i,jp))*dy2rei )*2*alpha(k) &
   + (p(i,j)-p(ip,j))*dxi                                                    *2*alpha(k) &
   +( hhdxi*(( u(im,j)+u(i,j) )**2                  -( u(i,j)+u(ip,j) )**2                ) &
    + hhdyi*(( u(i,jm)+u(i,j) )*( v(i,jm)+v(ip,jm) )-( u(i,j)+u(i,jp) )*( v(i,j)+v(ip,j) )) )*gamma(k) &
   +( hhdxi*((up(im,j)+up(i,j))**2                  -(up(i,j)+up(ip,j))**2                ) &
    + hhdyi*((up(i,jm)+up(i,j))*(vp(i,jm)+vp(ip,jm))-(up(i,j)+up(i,jp))*(vp(i,j)+vp(ip,j))) )*zeta(k)
  ut(i,j) =  u(i,j) + ut(i,j)*dt

  vt(i,j) =  &
    ( (v(im,j)-2*v(i,j)+v(ip,j))*dx2rei + (v(i,jm)-2*v(i,j)+v(i,jp))*dy2rei )*2*alpha(k) &
   + (p(i,j)-p(i,jp))*dyi                                                    *2*alpha(k) &
   +( hhdxi*(( u(im,j)+u(im,jp) )*( v(im,j)+v(i,j) )-( u(i,j)+u(i,jp) )*( v(i,j)+v(ip,j) )) &
    + hhdyi*(( v(i,jm)+v(i,j)   )**2                -( v(i,j)+v(i,jp) )**2                ) )*gamma(k) &
   +( hhdxi*((up(im,j)+up(im,jp))*(vp(im,j)+vp(i,j))-(up(i,j)+up(i,jp))*(vp(i,j)+vp(ip,j))) &
    + hhdyi*((vp(i,jm)+vp(i,j)  )**2                -(vp(i,j)+vp(i,jp))**2                ) )*zeta(k)
  vt(i,j) =  v(i,j) + vt(i,j)*dt
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(510)
#endif


!===(12 b)===!
#ifdef _TIMER_
  CALL timer_start(520)
#endif
!$OMP PARALLEL DO COLLAPSE(2)
DO m = 1, np
DO l = 1, nl
  Utl(l,m)=0
  Vtl(l,m)=0
ENDDO
ENDDO
!$OMP END PARALLEL DO
DO m = 1, np
DO l = 1, nl
DO j = ij(l,m,2), ij(l,m,2)+10
DO i = ij(l,m,1), ij(l,m,1)+10
    Utl(l,m)=Utl(l,m)+ut(i,j)*delta_x(i,j,l,m)*h**2
    Vtl(l,m)=Vtl(l,m)+vt(i,j)*delta_y(i,j,l,m)*h**2
ENDDO
ENDDO
ENDDO
ENDDO
#ifdef _TIMER_
  CALL timer_end(520)
#endif


!===(12 c)===!
#ifdef _TIMER_
  CALL timer_start(530)
#endif
DO m = 1, np
DO l = 1, nl
  Ful(l,m) = (Udl(l,m)-Utl(l,m))*dti
  Fvl(l,m) = (Vdl(l,m)-Vtl(l,m))*dti
ENDDO
ENDDO
#ifdef _TIMER_
  CALL timer_end(530)
#endif

!===(12 ')===!
#ifdef _TIMER_
  CALL timer_start(540)
#endif
DO m = 1, np
  sumFul(m) = 0
  sumFvl(m) = 0
ENDDO
DO m = 1, np
DO l = 1, nl
  sumFul(m) = sumFul(m) + Ful(l,m)*dv
  sumFvl(m) = sumFvl(m) + Fvl(l,m)*dv
ENDDO
ENDDO

!===(12 d)===!
!$OMP PARALLEL DO COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  fu(i,j)= sumFul(1)/(area*(rho-1)) - 2*alpha(k)*COS(theta)/(area*(rho-1))
  fv(i,j)= sumFvl(1)/(area*(rho-1)) - 2*alpha(k)*SIN(theta)/(area*(rho-1))
ENDDO
ENDDO
!$OMP END PARALLEL DO
DO m = 1, np
DO l = 1, nl
DO j = ij(l,m,2), ij(l,m,2)+10
DO i = ij(l,m,1), ij(l,m,1)+10
    fu(i,j)=fu(i,j)+Ful(l,m)*delta_x(i,j,l,m)*dv 
    fv(i,j)=fv(i,j)+Fvl(l,m)*delta_y(i,j,l,m)*dv 
ENDDO
ENDDO
ENDDO
ENDDO
#ifdef _TIMER_
  CALL timer_end(540)
#endif

!===(12 e)===!
#ifdef _TIMER_
  CALL timer_start(550)
#endif
CALL helmholtz(ua,ut,fu,u,k)
#ifdef _TIMER_
  CALL timer_end(550)
#endif
#ifdef _TIMER_
  CALL timer_start(550)
#endif
CALL helmholtz(va,vt,fv,v,k)
#ifdef _TIMER_
  CALL timer_end(550)
#endif



!===(12 f)===!
#ifdef _TIMER_
  CALL timer_start(560)
#endif
CALL poisson(phi,ua,va,k)
#ifdef _TIMER_
  CALL timer_end(560)
#endif


!===(12 g)===!
!===(12 h)===!
#ifdef _TIMER_
  CALL timer_start(570)
#endif
!$OMP PARALLEL DO PRIVATE(im,ip,jm,jp) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  im = MOD(i-2+nx,nx)+1 ; ip = MOD(i,nx)+1
  jm = MOD(j-2+ny,ny)+1 ; jp = MOD(j,ny)+1

  un(i,j) = ua(i,j) - 2*alpha(k)*dt*(phi(ip,j)-phi(i,j))*dxi
  vn(i,j) = va(i,j) - 2*alpha(k)*dt*(phi(i,jp)-phi(i,j))*dyi

  p(i,j)  = p(i,j) + phi(i,j)-alpha(k)*dt*rei* &
                  ( (phi(im,j)-2*phi(i,j)+phi(ip,j))* dx2i &
                   +(phi(i,jm)-2*phi(i,j)+phi(i,jp))* dy2i )
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(570)
#endif


!==円柱が受ける抗力==!
#ifdef _TIMER_
  CALL timer_start(580)
#endif
IF(k==1)THEN
  drag=0
ENDIF
drag(1) = drag(1) - sumFul(1)
drag(2) = drag(2) - sumFvl(1)

!===円柱の位置と速度の時間発展===!
DO m = 1, np
  pos(1,m) = pos(1,m) + alpha(k)*dt*vel(1,m)
  pos(2,m) = pos(2,m) + alpha(k)*dt*vel(2,m)
ENDDO
DO m = 1, np
  vel(1,m) = vel(1,m) - sumFul(m)*dt/(area*(rho-1)) + 2*alpha(k)*COS(theta)*dt/(area*(rho-1))
  vel(2,m) = vel(2,m) - sumFvl(m)*dt/(area*(rho-1)) + 2*alpha(k)*SIN(theta)*dt/(area*(rho-1))
ENDDO
DO m = 1, np
  pos(1,m) = pos(1,m) + alpha(k)*dt*vel(1,m)
  pos(2,m) = pos(2,m) + alpha(k)*dt*vel(2,m)
ENDDO

!===円柱の角度と角速度の時間発展===!
!DO m = 1, np
!  sumFl(m) = 0
!ENDDO
!DO m = 1, np
!DO l = 1, nl
!  sumFl(m) = sumFl(m) + ((Xl(l,m)-pos(1,m))*Fvl(l,m)-(Yl(l,m)-pos(2,m))*Ful(l,m))*dv
!ENDDO
!ENDDO
!
!DO m = 1, np
!  angle(m) = angle(m) + alpha(k)*dt*avel(m)
!ENDDO
!DO m = 1, np
!  avel(m) = avel(m) -sumFl(m)*dt/(inertia*(rho-1))
!ENDDO
!DO m = 1, np
!  angle(m) = angle(m) + alpha(k)*dt*avel(m)
!ENDDO
#ifdef _TIMER_
  CALL timer_end(580)
#endif

!==減衰==!
#ifdef _TIMER_
  CALL timer_start(590)
#endif
!$OMP PARALLEL DO PRIVATE(x,y) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  x=i*dx ; y=(j-0.5_nk)*dy                                !座標
  IF(x >= 45)THEN
    un(i,j) = -vel(1,1)
  ENDIF
ENDDO
ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(x,y) COLLAPSE(2)
DO j = 1, ny
DO i = 1, nx
  x=(i-0.5_nk)*dx ; y=j*dy                                !座標
  IF(x >= 45)THEN
    vn(i,j) = -vel(2,1)
  ENDIF
ENDDO
ENDDO
!$OMP END PARALLEL DO
#ifdef _TIMER_
  CALL timer_end(590)
#endif

END SUBROUTINE rg
