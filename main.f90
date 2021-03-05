!*****************************************************************************!
  PROGRAM main
!*****************************************************************************!
#ifdef _TIMER_
USE module_timer
#endif
USE parameter, ONLY:nk,nx,ny,nl,np,max_step,xo,yo,dt,alpha
IMPLICIT NONE
!==time==!
INTEGER :: n                                          !time step

!==fluid==!
REAL(nk) :: u(nx,ny),v(nx,ny)                         !velocity
REAL(nk) :: p(nx,ny)                                  !pressure
REAL(nk) :: u1(nx,ny),v1(nx,ny),u2(nx,ny),v2(nx,ny)   !velocity for RK method

!==cylinder==!
REAL(nk) :: pos(2,np)   !position
REAL(nk) :: angle(np)   !angular
REAL(nk) :: vel(2,np)   !velocity
REAL(nk) :: avel(np)    !angular velocity
REAL(nk) :: theta       !forcing direction
REAL(nk) :: drag(2)     !drag

!==IB method==!
REAL(nk) :: delta_x(nx,ny,nl,np),delta_y(nx,ny,nl,np) !delta function
REAL(nk) :: Udl(nl,np),Vdl(nl,np)                     !velocity at Lagrange points
REAL(nk) :: Xl(nl,np),Yl(nl,np)   !positions of  Lagrange points
INTEGER :: ij(nl,np,2)

!==time measurement==!
INTEGER(8) :: ns,t0,t1,t2,t_rate,t_max
INTEGER :: err


!==Python==!
CHARACTER(LEN=100) :: filename
INTEGER(8) :: file_size
INTEGER :: episode
REAL(nk) :: keep_state(2,0:1000)
INTEGER :: sample,i,j
INTEGER :: recl=8
!=============================================================================!
#ifdef _TIMER_
  CALL timer_init
  CALL timer_start(100)
  CALL timer_start(200)
#endif

!=================initial condition===================!

OPEN(200,FILE='INITIAL_u.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*nx*ny)
OPEN(201,FILE='INITIAL_v.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*nx*ny)
OPEN(202,FILE='INITIAL_pos.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*2*np)
OPEN(203,FILE='INITIAL_vel.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*2*np)
OPEN(204,FILE='INITIAL_keep.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*2*1001)
OPEN(205,FILE='INITIAL_angle.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*np)
OPEN(206,FILE='INITIAL_avel.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*np)
OPEN(207,FILE='INITIAL_p.bin',FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl*nx*ny)
  READ(200,REC=1) u
  READ(201,REC=1) v
  READ(202,REC=1) pos
  READ(203,REC=1) vel
  READ(204,REC=1) keep_state
  READ(205,REC=1) angle
  READ(206,REC=1) avel
  READ(207,REC=1) p
CLOSE(200);CLOSE(201);CLOSE(202);CLOSE(203)
CLOSE(204);CLOSE(205);CLOSE(206);CLOSE(207)


!==setting of fftw etc==!
CALL set_all
sample=1000
#ifdef _TIMER_
  CALL timer_end(200)
#endif

!==================IB==================!
#ifdef _TIMER_
  CALL timer_start(400)
#endif
CALL ibm(delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,pos,angle,vel,avel)  ! setting of IB method
#ifdef _TIMER_
  CALL timer_end(400)
#endif



#ifdef _TIMER_
  CALL timer_start(300)
#endif
CALL system_clock(t1) 

DO   ! iteration of episode

  file_size=0

  READ(*,*) episode
  !WRITE(filename,'(i4.4)') episode
!  OPEN(1000,FILE="action.bin",FORM='BINARY',STATUS='old',ACCESS='DIRECT',RECL=8)
!  OPEN(1001,FILE="observe.bin",FORM='BINARY',STATUS='replace')
!  OPEN(1000,FILE="action.bin",FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl)
  OPEN(1001,FILE="observe.bin",FORM='UNFORMATTED',ACCESS='DIRECT',RECL=recl*9,STATUS='replace') ; CLOSE(1001)

  WRITE(*,*) 'INPUT0'

  DO n = 1, max_step
    WRITE(*,*) "==========",n,"=========="
    
    !===action===!
    WRITE(*,*) 'INPUT'
    DO WHILE(file_size<8*n)
      !inquire(unit=1000, size=file_size)
      open(1000, file='action.bin', form='unformatted', access='stream',status='old', action='read', position='append')
      inquire(1000, pos=file_size )
      file_size = file_size - 1
      close(1000)
    ENDDO
    OPEN(1000,FILE="action.bin",FORM='UNFORMATTED',STATUS='old',ACCESS='DIRECT',RECL=recl, action='read')
    READ(1000,REC=n) theta ; CLOSE(1000)
    
    !===RK method===!
#ifdef _TIMER_
  CALL timer_start(500)
#endif
    CALL rg(u1,v1,p,u ,v ,u ,v ,pos,angle,vel,avel,theta,drag,delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,1)    !RK 1st step
    CALL rg(u2,v2,p,u1,v1,u ,v ,pos,angle,vel,avel,theta,drag,delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,2)    !RK 2nd step
    CALL rg(u ,v ,p,u2,v2,u1,v1,pos,angle,vel,avel,theta,drag,delta_x,delta_y,Udl,Vdl,Xl,Yl,ij,3)    !RK 3rd step
#ifdef _TIMER_
  CALL timer_end(500)
#endif

    !===for time-delay coodinate===!
    DO i = 0, sample-1
      keep_state(:,i) =  keep_state(:,i+1)
    ENDDO
    keep_state(1,sample)=vel(2,1)
    keep_state(2,sample)=drag(2)!SIN(angle(1))
    !keep_state(3,sample)=COS(angle(1))
    !keep_state(4,sample)=avel(1)

    !===send reward and next state===!

    !WRITE(1001) vel(1,1),vel(2,1),drag(2),keep_state(:,900),keep_state(:,800),keep_state(:,700)!,keep_state(1,600)
    OPEN(1001,FILE="observe.bin",FORM='UNFORMATTED',ACCESS='DIRECT',RECL=recl*9,STATUS='replace')
    WRITE(1001,REC=n) vel(1,1),vel(2,1),drag(2),keep_state(:,900),keep_state(:,800),keep_state(:,700)!,keep_state(1,600)
    CLOSE(1001)
    WRITE(*,*) 'OUTPUT'

  ENDDO


  WRITE(*,*) 'END'
  CLOSE(1000);CLOSE(1001);CLOSE(1002)
  CLOSE(100);CLOSE(101);CLOSE(102);CLOSE(103);CLOSE(104);CLOSE(105);CLOSE(106);CLOSE(107);CLOSE(108)

ENDDO

CALL system_clock(t2,t_rate,t_max)
#ifdef _TIMER_
  CALL timer_end(300)
  CALL timer_end(100)
  CALL timer_finalize
#endif
WRITE(*,*) (t2-t1)/dble(t_rate)


END
