!*****************************************************************************!
  MODULE parameter
!*****************************************************************************!
IMPLICIT NONE
!==precision==!
INTEGER,PARAMETER :: nk = 8

!==pi==!
REAL(16),PARAMETER :: piq = 3.1415926535897932384626433832795028_16
REAL(nk),PARAMETER :: pi = piq 

!==system parameters==!
INTEGER,PARAMETER :: nx = 500             !num. of grid for stream-wise direction(x)         
INTEGER,PARAMETER :: ny = 100             !num. of grid for the other direction(y)
REAL(16),PARAMETER :: lxq = 50            !system size for x 
REAL(16),PARAMETER :: lyq = 10            !system size for y
REAL(nk),PARAMETER :: lx=lxq, ly=lyq
REAL(nk),PARAMETER :: xo=40, yo=lyq/2     !cylinder posotion

!==time step==!
REAL(16),PARAMETER :: dtq = 0.01Q0!lyq/ny/10/2        !time increment
REAL(nk),PARAMETER :: dt = dtq, dti = 1 / dtq         !
INTEGER,PARAMETER :: total_time=50                  !total time
INTEGER,PARAMETER :: max_step =5000!total_time/dtq+1Q-8 !total time-step
INTEGER,PARAMETER :: file_step = 1/dtq                !

!==Reynolds number==!
REAL(16),PARAMETER :: req = 100,reiq = 1/req
REAL(nk),PARAMETER :: re = req, rei=1/req

!==for code optimaization==!
REAL(nk),PARAMETER :: dx=lxq/nx, dy=lyq/ny, h=dx
REAL(nk),PARAMETER :: dxi=nx/lxq, dyi=ny/lyq, hi=dxi
REAL(nk),PARAMETER :: dx2i=(nx/lxq)**2, dy2i=(ny/lyq)**2
REAL(nk),PARAMETER :: hhdxi=(nx/lxq)/4, hhdyi=(ny/lyq)/4
REAL(nk),PARAMETER :: dx2rei=(nx/lxq)**2/req, dy2rei=(ny/lyq)**2/req                                       

!==cylinder==!
REAL(nk),PARAMETER :: radius=1                   !radius
REAL(nk),PARAMETER :: rho=2                      !density ratio
REAL(nk),PARAMETER :: inertia=piq*radius**4/2    !moment of inertia/density
REAL(nk),PARAMETER :: area = piq*radius**2       !surface

!==IB method==!
INTEGER,PARAMETER :: np = 1                      !num of object
INTEGER,PARAMETER :: nl = 2*piq*radius/h         !num of surface element
REAL(nk),PARAMETER :: dtheta=2*piq/nl            !interval of surface element
REAL(nk),PARAMETER :: dv = h**2                  !area of surface element

!==RK method==!
REAL(nk) :: alpha(3),gamma(3),zeta(3)

!==FFTW==!
INTEGER(8) :: planf, plani                        !plan
REAL(nk) :: fftp(nx,ny)                           !
COMPLEX(nk) :: ffts(0:nx/2,0:ny-1)                !
COMPLEX(nk) :: laplacei(0:nx/2,0:ny-1,0:3)        !inverse Laplace
REAL(nk),PARAMETER :: nxnyi = 1.D0 / ( nx * ny )  !for normalization 

END MODULE parameter
