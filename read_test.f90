REAL(8) :: x(10) 
INTEGER :: i


DO i=1,10
  x(i)=i
ENDDO

OPEN(10,FILE='test.bin',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*10)
!OPEN(10,FILE='test.bin',FORM='BINARY',ACCESS='DIRECT',RECL=8*10)
WRITE(10,REC=1) x
CLOSE(10)
WRITE(*,*) x
x=0

OPEN(11,FILE='test.bin',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8, status='OLD')
!OPEN(11,FILE='test.bin',FORM='BINARY',ACCESS='DIRECT',RECL=8, status='OLD')
DO i=1,10
  READ(11,REC=i) x(i)
ENDDO

WRITE(*,*) x

END
