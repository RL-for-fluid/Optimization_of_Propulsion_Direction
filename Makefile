NAME=test

# use module_timer.f90
TIMER=off

#module switch PrgEnv-cray PrgEnv-intel
#

#FC=ifort
#FFLAGS=-O3  -axMIC-AVX512 -mcmodel=large -qopenmp -assume byterecl #-parallel
#FFLAGS=-O0   -CB -fpe0 -traceback -g -check uninit #-check noarg_temp_created #-heap-arrays -warn -assume byterecl
FC=gfortran
FFLAGS=-O2
FFLAGS2=$(FFLAGS)
MPIC=mpirun -np
FFTW_DIR= ./FFTW/lib

ifeq ("$(TIMER)","on")
  FFLAGS+=-D_TIMER_ -D_TIMER_NOMPI_
endif
FFLAGS2=$(FFLAGS)
AR="ar -r"
MOD_DIR = ./
LIB_DIR = ./
FFLAGS+=-cpp


OBJS = parameter.o main.o set_all.o poisson.o helmholtz.o ibm.o rg.o

ifeq ("$(TIMER)","on")
OBJS +=module_timer.o
endif

#.SUFFIXES : .f90
%.o : %.f90
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.f90
#.SUFFIXES : .f
%.o : %.f
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.f
#.SUFFIXES : .F
%.o : %.F
	$(FC) $(FFLAGS) -I$(MOD_DIR) -c  $*.F

#for ifort
#$(NAME).x :  $(OBJS)
#	$(FC) $(FFLAGS2) -o  $(NAME).x $(OBJS) -mkl
#for gfort
$(NAME).x :  $(OBJS)        \
	$(FFTW_DIR)/libfftw3.a
	$(FC) $(FFLAGS2) -o  $(NAME).x $(OBJS) \
	-L$(FFTW_DIR) -lfftw3

clean :
	rm -f $(OBJS)  *~ core*  *.mod *.o fort.* *.lst *.L  $(NAME).x

run  : $(NAME).x
	./$(NAME).x

ifeq ("$(TIMER)","on")
main.o:module_timer.o
march.o:module_timer.o
endif
