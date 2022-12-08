#FFLAGS=-O3 -fopenmp
# FFLAGS=-fbounds-check -fbacktrace -g
FFLAGS=-mcmodel=medium -O3 -fopenmp

F90SRCS =           \
bndset.f90 			\
calphi.f90			\
init.f90			\
nabla.f90			\
main.f90            \

SRCS  =  $(F90SRCS)

.SUFFIXES: .o .f .f90
# FCOBJS = $(FCSRCS:.f=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  =  $(F90OBJS)

.f90.o:
	gfortran -c $(FFLAGS) $<

.f.o:
	gfortran -c $(FFLAGS) $<

a.out: $(OBJS)
	gfortran $(OBJS) $(FFLAGS) -o a.out

clean:
	rm $(OBJS)
	rm a.out
