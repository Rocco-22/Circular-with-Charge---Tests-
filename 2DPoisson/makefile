.SUFFIXES:
.SUFFIXES: .f90 .c .o
#
O3= -O2

OPENMP= -qopenmp

CC= gcc

FOPT= -c -O2 -fno-alias  \
	-ftz -IPF_fma -IPF_fltacc \
	-complex-limited-range \
	-check noarg_temp_created \
	-qopenmp \
	-auto 

FDEB= -c -O0 -check all -check noarg_temp_created -auto -g -qopenmp

FC= ifx -ftz

LINK= $(FC) -g $(OPENMP)
F9FLAGS = $(FFLAGS) 

RBINS= ReadInputDataAndSolve.o Multigrid.o \
  SpecifyPhiAtxy.o GaussSeidel2D.o NaNs_etc.o \
  gcontr.o ModulePlots.o lszeic.o \
  isnanC.o

R2BINS= ReadPhiData.o gcontr.o ModulePlots.o lszeic.o

.f90.o:
	$(FC) $(FOPT) $< -o $@

.c.o:
	gcc -c $< -o $@

a.out: Kinds.o $(RBINS)
	$(LINK) Kinds.o $(RBINS) -o $@
	./a.out < Gitterfuellung

b.out: Kinds.o $(R2BINS)
	$(LINK) Kinds.o $(R2BINS) -o $@
	./b.out < PhiData

# Dependencies:
ReadInputDataAndSolve.o: Multigrid.o ModulePlots.o

Multigrid.o: NaNs_etc.o

MThreads.o: MaxCPUsETC.o

GaussSeidel2D.o: Multigrid.o
GaussSeidel2DN.o: Multigrid.o

lszeic.o: ModulePlots.o

# All which depend on Kinds.o:
$(RBINS): Kinds.o
$(R2BINS): Kinds.o
