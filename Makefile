#***h* Poisson_Solver/Makefile
# NAME
#   Makefile
#
# FUNCTION
#    Make everything
#
# SOURCE
#
#switch between a parallel and a serial approach

#sequential treatment
mpi_source = MPIfake.f90
mpi_include = mpif.h

F95 = g95 -fbounds-check -O 
F77 = g95 -fbounds-check -O 

#parallel treatment
#mpi_source = 
#mpi_include = 

#F95 = ompif90 -fbounds-check -O 
#F77 = ompif90 -fbounds-check -O 

all: PSolver

clean:
	\rm -f *.o ${PROGRAMS} *~ *.log *.aux fort.* mpif.h *.mod \
	ABINIT-common/*.f90 ABINIT-common/*.o ABINIT-common/*.mod

cleanall:
	make clean
	\rm -rf *.dat *.pdf gmon.out

tar:
	tar czvf solver.tgz $(SOURCES_ABINIT) $(SOURCES_F77) $(SOURCES_F90) $(INCLUDES) $(OTHERS)

PROGRAMS = \
	PSolver




SOURCES_ABINIT = ABINIT-common/defs_basis.F90 \
        ABINIT-common/defs_datatypes.F90 \
        ABINIT-common/defs_xc.F90 \
        ABINIT-common/drivexc.F90 \
        ABINIT-common/xctetr.F90 \
        ABINIT-common/xcwign.F90 \
        ABINIT-common/xcxalp.F90 \
        ABINIT-common/xchelu.F90 \
        ABINIT-common/xcpzca.F90 \
        ABINIT-common/xcspol.F90 \
        ABINIT-common/xchcth.F90 \
        ABINIT-common/xclb.F90 \
        ABINIT-common/xcpbe.F90 \
        ABINIT-common/invcb.F90 \
	ABINIT-common/size_dvxc.F90 

SOURCES_F90 = \
	Poisson_Solver.f90 \
	PSolver.f90 \
	timing.f90 \
	$(mpi_source)


SOURCES_F77 = \
	dcopy.f

INCLUDES = perfdata.inc \
	lazy_100.inc \
	lazy_16.inc \
	lazy_24.inc \
	lazy_40.inc \
	lazy_60.inc \
	lazy_14.inc \
	lazy_20.inc \
	lazy_30.inc \
	lazy_50.inc \
	lazy_8.inc

SOURCES_MOD = PSolver_Main.f90 \
	Build_Kernel.f90 \
	PSolver_Base.f90 \
	xcenergy.f90 \
	3Dgradient.f90 \
	fft3d.f90 \
	scaling_function.f90

OTHERS = Makefile \
	PSchk.f90 \
	MPIfake.f90 \
	bench.csh \
	perfdata.t \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	COPYING \
	TODO \
	$(SOURCES_MOD)

OBJECTS_ABINIT_PP = $(SOURCES_ABINIT:.F90=.f90)

$(OBJECTS_ABINIT_PP): %.f90: %.F90
	cpp -P $< $@

OBJECTS_ABINIT = $(OBJECTS_ABINIT_PP:.f90=.o)

OBJECTS_F90 = $(SOURCES_F90:.f90=.o) 

OBJECTS_F77 = $(SOURCES_F77:.f=.o) 


$(OBJECTS_F90): %.o: %.f90
	$(F95) -c $< -o $@

$(OBJECTS_F77): %.o: %.f
	$(F95) -c $< -o $@


$(OBJECTS_ABINIT): %.o: %.f90
	$(F95) -c $< -o $@

COMMON =  $(OBJECTS_ABINIT) \
	$(OBJECTS_F90) \
	dcopy.o

PSolver: $(COMMON) $(INCLUDES) 
	$(F95) -o $@ $(COMMON)

PSolver.o : Poisson_Solver.o

Poisson_Solver.o : $(mpi_include) $(SOURCES_MOD)

mpif.h:
	touch mpif.h &&\
	 echo "integer :: MPI_SUM, MPI_COMM_WORLD" >> mpif.h &&\
	 echo "integer :: MPI_DOUBLE_PRECISION" >> mpif.h &&\
	 echo "integer :: MPI_MIN, MPI_MAX" >> mpif.h

$(SOURCES_ABINIT:.f90=.o) : ABINIT-common/defs_basis.o

ABINIT-common/defs_datatypes.o: ABINIT-common/defs_basis.o

ABINIT-common/defs_xc.o: ABINIT-common/defs_datatypes.o

ABINIT-common/drivexc.o: ABINIT-common/defs_basis.o ABINIT-common/defs_xc.o
