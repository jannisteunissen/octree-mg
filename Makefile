# Makefile updated and simplified by Jannis Teunissen
F90 = mpif90
F90FLAGS = -O2 -cpp

# Create a static library
LIBNAME := poisfree
LIB := lib$(LIBNAME).a

PROGRAMS := PSolver

all: $(PROGRAMS)

clean:
	\rm -f *.o ${PROGRAMS} *~ *.log *.aux fort.* *.mod \
	ABINIT-common/*.f90 ABINIT-common/*.o ABINIT-common/*.mod

cleanall:
	make clean
	\rm -rf *.dat *.pdf gmon.out

SOURCES_ABINIT := ABINIT-common/defs_basis.F90 \
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

SOURCES_F90 := Poisson_Solver.f90 timing.f90

SOURCES_F77 := dcopy.f

OBJECTS := $(SOURCES_ABINIT:.F90=.o) $(SOURCES_F90:.f90=.o) $(SOURCES_F77:.f=.o)

# Rules for compiling object files
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(F90FLAGS) -c $< -o $@

# Rule for compiling programs
%: %.f90
	$(F90) $(F90FLAGS) $< -o $@ -l $(LIBNAME) -L .

$(LIB): $(OBJECTS) Poisson_Solver.o
	$(RM) $@
	$(AR) rcs $@ $^

# Dependencies
$(PROGRAMS): $(LIB)

# Avoid circular dependency
$(filter-out ABINIT-common/defs_basis.o, $(SOURCES_ABINIT:.F90=.o)): \
	ABINIT-common/defs_basis.o

ABINIT-common/defs_datatypes.o: ABINIT-common/defs_basis.o

ABINIT-common/defs_xc.o: ABINIT-common/defs_datatypes.o

ABINIT-common/drivexc.o: ABINIT-common/defs_basis.o ABINIT-common/defs_xc.o
