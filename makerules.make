# A custom compiler can be set with make F90C=...
F90C ?= mpif90

ifeq ($(F90C), mpif90)
FFLAGS = -O2 -std=f2008 -fopenmp -Wall -g -cpp -DNDIM=$(NDIM)	\
-Wno-unused-dummy-argument -Wno-unused-function

ifeq ($(DEBUG), 1)
FFLAGS += -fcheck=all -ffpe-trap=invalid,zero,overflow \
-pedantic -finit-real=snan
endif

ifeq ($(PROF), 1)
FFLAGS += -pg
endif
endif

ifeq ($(F90C), mpiifort)
FFLAGS	= -warn all -O2 -stand f08 -fpp -assume realloc-lhs -DNDIM=$(NDIM)
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90C) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(F90C) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(F90C) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

