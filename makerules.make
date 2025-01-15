F90 = mpif90

# Determine compiler brand
compiler_version = $(shell $(F90) --version)
compiler_brand = $(word 1, $(compiler_version))
$(info $(compiler_brand))

ifeq ($(compiler_brand), GNU)
	F90FLAGS ?= -O2 -std=f2008 -fopenmp -Wall -g -Wno-unused-dummy-argument	\
	-Wno-unused-function -cpp
# Add debugging flags when DEBUG=1
	ifeq ($(DEBUG), 1)
		F90FLAGS += -fcheck=all -ffpe-trap=invalid,zero,overflow \
		-pedantic -finit-real=snan
	endif
else ifeq ($(compiler_brand), ifort)
	F90FLAGS ?= -warn all -O2 -stand f08 -assume realloc-lhs -fpp
else ifeq ($(compiler_brand), nvfortran)
	F90FLAGS ?= -Wall -O2 -Mpreprocess
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(F90) -c -o $@ $< $(F90FLAGS) -DNDIM=$(NDIM) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(F90) -c -o $(@:.mod=.o) $< $(F90FLAGS) -DNDIM=$(NDIM) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=$(NDIM) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
