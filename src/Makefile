SRCS := m_octree_mg.f90 m_morton.f90
OBJS := $(SRCS:%.f90=%.o)
VPATH := morton_fortran
PROGS = test_basics_2d

.PHONY: all clean

all: $(PROGS)

clean:
	$(RM) $(OBJS) $(PROGS) *.mod

FC 	:= gfortran
FFLAGS	:= -O2 -std=f2008 -fopenmp -Wall -Wextra -Wimplicit-interface

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=all -g -ffpe-trap=invalid,zero,overflow \
	-pedantic -finit-real=snan
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Depedencies
$(PROGS:%=%.o): $(OBJS)
$(PROGS): $(OBJS)
m_octree_mg.o: m_morton.o
