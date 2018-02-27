SRCS := m_data_structures.f90 m_build_tree.f90 m_load_balance.f90
OBJS := $(SRCS:%.f90=%.o)
VPATH := morton_fortran
PROGS = test_one_level

.PHONY: all clean

all: $(PROGS)

clean:
	$(RM) $(OBJS) $(PROGS) *.mod

FC 	:= mpif90
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

m_build_tree.o: m_data_structures.o
m_load_balance.o: m_data_structures.o
