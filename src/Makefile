SRCS := m_data_structures.f90 m_build_tree.f90 m_load_balance.f90	\
m_ghost_cells.f90 m_allocate_storage.f90 m_mrgrnk.f90 m_restrict.f90	\
m_communication.f90 m_prolong.f90 m_multigrid.f90
OBJS := $(SRCS:%.f90=%.o)
VPATH := morton_fortran
PROGS = test_uniform_grid

.PHONY: all clean

all: $(PROGS)

clean:
	$(RM) $(OBJS) $(PROGS) *.mod

FC 	:= mpif90
FFLAGS	:= -O2 -std=f2008 -fopenmp -Wall -Wextra

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
m_ghost_cells.o: m_data_structures.o m_communication.o
m_prolong.o: m_data_structures.o m_communication.o
m_restrict.o: m_data_structures.o m_communication.o
m_allocate_storage.o: m_data_structures.o m_ghost_cells.o m_restrict.o m_prolong.o
m_communication.o: m_mrgrnk.o
m_multigrid.o: m_data_structures.o m_prolong.o m_restrict.o
