NDIM ?= 2

SRCS := m_data_structures.f90 m_build_tree.f90 m_load_balance.f90	\
m_ghost_cells.f90 m_allocate_storage.f90 m_mrgrnk.f90 m_restrict.f90	\
m_communication.f90 m_prolong.f90 m_multigrid.f90 m_octree_mg.f90
OBJS := $(SRCS:%.f90=%.o)
VPATH := morton_fortran
PROGS = test_uniform_grid

.PHONY: all clean

all: $(PROGS) liboctree-mg.a

clean:
	$(RM) $(OBJS) $(PROGS) *.mod

FC 	:= mpif90
FFLAGS	:= -O2 -std=f2008 -fopenmp -Wall -Wextra -g -cpp -DNDIM=$(NDIM)

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=all -ffpe-trap=invalid,zero,overflow \
	-pedantic -finit-real=snan
endif

ifeq ($(PROF), 1)
	FFLAGS += -pg -O0
endif

liboctree-mg.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Depedencies
$(PROGS:%=%.o): $(OBJS)
$(PROGS): $(OBJS)

m_allocate_storage.o: m_data_structures.mod
m_allocate_storage.o: m_ghost_cells.mod
m_allocate_storage.o: m_prolong.mod
m_allocate_storage.o: m_restrict.mod
m_build_tree.o: m_data_structures.mod
m_communication.o: m_data_structures.mod
m_communication.o: m_mrgrnk.mod
m_fishpack.o: m_data_structures.mod
m_ghost_cells.o: m_communication.mod
m_ghost_cells.o: m_data_structures.mod
m_load_balance.o: m_data_structures.mod
m_multigrid.o: m_data_structures.mod
m_multigrid.o: m_ghost_cells.mod
m_multigrid.o: m_prolong.mod
m_multigrid.o: m_restrict.mod
m_octree_mg.o: m_allocate_storage.mod
m_octree_mg.o: m_build_tree.mod
m_octree_mg.o: m_communication.mod
m_octree_mg.o: m_data_structures.mod
m_octree_mg.o: m_ghost_cells.mod
m_octree_mg.o: m_load_balance.mod
m_octree_mg.o: m_multigrid.mod
m_octree_mg.o: m_prolong.mod
m_octree_mg.o: m_restrict.mod
m_prolong.o: m_communication.mod
m_prolong.o: m_data_structures.mod
m_restrict.o: m_communication.mod
m_restrict.o: m_data_structures.mod
test_morton.o: m_morton.mod
test_uniform_grid.o: m_octree_mg.mod
