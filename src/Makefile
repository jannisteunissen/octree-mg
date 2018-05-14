# This makefile is included in the ../lib_2d and ../lib_3d folder makefiles
ifndef NDIM
$(error NDIM is not set, type make in top or the lib folders)
endif

SRCS := m_data_structures.f90 m_build_tree.f90 m_load_balance.f90	\
m_ghost_cells.f90 m_allocate_storage.f90 m_mrgrnk.f90 m_restrict.f90	\
m_communication.f90 m_prolong.f90 m_multigrid.f90 m_octree_mg.f90	\
m_laplacian.f90 m_vlaplacian.f90 m_helmholtz.f90 m_vhelmholtz.f90	\
m_diffusion.f90

OBJS := $(SRCS:%.f90=%.o)
LIB := libomg.a

include ../makerules.make

.PHONY: all clean

all: $(LIB)

clean:
	$(RM) $(OBJS) *.mod

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

m_allocate_storage.o: m_data_structures.mod
m_allocate_storage.o: m_ghost_cells.mod
m_allocate_storage.o: m_prolong.mod
m_allocate_storage.o: m_restrict.mod
m_build_tree.o: m_data_structures.mod
m_communication.o: m_data_structures.mod
m_communication.o: m_mrgrnk.mod
m_diffusion.o: m_data_structures.mod
m_diffusion.o: m_helmholtz.mod
m_diffusion.o: m_multigrid.mod
m_diffusion.o: m_vhelmholtz.mod
m_fishpack.o: m_data_structures.mod
m_ghost_cells.o: m_communication.mod
m_ghost_cells.o: m_data_structures.mod
m_helmholtz.o: m_data_structures.mod
m_laplacian.o: m_data_structures.mod
m_load_balance.o: m_data_structures.mod
m_multigrid.o: m_data_structures.mod
m_multigrid.o: m_ghost_cells.mod
m_multigrid.o: m_helmholtz.mod
m_multigrid.o: m_laplacian.mod
m_multigrid.o: m_prolong.mod
m_multigrid.o: m_restrict.mod
m_multigrid.o: m_vhelmholtz.mod
m_multigrid.o: m_vlaplacian.mod
m_octree_mg.o: m_allocate_storage.mod
m_octree_mg.o: m_build_tree.mod
m_octree_mg.o: m_communication.mod
m_octree_mg.o: m_data_structures.mod
m_octree_mg.o: m_ghost_cells.mod
m_octree_mg.o: m_helmholtz.mod
m_octree_mg.o: m_load_balance.mod
m_octree_mg.o: m_multigrid.mod
m_octree_mg.o: m_prolong.mod
m_octree_mg.o: m_restrict.mod
m_octree_mg.o: m_vhelmholtz.mod
m_prolong.o: m_communication.mod
m_prolong.o: m_data_structures.mod
m_restrict.o: m_communication.mod
m_restrict.o: m_data_structures.mod
m_vhelmholtz.o: m_data_structures.mod
m_vlaplacian.o: m_data_structures.mod
