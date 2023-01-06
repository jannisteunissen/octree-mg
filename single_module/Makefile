F90 := mpif90
F90FLAGS := -O2 -Wall -cpp
PROGS_1D := test_uniform_grid_1d test_refinement_1d
PROGS_2D := test_uniform_grid_2d test_refinement_2d
PROGS_3D := test_uniform_grid_3d test_refinement_3d

.PHONY: all

vpath %.f90 ../tests

all: $(PROGS_1D) $(PROGS_2D) $(PROGS_3D)

clean:
	$(RM) $(PROGS_2D) $(PROGS_3D) *.mod *.o

$(PROGS_1D): m_octree_mg_1d.o
$(PROGS_2D): m_octree_mg_2d.o
$(PROGS_3D): m_octree_mg_3d.o

%.o: %.f90
	$(F90) -c -o $@ $^ $(F90FLAGS)

%_1d: %.f90
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=1 -DSINGLE_MODULE

%_2d: %.f90
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=2 -DSINGLE_MODULE

%_3d: %.f90
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=3 -DSINGLE_MODULE
