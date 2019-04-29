F90 := mpif90
F90FLAGS := -O2 -Wall -cpp
PROGS_2D := test_uniform_grid_2d
PROGS_3D := test_uniform_grid_3d

.PHONY: all

all: $(PROGS_2D) $(PROGS_3D)

clean:
	$(RM) $(PROGS_2D) $(PROGS_3D) *.mod *.o

$(PROGS_2D): m_octree_mg_2d.o
$(PROGS_3D): m_octree_mg_3d.o

%.o: %.f90
	$(F90) -c -o $@ $^ $(F90FLAGS)

%_2d: %.f90
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=2

%_3d: %.f90
	$(F90) -o $@ $^ $(F90FLAGS) -DNDIM=3
