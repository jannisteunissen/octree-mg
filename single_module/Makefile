F90 := mpif90
PROGS_3D := test_uniform_grid_3d

.PHONY: all

all: $(PROGS_3D)

clean:
	$(RM) $(PROGS_3D) *.mod *.o

test_uniform_grid_3d: m_octree_mg.o

%.o: %.f90
	$(F90) -c -o $@ $^ -cpp -DNDIM=3 -Wall

%_3d: %.f90
	$(F90) -o $@ $^ -cpp -DNDIM=3 -Wall
