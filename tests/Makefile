PROGS_1D := test_uniform_grid_1d test_performance_1d test_refinement_1d
PROGS_2D := test_uniform_grid_2d test_performance_2d test_refinement_2d
PROGS_3D := test_uniform_grid_3d test_performance_3d test_free_space_3d test_refinement_3d
LIBS = omg

include ../makerules.make

.PHONY: all clean run

all: $(PROGS_1D) $(PROGS_2D) $(PROGS_3D)

clean:
	$(RM) $(PROGS_2D) $(PROGS_3D)

# Run with some combinations of settings
run: $(PROGS_2D) $(PROGS_3D)
	mpirun ./test_performance_2d 16 64 64 10 t
	mpirun ./test_performance_3d 16 64 64 64 10 t
	mpirun ./test_performance_2d 16 64 64 10 f
	mpirun ./test_performance_3d 16 64 64 64 10 f
	mpirun ./test_uniform_grid_2d 16 64 64 10 t
	mpirun ./test_uniform_grid_3d 16 64 64 64 10 t
	mpirun ./test_uniform_grid_2d 16 64 64 10 f
	mpirun ./test_uniform_grid_3d 16 64 64 64 10 f
	mpirun ./test_refinement_2d 3 16 64 64 10 t
	mpirun ./test_refinement_3d 3 16 64 64 64 10 t
	mpirun ./test_refinement_2d 3 16 64 64 10 f
	mpirun ./test_refinement_3d 3 16 64 64 64 10 f
	mpirun ./test_free_space_3d 16 64 64 64

$(PROGS_1D): LIBDIRS+=../lib_1d
$(PROGS_1D): INCDIRS+=../lib_1d
$(PROGS_1D): ../lib_1d/libomg.a

$(PROGS_2D): LIBDIRS+=../lib_2d
$(PROGS_2D): INCDIRS+=../lib_2d
$(PROGS_2D): ../lib_2d/libomg.a

$(PROGS_3D): LIBDIRS+=../lib_3d
$(PROGS_3D): INCDIRS+=../lib_3d
$(PROGS_3D): ../lib_3d/libomg.a

# Special rules to compile 2d/3d version of same program
%_1d: %.f90
	$(F90) -o $@ $^ -cpp -DNDIM=1 $(F90FLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_2d: %.f90
	$(F90) -o $@ $^ -cpp -DNDIM=2 $(F90FLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(F90) -o $@ $^ -cpp -DNDIM=3 $(F90FLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
