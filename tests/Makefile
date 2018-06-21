PROGS_2D := test_uniform_grid_2d test_performance_2d
PROGS_3D := test_uniform_grid_3d test_performance_3d test_free_space_3d
LIBS = omg

include ../makerules.make

test_free_space_3d: INCDIRS += $(HOME)/opt/sw/free_space_fft
test_free_space_3d: LIBDIRS += $(HOME)/opt/sw/free_space_fft
test_free_space_3d: LIBS += poisfft

.PHONY: all

all: $(PROGS_2D) $(PROGS_3D)

clean:
	$(RM) $(PROGS_2D) $(PROGS_3D)

# Special rules to compile 2d/3d version of same program
%_2d: %.f90
	$(F90C) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(F90C) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

$(PROGS_2D): LIBDIRS+=../lib_2d
$(PROGS_2D): INCDIRS+=../lib_2d
$(PROGS_2D): export NDIM:=2
$(PROGS_2D): ../lib_2d/libomg.a

$(PROGS_3D): LIBDIRS+=../lib_3d
$(PROGS_3D): INCDIRS+=../lib_3d
$(PROGS_3D): export NDIM:=3
$(PROGS_3D): ../lib_3d/libomg.a




