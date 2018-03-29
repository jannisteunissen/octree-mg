PROGS := test_uniform_grid_2d test_uniform_grid_3d
LIBS := omg

include ../makerules.make

.PHONY: all

all: $(PROGS)

clean:
	$(RM) $(PROGS)

# Special rules to compile 2d/3d version of same program
%_2d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

test_uniform_grid_2d: LIBDIRS:=../lib_2d
test_uniform_grid_2d: INCDIRS:=../lib_2d
test_uniform_grid_2d: export NDIM:=2

test_uniform_grid_3d: LIBDIRS:=../lib_3d
test_uniform_grid_3d: INCDIRS:=../lib_3d
test_uniform_grid_3d: export NDIM:=3




