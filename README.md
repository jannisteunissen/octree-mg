Octree-mg
==

The octree-mg library implements parallel geometric multigrid methods on
quadtree/octree grids, which can be used to solve elliptic PDEs such as
Poissons's equation. The provided solvers can be used in existing
adaptive-mesh-refinement (AMR) frameworks that employ quadtree/octree grids.

## Usage

Type `make` in the top folder, and run the programs in the `tests` folder.

## Including the library in other projects

You can either include a full copy of `octree-mg` using e.g. git submodules of subtrees, or include it as a single Fortran module, see the folder `single_module`.

## Features

* MPI parallelization
* The same code can be used in 2D/3D (compiled into two different libraries)
* Support for adaptively refined grids, with a consistent discretization near refinement boundaries
* Rectangular grids can easily be created (e.g. 512 x 256 x 256 cells).
* Operators with sparse stencils (5/7-point in 2D/3D) are supported
* Support for periodic, Dirichlet, Neumann, and continuous boundary conditions
* Coarse grids are automatically created

## Restrictions

* One layer of ghost cells is used, and diagonal ghost cells are currently not set. This
  means stencils with diagonal elements are not possible.
* Point-wise smoothers are employed, so currently a requirement is that dx, dy and dz are similar

## Using the free-space solver in MPI-AMRVAC

The "single-module" version of this library is included in [MPI-AMRVAC](http://amrvac.org/), but it does not include the free-space module. An example of how to do this can be found under [tests/multigrid/free_space_3d](https://github.com/amrvac/amrvac/tree/master/tests/multigrid/free_space_3d) in the MPI-AMRVAC folder, see in particular the file `local.make`.
