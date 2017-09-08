octree-mg
==

Geometric multigrid framework for octree/quadtree grids, using MPI. Focus on simple and efficient methods:

* 2nd order operators with sparse stencil (5/7-point)
* Gauss-Seidel red-black
* Linear prolongation

## Design ideas

What type of input to create tree? Store rank and block index for neighbor,
perhaps separate integer indicating kind of neighbor

Specified beforehand:
* Block size
* Operators
* Boundary conditions

Create octree-mg data structure:
* For each block: refinement level, spatial index, mpi rank
* Parents assigned to rank with most children

Parallel solve:
* rhs for each block (no ghost cells)
* guess for phi
* Output: phi without ghost cells

