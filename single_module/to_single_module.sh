#!/bin/bash

# This script generates a single fortran module from all the octree-mg modules,
# which can make it easier to include the library in another project.

# All modules to include
modules="../src/m_data_structures.f90 ../src/m_allocate_storage.f90
../src/m_diffusion.f90 ../src/m_laplacian.f90 ../src/m_vhelmholtz.f90
../src/m_build_tree.f90 ../src/m_load_balance.f90 ../src/m_vlaplacian.f90
../src/m_communication.f90 ../src/m_ghost_cells.f90 ../src/m_prolong.f90
../src/m_helmholtz.f90 ../src/m_multigrid.f90 ../src/m_restrict.f90
../src/m_mrgrnk.f90"

# Name of the new modules
module_xd="m_octree_mg.f90"
module_2d="m_octree_mg_2d.f90"
module_3d="m_octree_mg_3d.f90"
newline=$'\n'
header=""
code=""

# Split existing modules into a header and a code section
for mod in $modules; do
    modname=$"${newline}  !! File ${mod}${newline}"
    # Grep the lines before/after the 'contains' line
    tmp=$(grep -E '^ *contains *$' -B 5000 "${mod}")
    if [ -n "$tmp" ]; then
        header="${header}${modname}${tmp}"
    fi

    tmp=$(grep -E ' *contains *$' -A 5000 "${mod}")
    if [ -n "$tmp" ]; then
        code="${code}${modname}${tmp}"
    fi
done

# Remove certain lines
header=$(echo "$header" | sed 's/^ *implicit none.*$//')
header=$(echo "$header" | sed 's/^ *private *$//')
header=$(echo "$header" | sed 's/^ *public *$//')
header=$(echo "$header" | sed 's/^ *contains.*$//')
header=$(echo "$header" | sed 's/^ *module m_.*$//')
header=$(echo "$header" | sed 's/^ *Module m_.*$//')
header=$(echo "$header" | sed 's/^ *#include .*$//')
header=$(echo "$header" | sed 's/^ *use m_.*$//')

# Merge blank lines
header=$(echo "$header" | cat -s)

# Remove certain lines
code=$(echo "$code" | sed 's/^ *contains.*$//')
code=$(echo "$code" | sed 's/^ *end module.*$//')
code=$(echo "$code" | sed 's/^ *use m_.*$//')

# Merge blank lines
code=$(echo "$code" | cat -s)

# Write to the new module
echo '! Single module generated from the octree-mg sources.
! This file can be easier to include in existing projects.
!
! Notes:
! 1. The module name is here extended by _2d or _3d
! 2. The free space Poisson solver is not included here.
! 3. It is best to make changes in the original repository at
!    https://github.com/jannisteunissen/octree-mg
! 4. This file can be generated as follows:
!    cd octree-mg/single_module
!    ./to_single_module.sh
!
! The modules can be compiled with:
! mpif90 -c m_octree_mg_2d.f90 [other options]
! mpif90 -c m_octree_mg_3d.f90 [other options]
! mpif90 -c m_octree_mg.f90 -cpp -DNDIM=3 [other options]
' > "$module_xd"

# Copy in preprocessor directives
cat "../src/cpp_macros.h" >> "$module_xd"

# Beginning of the module
echo '
#if NDIM == 2
module m_octree_mg_2d
#elif NDIM == 3
module m_octree_mg_3d
#endif
  use mpi
  implicit none
  private' >> "$module_xd"

# Copy the combined header in
echo "$header" >> "$module_xd"

echo "contains" >> "$module_xd"

# Copy the subroutines and functions in
echo "$code" >> "$module_xd"

# Close the module
echo '
#if NDIM == 2
end module m_octree_mg_2d
#elif NDIM == 3
end module m_octree_mg_3d
#endif' >> "$module_xd"
echo "Generated $module_xd"

# Generate 2D and 3D versions
gfortran -E -cpp -P -DNDIM=2 "$module_xd" -o "$module_2d"
echo "Generated $module_2d"

gfortran -E -cpp -P -DNDIM=3 "$module_xd" -o "$module_3d"
echo "Generated $module_3d"
