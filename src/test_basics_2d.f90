program test_basics
  use m_octree_mg

  implicit none

  type(tree_2d) :: tree
  integer :: id

  integer, parameter :: n_leaves = 7
  integer            :: lvls(n_leaves)
  integer            :: ranks(n_leaves)
  integer            :: ixs(2, n_leaves)

  lvls(1:3) = 2
  lvls(4:7) = 3
  ranks(1:3) = 0
  ranks(4:7) = 1
  ixs(:, 3) = [2,1]
  ixs(:, 2) = [1,2]
  ixs(:, 1) = [2,2]
  ixs(:, 4) = [1,1]
  ixs(:, 5) = [2,2]
  ixs(:, 6) = [1,2]
  ixs(:, 7) = [2,1]

  print *, "Test basics"

  call create_tree_2d(n_leaves, lvls, ranks, ixs, tree)

  print *, "highest lvl", tree%highest_lvl

  do id = 1, tree%n_boxes
     print *, "lvl", tree%boxes(id)%lvl
     print *, "ix", tree%boxes(id)%ix
     print *, "rank", tree%boxes(id)%rank
     print *, "parent", tree%boxes(id)%parent
     print *, "children", tree%boxes(id)%children
     print *, "neighbors", tree%boxes(id)%neighbors
     print *, ""
  end do

end program test_basics
