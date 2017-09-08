module m_octree_mg
  use m_morton

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: no_box = 0

  integer, parameter :: nb_offset_2d(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])

  !> Lists of blocks per refinement level
  type level_t
     integer, allocatable :: leaves(:)
     integer, allocatable :: parents(:)
     integer, allocatable :: ids(:)
  end type level_t

  !> Block data structure
  type block_2d_t
     integer               :: lvl
     integer               :: ix(2)
     integer               :: rank
     integer               :: children(4) = no_box
     integer               :: neighbors(4)
     integer               :: parent = no_box
     integer(int64)        :: morton
     real(dp), allocatable :: phi(:, :)
     real(dp), allocatable :: rhs(:, :)
     real(dp), allocatable :: tmp(:, :)
  end type block_2d_t

  type tree_2d
     integer                       :: block_size(2) = -1
     integer                       :: highest_lvl   = 0
     integer                       :: n_boxes
     type(level_t), allocatable    :: lvls(:)
     type(block_2d_t), allocatable :: boxes(:)
  end type tree_2d

  ! Public methods
  public :: create_tree_2d

  ! Public types
  public :: tree_2d

contains

  subroutine create_tree_2d(n_leaves, blvls, branks, bixs, tree)
    integer, intent(in) :: n_leaves
    integer, intent(in) :: blvls(n_leaves)
    integer, intent(in) :: branks(n_leaves)
    integer, intent(in) :: bixs(2, n_leaves)
    type(tree_2d), intent(inout) :: tree

    type(tree_2d) :: copy
    integer, allocatable :: ids(:), ix_sort(:)

    integer(int64) :: mn
    integer(int64), allocatable :: mortons(:)
    integer :: i_ch, ix(2), p_id, nb, nb_id, nb_rev
    integer :: n_boxes, n_done, search_val
    integer :: lvl, highest_lvl, i, n, id, i_box
    integer, allocatable :: leaf_cnt(:), parent_cnt(:), total_cnt(:)

    highest_lvl = maxval(blvls)

    tree%highest_lvl = highest_lvl
    allocate(tree%lvls(highest_lvl))
    allocate(leaf_cnt(highest_lvl))
    allocate(parent_cnt(highest_lvl))
    allocate(total_cnt(highest_lvl))
    leaf_cnt(:)  = 0
    parent_cnt(:) = 0
    total_cnt(:) = 0

    ! Create tree structure

    ! First determine number of leaves and parents per level
    do n = 1, n_leaves
       lvl = blvls(n)
       leaf_cnt(lvl) = leaf_cnt(lvl) + 1
    end do

    do lvl = highest_lvl, 2, -1
       parent_cnt(lvl-1) = (leaf_cnt(lvl) + parent_cnt(lvl))/4
    end do

    ! Now allocate storage per level
    n_boxes = 0
    do lvl = 1, highest_lvl
       allocate(tree%lvls(lvl)%leaves(leaf_cnt(lvl)))
       allocate(tree%lvls(lvl)%parents(parent_cnt(lvl)))
       n        = parent_cnt(lvl) + leaf_cnt(lvl)
       allocate(tree%lvls(lvl)%ids(n))
       n_boxes = n_boxes + n
    end do

    tree%n_boxes = n_boxes
    allocate(tree%boxes(n_boxes))

    ! First add all leaves
    total_cnt(:)  = 0
    i_box = 0

    do n = 1, n_leaves
       i_box = i_box + 1
       lvl   = blvls(n)

       total_cnt(lvl) = total_cnt(lvl) + 1
       tree%lvls(lvl)%ids(total_cnt(lvl)) = i_box

       tree%boxes(i_box)%lvl  = blvls(n)
       tree%boxes(i_box)%ix   = bixs(:, n)
       tree%boxes(i_box)%rank = branks(n)
       tree%boxes(i_box)%morton = &
            morton_from_ix2(tree%boxes(i_box)%ix-1)
    end do

    ! Then add all the parents
    do lvl = highest_lvl, 2, -1
       do n = 1, total_cnt(lvl)
          id = tree%lvls(lvl)%ids(n)

          if (first_child(tree%boxes(id)%ix)) then
             ! Add parent
             i_box = i_box + 1

             total_cnt(lvl-1) = total_cnt(lvl-1) + 1
             tree%lvls(lvl-1)%ids(total_cnt(lvl-1)) = i_box

             tree%boxes(i_box)%lvl = lvl-1
             tree%boxes(i_box)%ix = (tree%boxes(id)%ix + 1) / 2
             tree%boxes(i_box)%rank = tree%boxes(id)%rank
             tree%boxes(i_box)%morton = &
                  morton_from_ix2(tree%boxes(i_box)%ix-1)
          end if
       end do
    end do

    copy = tree
    n_done = 0

    ! Sort all levels according to their Morton order
    do lvl = 1, highest_lvl
       n = total_cnt(lvl)
       ids = copy%lvls(lvl)%ids
       ix_sort = ids

       call morton_rank(copy%boxes(ids)%morton, ix_sort)

       tree%boxes(n_done+1:n_done+n) = copy%boxes(ids(ix_sort))
       tree%lvls(lvl)%ids = [(n_done+i, i=1,n)]

       n_done = n_done + n
    end do

    ! Set parents / children
    do lvl = 2, highest_lvl
       mortons = tree%boxes(tree%lvls(lvl-1)%ids)%morton

       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          ix = tree%boxes(id)%ix

          ! Morton number of parent
          mn = morton_from_ix2((ix-1)/2)

          ! Find parent in morton list
          p_id = tree%lvls(lvl-1)%ids(morton_bsearch(mortons, mn))
          tree%boxes(id)%parent = p_id

          ! Set child on parent
          i_ch = ix2_to_ichild(ix)
          tree%boxes(p_id)%children(i_ch) = id
       end do
    end do

    ! Set neighbors
    do lvl = 1, highest_lvl
       mortons = tree%boxes(tree%lvls(lvl)%ids)%morton

       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          ix = tree%boxes(id)%ix

          do nb = 1, 4
             mn = morton_from_ix2(ix + nb_offset_2d(:, nb) - 1)
             search_val = morton_bsearch(mortons, mn)

             if (search_val /= -1) then
                nb_id = tree%lvls(lvl)%ids(search_val)
                tree%boxes(id)%neighbors(nb) = nb_id
                ! Reverse neighbor direction (odd -> +1, even -> -1)
                nb_rev = nb - 1 + 2 * iand(nb, 1)
                tree%boxes(nb_id)%neighbors(nb_rev) = id
             else
                tree%boxes(id)%neighbors(nb) = no_box
             end if
          end do
       end do
    end do

  end subroutine create_tree_2d

  pure logical function first_child(ix)
    integer, intent(in) :: ix(:)
    first_child = all(iand(ix, 1) == 1)
  end function first_child

  !> Compute the 'child index' for a box with spatial index ix. With 'child
  !> index' we mean the index in the children(:) array of its parent.
  pure integer function ix2_to_ichild(ix)
    integer, intent(in) :: ix(2) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
    ix2_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
  end function ix2_to_ichild

end module m_octree_mg
