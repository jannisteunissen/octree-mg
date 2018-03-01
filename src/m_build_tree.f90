module m_build_tree
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: build_uniform_tree

contains

  subroutine build_uniform_tree(mg, block_size, domain_size, dr_coarse)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in) :: block_size
    integer, intent(in) :: domain_size
    real(dp), intent(in) :: dr_coarse

    integer              :: i, id, lvl, n
    integer              :: max_lvl, n_boxes
    logical, allocatable :: has_children(:)

    max_lvl = 1 + nint(log(real(domain_size, dp)/block_size) / log(2.0_dp))

    if (block_size * 2**(max_lvl-1) /= domain_size) then
       error stop "Invalid domain_size should be 2^N * block_size"
    end if

    mg%box_size = block_size
    mg%highest_lvl = max_lvl

    n_boxes = 0
    do lvl = 1, max_lvl
       n_boxes = n_boxes + 2**((lvl-1) * 2)
    end do
    allocate(mg%boxes(n_boxes))

    mg%dr = [(dr_coarse * 0.5**(n-1), n = 1, max_lvl)]

    mg%n_boxes               = 1
    mg%boxes(1)%rank         = 0
    mg%boxes(1)%id           = 1
    mg%boxes(1)%lvl          = 0
    mg%boxes(1)%ix           = [1, 1]
    mg%boxes(1)%r_min        = [0.0_dp, 0.0_dp]
    mg%boxes(1)%parent       = no_box
    ! TODO
    mg%boxes(1)%neighbors(:) = physical_boundary
    mg%boxes(1)%children(:)  = no_box

    ! Boxes on finest level
    n = 2**((max_lvl-1) * 2)
    allocate(has_children(n_boxes))
    has_children(1:n_boxes-n) = .true.
    has_children(n_boxes-n+1:n_boxes) = .false.

    allocate(mg%lvls(max_lvl))
    allocate(mg%lvls(1)%ids(1))
    mg%lvls(1)%ids(1) = 1

    do lvl = 1, max_lvl-1
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          call add_children(mg, id)
       end do

       call set_leaves_parents(mg%boxes, mg%lvls(lvl))

       ! Set next level ids to children of this level
       n = 4 * size(mg%lvls(lvl)%parents)
       allocate(mg%lvls(lvl+1)%ids(n))

       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)
          mg%lvls(lvl+1)%ids(4*i-3:4*i) = mg%boxes(id)%children
       end do
    end do

    call set_leaves_parents(mg%boxes, mg%lvls(max_lvl))

    ! Update connectivity of new boxes
    do lvl = 2, max_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          call set_neighbs_2d(mg%boxes, id)
       end do
    end do

  end subroutine build_uniform_tree

  ! Set the neighbors of id (using their parent)
  subroutine set_neighbs_2d(boxes, id)
    type(box_2d_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: nb, nb_id

    do nb = 1, num_neighbors
       if (boxes(id)%neighbors(nb) == no_box) then
          nb_id = find_neighb_2d(boxes, id, nb)
          if (nb_id > no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(neighb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_neighbs_2d

  !> Get the id of neighbor nb of boxes(id), through its parent
  function find_neighb_2d(boxes, id, nb) result(nb_id)
    type(box_2d_t), intent(in) :: boxes(:) !< List with all the boxes
    integer, intent(in)      :: id       !< Box whose neighbor we are looking for
    integer, intent(in)      :: nb       !< Neighbor index
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id    = boxes(id)%parent
    old_pid = p_id
    c_ix    = ix_to_ichild(boxes(id)%ix)
    d       = neighb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (child_low(c_ix, d) .eqv. neighb_low(nb)) &
         p_id = boxes(p_id)%neighbors(nb)

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(child_rev(c_ix, d))
  end function find_neighb_2d

  !> Return .true. if a box has children
  elemental logical function has_children(box)
    type(box_2d_t), intent(in) :: box

    ! Boxes are either fully refined or not, so we only need to check one of the
    ! children
    has_children = (box%children(1) /= no_box)
  end function has_children

  !> Create a list of leaves and a list of parents for a level
  subroutine set_leaves_parents(boxes, level)
    type(box_2d_t), intent(in)   :: boxes(:) !< List of boxes
    type(lvl_t), intent(inout) :: level !< Level type which contains the indices of boxes
    integer                    :: i, id, i_leaf, i_parent
    integer                    :: n_parents, n_leaves

    n_parents = count(has_children(boxes(level%ids)))
    n_leaves = size(level%ids) - n_parents

    if (.not. allocated(level%parents)) then
       allocate(level%parents(n_parents))
    else if (n_parents /= size(level%parents)) then
       deallocate(level%parents)
       allocate(level%parents(n_parents))
    end if

    if (.not. allocated(level%leaves)) then
       allocate(level%leaves(n_leaves))
    else if (n_leaves /= size(level%leaves)) then
       deallocate(level%leaves)
       allocate(level%leaves(n_leaves))
    end if

    i_leaf   = 0
    i_parent = 0
    do i = 1, size(level%ids)
       id = level%ids(i)
       if (has_children(boxes(id))) then
          i_parent                = i_parent + 1
          level%parents(i_parent) = id
       else
          i_leaf               = i_leaf + 1
          level%leaves(i_leaf) = id
       end if
    end do
  end subroutine set_leaves_parents

  subroutine add_children(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id      !< Id of box that gets children
    integer                      :: lvl, i, nb, child_nb(2**(2-1))
    integer                      :: c_ids(4), c_id, c_ix_base(2)

    c_ids                 = [(mg%n_boxes+i, i=1,4)]
    mg%n_boxes            = mg%n_boxes + 4
    mg%boxes(id)%children = c_ids
    c_ix_base             = 2 * mg%boxes(id)%ix - 1
    lvl                   = mg%boxes(id)%lvl+1

    do i = 1, 4
       c_id                     = c_ids(i)
       mg%boxes(c_id)%rank      = mg%boxes(id)%rank
       mg%boxes(c_id)%ix        = c_ix_base + child_dix(:, i)
       mg%boxes(c_id)%lvl       = lvl
       mg%boxes(c_id)%parent    = id
       mg%boxes(c_id)%children  = no_box
       mg%boxes(c_id)%neighbors = no_box
       mg%boxes(c_id)%r_min     = mg%boxes(id)%r_min + &
            0.5_dp * mg%dr(lvl) * child_dix(:, i) * mg%box_size
    end do

    ! Set boundary conditions at children
    do nb = 1, 4
       if (mg%boxes(id)%neighbors(nb) < no_box) then
          child_nb = c_ids(child_adj_nb(:, nb)) ! Neighboring children
          mg%boxes(child_nb)%neighbors(nb) = mg%boxes(id)%neighbors(nb)
       end if
    end do
  end subroutine add_children

  ! subroutine build_tree_2d(n_leaves, blvls, branks, bixs, tree)
  !   integer, intent(in) :: n_leaves
  !   integer, intent(in) :: blvls(n_leaves)
  !   integer, intent(in) :: branks(n_leaves)
  !   integer, intent(in) :: bixs(2, n_leaves)
  !   type(tree_2d), intent(inout) :: tree

  !   type(tree_2d) :: copy
  !   integer, allocatable :: ids(:), ix_sort(:)

  !   integer(int64) :: mn
  !   integer(int64), allocatable :: mortons(:)
  !   integer :: i_ch, ix(2), p_id, nb, nb_id, nb_rev
  !   integer :: n_boxes, n_done, search_val
  !   integer :: lvl, highest_lvl, i, n, id, i_box
  !   integer :: child_ranks(4)
  !   integer, allocatable :: leaf_cnt(:), parent_cnt(:), total_cnt(:)

  !   highest_lvl = maxval(blvls)

  !   tree%highest_lvl = highest_lvl
  !   allocate(tree%lvls(highest_lvl))
  !   allocate(leaf_cnt(highest_lvl))
  !   allocate(parent_cnt(highest_lvl))
  !   allocate(total_cnt(highest_lvl))
  !   leaf_cnt(:)  = 0
  !   parent_cnt(:) = 0
  !   total_cnt(:) = 0

  !   ! Create tree structure

  !   ! First determine number of leaves and parents per level
  !   do n = 1, n_leaves
  !      lvl = blvls(n)
  !      leaf_cnt(lvl) = leaf_cnt(lvl) + 1
  !   end do

  !   do lvl = highest_lvl, 2, -1
  !      parent_cnt(lvl-1) = (leaf_cnt(lvl) + parent_cnt(lvl))/4
  !   end do

  !   ! Now allocate storage per level
  !   n_boxes = 0
  !   do lvl = 1, highest_lvl
  !      allocate(tree%lvls(lvl)%leaves(leaf_cnt(lvl)))
  !      allocate(tree%lvls(lvl)%parents(parent_cnt(lvl)))
  !      n        = parent_cnt(lvl) + leaf_cnt(lvl)
  !      allocate(tree%lvls(lvl)%ids(n))
  !      n_boxes = n_boxes + n
  !   end do

  !   tree%n_boxes = n_boxes
  !   allocate(tree%boxes(n_boxes))

  !   ! First add all leaves
  !   total_cnt(:)  = 0
  !   i_box = 0

  !   do n = 1, n_leaves
  !      i_box = i_box + 1
  !      lvl   = blvls(n)

  !      total_cnt(lvl) = total_cnt(lvl) + 1
  !      tree%lvls(lvl)%ids(total_cnt(lvl)) = i_box

  !      tree%boxes(i_box)%lvl  = blvls(n)
  !      tree%boxes(i_box)%ix   = bixs(:, n)
  !      tree%boxes(i_box)%rank = branks(n)
  !      tree%boxes(i_box)%morton = &
  !           morton_from_ix2(tree%boxes(i_box)%ix-1)
  !   end do

  !   ! Then add all the parents
  !   do lvl = highest_lvl, 2, -1
  !      do n = 1, total_cnt(lvl)
  !         id = tree%lvls(lvl)%ids(n)

  !         if (first_child(tree%boxes(id)%ix)) then
  !            ! Add parent
  !            i_box = i_box + 1

  !            total_cnt(lvl-1) = total_cnt(lvl-1) + 1
  !            tree%lvls(lvl-1)%ids(total_cnt(lvl-1)) = i_box

  !            tree%boxes(i_box)%lvl = lvl-1
  !            tree%boxes(i_box)%ix = (tree%boxes(id)%ix + 1) / 2
  !            tree%boxes(i_box)%morton = &
  !                 morton_from_ix2(tree%boxes(i_box)%ix-1)
  !         end if
  !      end do
  !   end do

  !   copy = tree
  !   n_done = 0

  !   ! Sort all levels according to their Morton order
  !   do lvl = 1, highest_lvl
  !      n = total_cnt(lvl)
  !      ids = copy%lvls(lvl)%ids
  !      ix_sort = ids

  !      call morton_rank(copy%boxes(ids)%morton, ix_sort)

  !      tree%boxes(n_done+1:n_done+n) = copy%boxes(ids(ix_sort))
  !      tree%lvls(lvl)%ids = [(n_done+i, i=1,n)]

  !      n_done = n_done + n
  !   end do

  !   ! Set parents / children
  !   do lvl = 2, highest_lvl
  !      mortons = tree%boxes(tree%lvls(lvl-1)%ids)%morton

  !      do i = 1, size(tree%lvls(lvl)%ids)
  !         id = tree%lvls(lvl)%ids(i)
  !         ix = tree%boxes(id)%ix

  !         ! Morton number of parent
  !         mn = morton_from_ix2((ix-1)/2)

  !         ! Find parent in morton list
  !         p_id = tree%lvls(lvl-1)%ids(morton_bsearch(mortons, mn))
  !         tree%boxes(id)%parent = p_id

  !         ! Set child on parent
  !         i_ch = ix2_to_ichild(ix)
  !         tree%boxes(p_id)%children(i_ch) = id
  !      end do
  !   end do

  !   ! Set neighbors (can optimize this later)
  !   do lvl = 1, highest_lvl
  !      mortons = tree%boxes(tree%lvls(lvl)%ids)%morton

  !      do i = 1, size(tree%lvls(lvl)%ids)
  !         id = tree%lvls(lvl)%ids(i)
  !         ix = tree%boxes(id)%ix

  !         do nb = 1, 4
  !            mn = morton_from_ix2(ix + nb_offset_2d(:, nb) - 1)
  !            search_val = morton_bsearch(mortons, mn)

  !            if (search_val /= -1) then
  !               nb_id = tree%lvls(lvl)%ids(search_val)
  !               tree%boxes(id)%neighbors(nb) = nb_id
  !               ! Reverse neighbor direction (odd -> +1, even -> -1)
  !               nb_rev = nb - 1 + 2 * iand(nb, 1)
  !               tree%boxes(nb_id)%neighbors(nb_rev) = id
  !            else
  !               tree%boxes(id)%neighbors(nb) = no_box
  !            end if
  !         end do
  !      end do
  !   end do

  !   ! Fill arrays of parents/leaves, and set ranks of parents
  !   leaf_cnt(:)  = 0
  !   parent_cnt(:) = 0

  !   do lvl = highest_lvl, 1, -1
  !      do i = 1, size(tree%lvls(lvl)%ids)
  !         id = tree%lvls(lvl)%ids(i)
  !         if (tree%boxes(id)%children(1) > no_box) then
  !            parent_cnt(lvl) = parent_cnt(lvl) + 1
  !            tree%lvls(lvl)%parents(parent_cnt(lvl)) = id

  !            child_ranks = tree%boxes(tree%boxes(id)%children)%rank
  !            tree%boxes(id)%rank = most_popular(child_ranks)
  !         else
  !            leaf_cnt(lvl) = leaf_cnt(lvl) + 1
  !            tree%lvls(lvl)%leaves(leaf_cnt(lvl)) = id
  !         end if
  !      end do
  !   end do

  ! end subroutine build_tree_2d

  ! pure logical function first_child(ix)
  !   integer, intent(in) :: ix(:)
  !   first_child = all(iand(ix, 1) == 1)
  ! end function first_child

  ! pure integer function most_popular(list)
  !   integer, intent(in) :: list(:)
  !   integer             :: i, best_count, current_count

  !   best_count   = 0
  !   most_popular = -1

  !   do i = 1, size(list)
  !      current_count = count(list == list(i))

  !      if (current_count > best_count) then
  !         most_popular = list(i)
  !      end if
  !   end do

  ! end function most_popular

  ! !> Compute the 'child index' for a box with spatial index ix. With 'child
  ! !> index' we mean the index in the children(:) array of its parent.
  ! pure integer function ix2_to_ichild(ix)
  !   integer, intent(in) :: ix(2) !< Spatial index of the box
  !   ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
  !   ix2_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
  ! end function ix2_to_ichild

end module m_build_tree
