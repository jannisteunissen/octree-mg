#include "cpp_macros.h"
module m_build_tree
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: mg_build_rectangle
  public :: mg_add_children
  public :: mg_set_leaves_parents
  public :: mg_set_next_level_ids
  public :: mg_set_refinement_boundaries
  public :: mg_set_neighbors_lvl

contains

  subroutine mg_build_rectangle(mg, domain_size, box_size, dx, r_min, &
       periodic, n_finer)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: domain_size(NDIM)
    integer, intent(in)       :: box_size
    real(dp), intent(in)      :: dx(NDIM)
    real(dp), intent(in)      :: r_min(NDIM)
    logical, intent(in)       :: periodic(NDIM)
    integer, intent(in)       :: n_finer
    integer                   :: IJK, lvl, n, id, nx(NDIM)
    integer                   :: boxes_per_dim(NDIM, lvl_lo_bnd:1)
    integer                   :: periodic_offset(NDIM)

    if (modulo(box_size, 2) /= 0) &
         error stop "box_size should be even"
    if (any(modulo(domain_size, box_size) /= 0)) &
         error stop "box_size does not divide domain_size"

    if (all(periodic)) then
       mg%subtract_mean = .true.
    end if

    nx                  = domain_size
    mg%box_size         = box_size
    mg%box_size_lvl(1)  = box_size
    mg%first_normal_lvl = 1
    mg%dr(:, 1)         = dx
    boxes_per_dim(:, :) = 0
    boxes_per_dim(:, 1) = domain_size / box_size

    do lvl = 1, lvl_lo_bnd+1, -1
       if (any(modulo(nx, 2) == 1 .or. nx == mg%coarsest_grid)) exit

       if (all(modulo(nx/mg%box_size_lvl(lvl), 2) == 0)) then
          mg%box_size_lvl(lvl-1) = mg%box_size_lvl(lvl)
          boxes_per_dim(:, lvl-1) = boxes_per_dim(:, lvl)/2
          mg%first_normal_lvl = lvl-1
       else
          mg%box_size_lvl(lvl-1) = mg%box_size_lvl(lvl)/2
          boxes_per_dim(:, lvl-1) = boxes_per_dim(:, lvl)
       end if

       mg%dr(:, lvl-1) = mg%dr(:, lvl) * 2
       nx = nx / 2
    end do

    mg%lowest_lvl = lvl
    mg%highest_lvl = 1

    do lvl = 2, lvl_hi_bnd
       mg%dr(:, lvl) = mg%dr(:, lvl-1) * 0.5_dp
       mg%box_size_lvl(lvl) = box_size
    end do

    n = sum(product(boxes_per_dim, dim=1)) + n_finer
    allocate(mg%boxes(n))

    ! Create lowest level
    nx = boxes_per_dim(:, mg%lowest_lvl)
#if NDIM == 2
    periodic_offset = [nx(1)-1, (nx(2)-1)*nx(1)]
#elif NDIM == 3
    periodic_offset = [nx(1)-1, (nx(2)-1)*nx(1), &
         (nx(3)-1) * nx(2) * nx(1)]
#endif

    do KJI_DO_VEC(nx)
       mg%n_boxes = mg%n_boxes + 1
       n          = mg%n_boxes

       mg%boxes(n)%rank        = 0
       mg%boxes(n)%id          = n
       mg%boxes(n)%lvl         = mg%lowest_lvl
       mg%boxes(n)%ix(:)       = [IJK]
       mg%boxes(n)%r_min(:)    = r_min + (mg%boxes(n)%ix(:) - 1) * &
            mg%box_size_lvl(mg%lowest_lvl) * mg%dr(:, mg%lowest_lvl)
       mg%boxes(n)%parent      = no_box
       mg%boxes(n)%children(:) = no_box

       ! Set default neighbors
#if NDIM == 2
       mg%boxes(n)%neighbors(:) = [n-1, n+1, n-nx(1), n+nx(1)]
#elif NDIM == 3
       mg%boxes(n)%neighbors(:) = [n-1, n+1, n-nx(1), n+nx(1), &
            n-nx(1)*nx(2), n+nx(1)*nx(2)]
#endif

       ! Handle boundaries
       where ([IJK] == 1 .and. .not. periodic)
          mg%boxes(n)%neighbors(1:num_neighbors:2) = &
               physical_boundary
       end where
       where ([IJK] == 1 .and. periodic)
          mg%boxes(n)%neighbors(1:num_neighbors:2) = &
               n + periodic_offset
       end where

       where ([IJK] == nx .and. .not. periodic)
          mg%boxes(n)%neighbors(2:num_neighbors:2) = &
               physical_boundary
       end where
       where ([IJK] == nx .and. periodic)
          mg%boxes(n)%neighbors(2:num_neighbors:2) = &
               n - periodic_offset
       end where
    end do; CLOSE_DO

    mg%lvls(mg%lowest_lvl)%ids = [(n, n=1, mg%n_boxes)]

    ! Add higher levels
    do lvl = mg%lowest_lvl, 0
       if (mg%box_size_lvl(lvl+1) == mg%box_size_lvl(lvl)) then
          do i = 1, size(mg%lvls(lvl)%ids)
             id = mg%lvls(lvl)%ids(i)
             call mg_add_children(mg, id)
          end do

          call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))
          call mg_set_next_level_ids(mg, lvl)
          call mg_set_neighbors_lvl(mg, lvl+1)
       else
          do i = 1, size(mg%lvls(lvl)%ids)
             id = mg%lvls(lvl)%ids(i)
             call add_single_child(mg, id, size(mg%lvls(lvl)%ids))
          end do

          call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))
          call mg_set_next_level_ids(mg, lvl)
       end if
    end do

    call mg_set_leaves_parents(mg%boxes, mg%lvls(1))

    ! No refinement boundaries
    do lvl = mg%lowest_lvl, 1
       if (allocated(mg%lvls(lvl)%ref_bnds)) &
            deallocate(mg%lvls(lvl)%ref_bnds)
       allocate(mg%lvls(lvl)%ref_bnds(0))
    end do

    mg%tree_created = .true.
  end subroutine mg_build_rectangle

  subroutine mg_set_neighbors_lvl(mg, lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id

    do i = 1, size(mg%lvls(lvl)%ids)
       id = mg%lvls(lvl)%ids(i)
       call set_neighbs(mg%boxes, id)
    end do
  end subroutine mg_set_neighbors_lvl

  subroutine mg_set_next_level_ids(mg, lvl)
    type(mg_t), intent(inout)  :: mg
    integer, intent(in)        :: lvl
    integer                    :: n, i, id

    if (allocated(mg%lvls(lvl+1)%ids)) &
         deallocate(mg%lvls(lvl+1)%ids)

    ! Set next level ids to children of this level
    if (mg%box_size_lvl(lvl+1) == mg%box_size_lvl(lvl)) then
       n = num_children * size(mg%lvls(lvl)%parents)
       allocate(mg%lvls(lvl+1)%ids(n))

       n = num_children
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)
          mg%lvls(lvl+1)%ids(n*(i-1)+1:n*i) = mg%boxes(id)%children
       end do
    else
       n = size(mg%lvls(lvl)%parents)
       allocate(mg%lvls(lvl+1)%ids(n))

       n = 1
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)
          mg%lvls(lvl+1)%ids(i) = mg%boxes(id)%children(1)
       end do
    end if

  end subroutine mg_set_next_level_ids

  ! Set the neighbors of id (using their parent)
  subroutine set_neighbs(boxes, id)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: nb, nb_id

    do nb = 1, num_neighbors
       if (boxes(id)%neighbors(nb) == no_box) then
          nb_id = find_neighb(boxes, id, nb)
          if (nb_id > no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(neighb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_neighbs

  !> Get the id of neighbor nb of boxes(id), through its parent
  function find_neighb(boxes, id, nb) result(nb_id)
    type(box_t), intent(in) :: boxes(:) !< List with all the boxes
    integer, intent(in)      :: id       !< Box whose neighbor we are looking for
    integer, intent(in)      :: nb       !< Neighbor index
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id    = boxes(id)%parent
    old_pid = p_id
    c_ix    = ix_to_ichild(boxes(id)%ix)
    d       = neighb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (child_low(d, c_ix) .eqv. neighb_low(nb)) then
       p_id = boxes(p_id)%neighbors(nb)
    end if

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(child_rev(c_ix, d))
  end function find_neighb

  !> Create a list of leaves and a list of parents for a level
  subroutine mg_set_leaves_parents(boxes, level)
    type(box_t), intent(in)   :: boxes(:) !< List of boxes
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
  end subroutine mg_set_leaves_parents

  !> Create a list of refinement boundaries (from the coarse side)
  subroutine mg_set_refinement_boundaries(boxes, level)
    type(box_t), intent(in)    :: boxes(:)
    type(lvl_t), intent(inout) :: level
    integer, allocatable       :: tmp(:)
    integer                    :: i, id, nb, nb_id, ix

    if (allocated(level%ref_bnds)) deallocate(level%ref_bnds)

    if (size(level%parents) == 0) then
       ! There are no refinement boundaries
       allocate(level%ref_bnds(0))
    else
       allocate(tmp(size(level%leaves)))
       ix = 0
       do i = 1, size(level%leaves)
          id = level%leaves(i)

          do nb = 1, num_neighbors
             nb_id = boxes(id)%neighbors(nb)
             if (nb_id > no_box) then
                if (has_children(boxes(nb_id))) then
                   ix = ix + 1
                   tmp(ix) = id
                   exit
                end if
             end if
          end do
       end do

       allocate(level%ref_bnds(ix))
       level%ref_bnds(:) = tmp(1:ix)
    end if
  end subroutine mg_set_refinement_boundaries

  subroutine mg_add_children(mg, id)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)          :: id      !< Id of box that gets children
    integer                      :: lvl, i, nb, child_nb(2**(NDIM-1))
    integer                      :: c_ids(num_children), c_id, c_ix_base(NDIM)

    c_ids                 = [(mg%n_boxes+i, i=1,num_children)]
    mg%n_boxes            = mg%n_boxes + num_children
    mg%boxes(id)%children = c_ids
    c_ix_base             = 2 * mg%boxes(id)%ix - 1
    lvl                   = mg%boxes(id)%lvl+1

    do i = 1, num_children
       c_id                     = c_ids(i)
       mg%boxes(c_id)%rank      = mg%boxes(id)%rank
       mg%boxes(c_id)%ix        = c_ix_base + child_dix(:, i)
       mg%boxes(c_id)%lvl       = lvl
       mg%boxes(c_id)%parent    = id
       mg%boxes(c_id)%children  = no_box
       mg%boxes(c_id)%neighbors = no_box
       mg%boxes(c_id)%r_min     = mg%boxes(id)%r_min + &
            mg%dr(:, lvl) * child_dix(:, i) * mg%box_size
    end do

    ! Set boundary conditions at children
    do nb = 1, num_neighbors
       if (mg%boxes(id)%neighbors(nb) < no_box) then
          child_nb = c_ids(child_adj_nb(:, nb)) ! Neighboring children
          mg%boxes(child_nb)%neighbors(nb) = mg%boxes(id)%neighbors(nb)
       end if
    end do
  end subroutine mg_add_children

  subroutine add_single_child(mg, id, n_boxes_lvl)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id !< Id of box that gets children
    integer, intent(in)       :: n_boxes_lvl
    integer                   :: lvl, c_id

    c_id                     = mg%n_boxes + 1
    mg%n_boxes               = mg%n_boxes + 1
    mg%boxes(id)%children(1) = c_id
    lvl                      = mg%boxes(id)%lvl+1

    mg%boxes(c_id)%rank      = mg%boxes(id)%rank
    mg%boxes(c_id)%ix        = mg%boxes(id)%ix
    mg%boxes(c_id)%lvl       = lvl
    mg%boxes(c_id)%parent    = id
    mg%boxes(c_id)%children  = no_box
    where (mg%boxes(id)%neighbors == physical_boundary)
       mg%boxes(c_id)%neighbors = mg%boxes(id)%neighbors
    elsewhere
       mg%boxes(c_id)%neighbors = mg%boxes(id)%neighbors + n_boxes_lvl
    end where
    mg%boxes(c_id)%r_min = mg%boxes(id)%r_min

  end subroutine add_single_child

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

end module m_build_tree
