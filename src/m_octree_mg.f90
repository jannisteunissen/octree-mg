module m_octree_mg
  use m_morton

  implicit none
  private

  ! Public methods
  public :: create_tree_2d

  ! Public types
  public :: tree_2d

contains

  ! subroutine tree_per_rank(tree, my_tree)

  subroutine mg2d_fas_vcycle(tree, mg, set_residual, highest_lvl)
    type(tree_2d_t), intent(inout)    :: tree         !< Tree to do multigrid on
    type(mg2d_t), intent(in)      :: mg           !< Multigrid options
    logical, intent(in)           :: set_residual !< If true, store residual in i_tmp
    integer, intent(in), optional :: highest_lvl  !< Maximum level for V-cycle
    integer                       :: lvl, min_lvl, i, id, max_lvl

    call check_mg(mg)           ! Check whether mg options are set
    min_lvl = lbound(tree%lvls, 1)
    max_lvl = tree%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
       ! correction later
       call update_coarse(tree, lvl, mg)
    end do

    lvl = min_lvl
    call gsrb_boxes(tree, tree%lvls(lvl)%ids, mg, mg%n_cycle_base)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

       ! Have to fill ghost cells after correction
       call tree_2d_gc_ids(tree, tree%lvls(lvl)%ids, mg%i_phi, &
            mg%sides_rb, mg%sides_bc)

       ! Upwards relaxation
       call gsrb_boxes(tree, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do

    if (set_residual) then
       !$omp parallel private(lvl, i, id)
       do lvl = min_lvl, max_lvl
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call residual_box(tree%boxes(id), mg)
          end do
          !$omd end do nowait
       end do
       !$omp end parallel
    end if
  end subroutine mg2d_fas_vcycle

    !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg2d_box_gsrb_lpl(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2d_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: dx2
#if $D == 3
    integer                     :: k
    real(dp), parameter         :: sixth = 1/6.0_dp
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             box%cc(i, j, k, i_phi) = sixth * ( &
                  box%cc(i+1, j, k, i_phi) + box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + box%cc(i, j, k-1, i_phi) - &
                  dx2 * box%cc(i, j, k, i_rhs))
          end do
       end do
    end do
#endif
  end subroutine mg2d_box_gsrb_lpl

  !> Perform Laplacian operator on a box
  subroutine mg2d_box_lpl(box, i_out, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2d_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: inv_dr_sq
#if $D == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi = mg%i_phi

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
               box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
               box%cc(i, j+1, i_phi) - 4 * box%cc(i, j, i_phi))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, k, i_out) = inv_dr_sq * (box%cc(i-1, j, k, i_phi) + &
                  box%cc(i+1, j, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j, k-1, i_phi) + &
                  box%cc(i, j, k+1, i_phi) - 6 * box%cc(i, j, k, i_phi))
          end do
       end do
    end do
#endif
  end subroutine mg2d_box_lpl

    ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
  ! correction later
  subroutine update_coarse(tree, lvl, mg)
    use m_a$D_utils, only: a$D_n_cell, a$D_box_add_cc, a$D_box_copy_cc
    use m_a$D_ghostcell, only: a$D_gc_ids
    type(a$D_t), intent(inout) :: tree !< Tree containing full grid
    integer, intent(in)        :: lvl !< Update coarse values at lvl-1
    type(mg$D_t), intent(in)   :: mg !< Multigrid options
    integer                    :: i, id, p_id, nc
#if $D == 2
    real(dp), allocatable :: tmp(:,:)
#elif $D == 3
    real(dp), allocatable :: tmp(:,:,:)
#endif

    id = tree%lvls(lvl)%ids(1)
    nc = a$D_n_cell(tree, lvl)
#if $D == 2
    allocate(tmp(1:nc, 1:nc))
#elif $D == 3
    allocate(tmp(1:nc, 1:nc, 1:nc))
#endif

    ! Restrict phi and the residual
    !$omp parallel do private(id, p_id, tmp)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       p_id = tree%boxes(id)%parent

       ! Copy the data currently in i_tmp, and restore it later (i_tmp holds the
       ! previous state of i_phi)
#if $D == 2
       tmp = tree%boxes(id)%cc(1:nc, 1:nc, mg%i_tmp)
#elif $D == 3
       tmp = tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg%i_tmp)
#endif
       call residual_box(tree%boxes(id), mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp, mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi, mg)
#if $D == 2
       tree%boxes(id)%cc(1:nc, 1:nc, mg%i_tmp) = tmp
#elif $D == 3
       tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg%i_tmp) = tmp
#endif
    end do
    !$omp end parallel do

    call a$D_gc_ids(tree, tree%lvls(lvl-1)%ids, mg%i_phi, &
         mg%sides_rb, mg%sides_bc)

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
    ! store current coarse phi in tmp.

    !$omp parallel do private(id)
    do i = 1, size(tree%lvls(lvl-1)%parents)
       id = tree%lvls(lvl-1)%parents(i)

       ! Set rhs = L phi
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)

       ! Add tmp (the fine grid residual) to rhs
       call a$D_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)

       ! Story a copy of phi in tmp
       call a$D_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_tmp)
    end do
    !$omp end parallel do
  end subroutine update_coarse

    subroutine gsrb_boxes(tree, ids, mg, n_cycle)
    use m_a$D_ghostcell, only: a$D_gc_box
    type(a$D_t), intent(inout) :: tree    !< Tree containing full grid
    type(mg$D_t), intent(in)   :: mg      !< Multigrid options
    integer, intent(in)        :: ids(:)  !< Operate on these boxes
    integer, intent(in)        :: n_cycle !< Number of cycles to perform
    integer                    :: n, i

    do n = 1, 2 * n_cycle
       do i = 1, size(ids)
          call mg%box_gsrb(tree%boxes(ids(i)), n, mg)
       end do

       do i = 1, size(ids)
          call a$D_gc_box(tree, ids(i), mg%i_phi, mg%sides_rb, &
               mg%sides_bc)
       end do
    end do
  end subroutine gsrb_boxes

end module m_octree_mg
