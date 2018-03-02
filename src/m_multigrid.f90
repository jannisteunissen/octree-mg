module m_multigrid
  use m_data_structures
  use m_prolong
  use m_restrict
  use m_ghost_cells

  implicit none
  private

  ! Public methods
  public :: mg_fas_vcycle

contains

  subroutine mg_fas_vcycle(mg, set_residual, highest_lvl)
    type(mg_2d_t), intent(inout)  :: mg
    logical, intent(in)           :: set_residual !< If true, store residual in i_res
    integer, intent(in), optional :: highest_lvl  !< Maximum level for V-cycle
    integer                       :: lvl, min_lvl, i, id, max_lvl

    min_lvl = 1
    max_lvl = mg%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(mg, lvl, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_old for the
       ! correction later
       call update_coarse(mg, lvl)
    end do

    lvl = min_lvl
    call gsrb_boxes(mg, lvl, mg%n_cycle_base)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(mg, lvl-1)

       ! Have to fill ghost cells after correction
       call fill_ghost_cells_lvl(mg, lvl)

       ! Upwards relaxation
       call gsrb_boxes(mg, lvl, mg%n_cycle_up)
    end do

    if (set_residual) then
       do lvl = min_lvl, max_lvl
          do i = 1, size(mg%lvls(lvl)%my_ids)
             id = mg%lvls(lvl)%my_ids(i)
             call residual_box(mg, id)
          end do
       end do
    end if
  end subroutine mg_fas_vcycle

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine box_gsrb_lpl(mg, id, redblack_cntr)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: redblack_cntr !< Iteration counter
    integer                      :: i, i0, j, nc
    real(dp)                     :: dx2

    dx2   = mg%dr(mg%boxes(id)%lvl)**2
    nc    = mg%box_size

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (box => mg%boxes(id))
      do j = 1, nc
         i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, 2
            box%cc(i, j, i_phi) = 0.25_dp * ( &
                 box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
                 box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
                 dx2 * box%cc(i, j, i_rhs))
         end do
      end do
    end associate
  end subroutine box_gsrb_lpl

  !> Perform Laplacian operator on a box
  subroutine box_lpl(mg, id, i_out)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: i_out !< Index of variable to store Laplacian in
    integer                      :: i, j, nc
    real(dp)                     :: inv_dr_sq

    nc    = mg%box_size
    inv_dr_sq = 1 / mg%dr(mg%boxes(id)%lvl)**2

    associate (box => mg%boxes(id))
      do j = 1, nc
         do i = 1, nc
            box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
                 box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
                 box%cc(i, j+1, i_phi) - 4 * box%cc(i, j, i_phi))
         end do
      end do
    end associate
  end subroutine box_lpl

  ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_old for the
  ! correction later
  subroutine update_coarse(mg, lvl)
    type(mg_2d_t), intent(inout) :: mg  !< Tree containing full grid
    integer, intent(in)          :: lvl !< Update coarse values at lvl-1
    integer                      :: i, id, nc

    nc = mg%box_size

    ! Compute residual
    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call residual_box(mg, id)
    end do

    ! Restrict phi and the residual
    call restrict(mg, i_phi, lvl)
    call restrict(mg, i_res, lvl)

    call fill_ghost_cells_lvl(mg, lvl-1)

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
    ! store current coarse phi in old.
    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)

       ! Set rhs = L phi
       call box_lpl(mg, id, i_rhs)

       ! Add tmp (the fine grid residual) to rhs
       mg%boxes(id)%cc(1:nc, 1:nc, i_rhs) = &
            mg%boxes(id)%cc(1:nc, 1:nc, i_rhs) + &
            mg%boxes(id)%cc(1:nc, 1:nc, i_res)

       ! Story a copy of phi in tmp
       mg%boxes(id)%cc(:, :, i_old) = &
            mg%boxes(id)%cc(:, :, i_phi)
    end do
  end subroutine update_coarse

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(mg, lvl)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: lvl
    integer                      :: i, id

    do i = 1, size(mg%lvls(lvl)%my_parents)
       id = mg%lvls(lvl)%my_parents(i)

       ! Store the correction in i_res
       mg%boxes(id)%cc(:, :, i_res) = &
            mg%boxes(id)%cc(:, :, i_phi) - &
            mg%boxes(id)%cc(:, :, i_old)
    end do

    call prolong(mg, lvl, i_res, i_phi, add=.true.)
  end subroutine correct_children

  subroutine gsrb_boxes(mg, lvl, n_cycle)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: lvl
    integer, intent(in)          :: n_cycle !< Number of cycles to perform
    integer                      :: n, i, id

    do n = 1, 2 * n_cycle
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call box_gsrb_lpl(mg, id, n)
       end do

       call fill_ghost_cells_lvl(mg, lvl)
    end do
  end subroutine gsrb_boxes

  subroutine residual_box(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: nc

    nc = mg%box_size

    call box_lpl(mg, id, i_res)

    mg%boxes(id)%cc(1:nc, 1:nc, i_res) = &
         mg%boxes(id)%cc(1:nc, 1:nc, i_rhs) &
         - mg%boxes(id)%cc(1:nc, 1:nc, i_res)
  end subroutine residual_box

end module m_multigrid
