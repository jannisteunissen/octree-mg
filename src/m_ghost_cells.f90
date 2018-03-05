module m_ghost_cells
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: ghost_cell_buffer_size
  public :: fill_ghost_cells_lvl
  public :: set_bc_continuous
  public :: set_bc_dirichlet_zero
  public :: set_bc_neumann_zero

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine ghost_cell_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(out)         :: n_send(0:mg%n_cpu-1)
    integer, intent(out)         :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)         :: dsize
    integer                      :: i, id, lvl

    allocate(mg%comm_ghostcell%n_send(0:mg%n_cpu-1, mg%highest_lvl))
    allocate(mg%comm_ghostcell%n_recv(0:mg%n_cpu-1, mg%highest_lvl))

    dsize = mg%box_size

    do lvl = 1, mg%highest_lvl
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_recv = 0
       mg%buf(:)%i_ix = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call buffer_ghost_cells(mg, id, dry_run=.true.)
       end do

       if (lvl > 1) then
          do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
             id = mg%lvls(lvl-1)%my_ref_bnds(i)
             call buffer_refinement_boundaries(mg, id, dry_run=.true.)
          end do
       end if

       ! Set ghost cells to received data
       mg%buf(:)%i_recv = 0
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call set_ghost_cells(mg, id, dry_run=.true.)
       end do

       mg%comm_ghostcell%n_send(:, lvl) = mg%buf(:)%i_send/dsize
       mg%comm_ghostcell%n_recv(:, lvl) = mg%buf(:)%i_recv/dsize
       ! print *, mg%my_rank, "send", lvl, mg%comm_ghostcell%n_send(:, lvl)
       ! print *, mg%my_rank, "recv", lvl, mg%comm_ghostcell%n_recv(:, lvl)
    end do

    n_send = maxval(mg%comm_ghostcell%n_send, dim=2)
    n_recv = maxval(mg%comm_ghostcell%n_recv, dim=2)
  end subroutine ghost_cell_buffer_size

  subroutine fill_ghost_cells_lvl(mg, lvl)
    use m_communication
    type(mg_2d_t)        :: mg
    integer, intent(in)  :: lvl
    integer              :: i, id

    if (lvl < 1) error stop "fill_ghost_cells_lvl: lvl < 1"
    if (lvl > mg%highest_lvl) error stop "fill_ghost_cells_lvl: lvl > highest_lvl"

    mg%buf(:)%i_send = 0
    mg%buf(:)%i_recv = 0
    mg%buf(:)%i_ix = 0

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call buffer_ghost_cells(mg, id, .false.)
    end do

    if (lvl > 1) then
       do i = 1, size(mg%lvls(lvl-1)%my_ref_bnds)
          id = mg%lvls(lvl-1)%my_ref_bnds(i)
          call buffer_refinement_boundaries(mg, id, .false.)
       end do
    end if

    ! Transfer data between processes
    mg%buf(:)%i_recv = mg%comm_ghostcell%n_recv(:, lvl) * mg%box_size
    call sort_and_transfer_buffers(mg, mg%box_size)

    ! Set ghost cells to received data
    mg%buf(:)%i_recv = 0
    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call set_ghost_cells(mg, id, .false.)
    end do

  end subroutine fill_ghost_cells_lvl

  subroutine buffer_ghost_cells(mg, id, dry_run)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    logical, intent(in)          :: dry_run
    integer                      :: nb, nb_id, nb_rank

    do nb = 1, num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call buffer_for_nb(mg, mg%boxes(id), nb_id, nb_rank, nb, dry_run)
          end if
       end if
    end do
  end subroutine buffer_ghost_cells

  subroutine buffer_refinement_boundaries(mg, id, dry_run)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    logical, intent(in)          :: dry_run
    integer                      :: nb, nb_id, c_ids(2), n, c_id, c_rank

    do nb = 1, num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)
       if (nb_id > no_box) then
          if (has_children(mg%boxes(nb_id))) then
             c_ids = mg%boxes(nb_id)%children(child_adj_nb(:, neighb_rev(nb)))

             do n = 1, num_children/2
                c_id = c_ids(n)
                c_rank = mg%boxes(c_id)%rank

                if (c_rank /= mg%my_rank) then
                   ! Send all coarse ghost cells
                   call buffer_for_nb(mg, mg%boxes(id), c_id, c_rank, nb, dry_run)
                end if
             end do
          end if
       end if
    end do
  end subroutine buffer_refinement_boundaries

  subroutine set_ghost_cells(mg, id, dry_run)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    logical, intent(in)          :: dry_run
    integer                      :: nb, nb_id, nb_rank, bc_type

    do nb = 1, num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call fill_buffered_nb(mg, mg%boxes(id), nb_rank, nb, dry_run)
          else if (.not. dry_run) then
             call copy_from_nb(mg, mg%boxes(id), mg%boxes(nb_id), nb)
          end if
       else if (nb_id == no_box) then
          ! Refinement boundary
          call fill_refinement_bnd(mg, id, nb, dry_run)
       else if (.not. dry_run) then
          ! Physical boundary
          if (.not. associated(mg%boundary_cond)) then
             error stop "Physical boundary but mg%boundary_cond not set"
          end if
          call mg%boundary_cond(mg, id, nb, bc_type)
          call bc_to_gc(mg, id, nb, bc_type)
       end if
    end do
  end subroutine set_ghost_cells

  subroutine fill_refinement_bnd(mg, id, nb, dry_run)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb
    logical, intent(in)          :: dry_run
    real(dp)                     :: cgc(mg%box_size)
    integer                      :: p_id, p_nb_id
    integer                      :: i, dsize, p_nb_rank

    dsize     = mg%box_size
    p_id      = mg%boxes(id)%parent
    p_nb_id   = mg%boxes(p_id)%neighbors(nb)
    p_nb_rank = mg%boxes(p_nb_id)%rank

    if (p_nb_rank /= mg%my_rank) then
       i = mg%buf(p_nb_rank)%i_recv
       if (.not. dry_run) then
          cgc = mg%buf(p_nb_rank)%recv(i+1:i+dsize)
       end if
       mg%buf(p_nb_rank)%i_recv = mg%buf(p_nb_rank)%i_recv + dsize
    else if (.not. dry_run) then
       call box_gc_for_neighbor(mg%boxes(p_nb_id), neighb_rev(nb), &
            mg%box_size, cgc)
    end if

    if (.not. dry_run) then
       call sides_rb(mg, id, nb, cgc)
    end if
  end subroutine fill_refinement_bnd

  subroutine copy_from_nb(mg, box, box_nb, nb)
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    type(box_2d_t), intent(in)    :: box_nb
    integer, intent(in)           :: nb
    real(dp)                      :: gc(mg%box_size)

    call box_gc_for_neighbor(box_nb, neighb_rev(nb), mg%box_size, gc)
    call box_set_gc(box, nb, mg%box_size, gc)
  end subroutine copy_from_nb

  subroutine buffer_for_nb(mg, box, nb_id, nb_rank, nb, dry_run)
    use mpi
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    integer, intent(in)           :: nb_id
    integer, intent(in)           :: nb_rank
    integer, intent(in)           :: nb
    logical, intent(in)           :: dry_run
    integer                       :: i

    i = mg%buf(nb_rank)%i_send
    if (.not. dry_run) then
       call box_gc_for_neighbor(box, nb, mg%box_size, &
            mg%buf(nb_rank)%send(i+1:i+mg%box_size))
    end if

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%buf(nb_rank)%i_ix
    if (.not. dry_run) then
       mg%buf(nb_rank)%ix(i+1) = num_neighbors * nb_id + neighb_rev(nb)
    end if

    mg%buf(nb_rank)%i_send = mg%buf(nb_rank)%i_send + mg%box_size
    mg%buf(nb_rank)%i_ix   = mg%buf(nb_rank)%i_ix + 1
  end subroutine buffer_for_nb

  subroutine fill_buffered_nb(mg, box, nb_rank, nb, dry_run)
    use mpi
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    integer, intent(in)           :: nb_rank
    integer, intent(in)           :: nb
    logical, intent(in)           :: dry_run
    integer                       :: i

    i = mg%buf(nb_rank)%i_recv
    if (.not. dry_run) then
       call box_set_gc(box, nb, mg%box_size, &
            mg%buf(nb_rank)%recv(i+1:i+mg%box_size))
    end if
    mg%buf(nb_rank)%i_recv = mg%buf(nb_rank)%i_recv + mg%box_size

  end subroutine fill_buffered_nb

  subroutine box_gc_for_neighbor(box, nb, nc, gc)
    type(box_2d_t), intent(in) :: box
    integer, intent(in)        :: nb, nc
    real(dp), intent(out)      :: gc(nc)

    select case (nb)
    case (neighb_lowx)
       gc = box%cc(1, 1:nc, i_phi)
    case (neighb_highx)
       gc = box%cc(nc, 1:nc, i_phi)
    case (neighb_lowy)
       gc = box%cc(1:nc, 1, i_phi)
    case (neighb_highy)
       gc = box%cc(1:nc, nc, i_phi)
    end select
  end subroutine box_gc_for_neighbor

  subroutine box_set_gc(box, nb, nc, gc)
    type(box_2d_t), intent(inout) :: box
    integer, intent(in)           :: nb, nc
    real(dp), intent(in)          :: gc(nc)

    select case (nb)
    case (neighb_lowx)
       box%cc(0, 1:nc, i_phi) = gc
    case (neighb_highx)
       box%cc(nc+1, 1:nc, i_phi) = gc
    case (neighb_lowy)
       box%cc(1:nc, 0, i_phi) = gc
    case (neighb_highy)
       box%cc(1:nc, nc+1, i_phi) = gc
    end select
  end subroutine box_set_gc

  subroutine bc_to_gc(mg, id, nb, bc_type)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb      !< Neighbor direction
    integer, intent(in)          :: bc_type !< Type of b.c.
    real(dp)                     :: c0, c1, c2
    integer                      :: nc

    nc = mg%box_size

    ! If we call the interior point x1, x2 and the ghost point x0, then a
    ! Dirichlet boundary value b can be imposed as:
    ! x0 = -x1 + 2*b
    ! A Neumann b.c. can be imposed as:
    ! x0 = x1 +/- dx * b
    ! A continuous boundary (same slope) as:
    ! x0 = 2 * x1 - x2
    ! Below, we set coefficients to handle these cases
    select case (bc_type)
    case (bc_dirichlet)
       c0 = 2
       c1 = -1
       c2 = 0
    case (bc_neumann)
       c0 = mg%dr(mg%boxes(id)%lvl) * neighb_high_pm(nb) ! This gives a + or - sign
       c1 = 1
       c2 = 0
    case (bc_continuous)
       c0 = 0
       c1 = 2
       c2 = -1
    case default
       error stop "bc_to_gc: unknown boundary condition"
    end select

    select case (nb)
    case (neighb_lowx)
       mg%boxes(id)%cc(0, 1:nc, i_phi) = &
            c0 * mg%boxes(id)%cc(0, 1:nc, i_phi) + &
            c1 * mg%boxes(id)%cc(1, 1:nc, i_phi) + &
            c2 * mg%boxes(id)%cc(2, 1:nc, i_phi)
    case (neighb_highx)
       mg%boxes(id)%cc(nc+1, 1:nc, i_phi) = &
            c0 * mg%boxes(id)%cc(nc+1, 1:nc, i_phi) + &
            c1 * mg%boxes(id)%cc(nc, 1:nc, i_phi) + &
            c2 * mg%boxes(id)%cc(nc-1, 1:nc, i_phi)
    case (neighb_lowy)
       mg%boxes(id)%cc(1:nc, 0, i_phi) = &
            c0 * mg%boxes(id)%cc(1:nc, 0, i_phi) + &
            c1 * mg%boxes(id)%cc(1:nc, 1, i_phi) + &
            c2 * mg%boxes(id)%cc(1:nc, 2, i_phi)
    case (neighb_highy)
       mg%boxes(id)%cc(1:nc, nc+1, i_phi) = &
            c0 * mg%boxes(id)%cc(1:nc, nc+1, i_phi) + &
            c1 * mg%boxes(id)%cc(1:nc, nc, i_phi) + &
            c2 * mg%boxes(id)%cc(1:nc, nc-1, i_phi)
    end select
  end subroutine bc_to_gc

    ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine set_bc_neumann_zero(mg, id, nb, bc_type)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb
    integer, intent(out)         :: bc_type
    integer                      :: i

    bc_type = bc_neumann
    call box_set_gc(mg%boxes(id), nb, mg%box_size, &
         [(0.0_dp, i=1, mg%box_size)])
  end subroutine set_bc_neumann_zero

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine set_bc_dirichlet_zero(mg, id, nb, bc_type)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb
    integer, intent(out)         :: bc_type
    integer                      :: i

    bc_type = bc_dirichlet
    call box_set_gc(mg%boxes(id), nb, mg%box_size, &
         [(0.0_dp, i=1, mg%box_size)])
  end subroutine set_bc_dirichlet_zero

  ! This fills ghost cells near physical boundaries using the same slope
  subroutine set_bc_continuous(mg, id, nb, bc_type)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb
    integer, intent(out)         :: bc_type
    integer                      :: i

    bc_type = bc_continuous
    ! Set values to zero (to prevent problems with NaN)
    call box_set_gc(mg%boxes(id), nb, mg%box_size, &
         [(0.0_dp, i=1, mg%box_size)])
  end subroutine set_bc_continuous

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  !>
  !> Basically, we extrapolate from the fine cells to a corner point, and then
  !> take the average between this corner point and a coarse neighbor to fill
  !> ghost cells for the fine cells.
  subroutine sides_rb(mg, id, nb, cgc)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id !< Id of box
    integer, intent(in)          :: nb !< Ghost cell direction
    !> Unmodified coarse grid ghost cells (including data for neighbors)
    real(dp), intent(in)         :: cgc(mg%box_size)
    integer                      :: nc, ix, dix, i, di, j, dj
    integer                      :: ix_off(2)

    nc     = mg%box_size
    ix_off = get_child_offset(mg, id)

    if (neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    select case (neighb_dim(nb))
    case (1)
       i = ix
       di = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          ! Extrapolation using 3 points
          mg%boxes(id)%cc(i-di, j, i_phi) = 0.5_dp * cgc(ix_off(2)+(j+1)/2) + &
               mg%boxes(id)%cc(i, j, i_phi) - 0.25_dp * &
               (mg%boxes(id)%cc(i+di, j, i_phi) + mg%boxes(id)%cc(i, j+dj, i_phi))

          ! Extrapolation using 2 points
          ! mg%boxes(id)%cc(i-di, j, i_phi) = 0.5_dp * mg%boxes(id)%cc(i-di, j, i_phi) + &
          !      0.75_dp * mg%boxes(id)%cc(i, j, i_phi) - 0.25_dp * &
          !      mg%boxes(id)%cc(i+di, j+dj, i_phi)
       end do
    case (2)
       j = ix
       dj = dix
       do i = 1, nc
          di = -1 + 2 * iand(i, 1)
          ! Extrapolation using 3 points
          mg%boxes(id)%cc(i, j-dj, i_phi) = 0.5_dp * cgc(ix_off(1)+(i+1)/2) + &
               mg%boxes(id)%cc(i, j, i_phi) - 0.25_dp * &
               (mg%boxes(id)%cc(i, j+dj, i_phi) + mg%boxes(id)%cc(i+di, j, i_phi))

          ! Extrapolation using 2 points
          ! mg%boxes(id)%cc(i, j-dj, i_phi) = 0.5_dp * mg%boxes(id)%cc(i, j-dj, i_phi) + &
          !      0.75_dp * mg%boxes(id)%cc(i, j, i_phi) - 0.25_dp * &
          !      mg%boxes(id)%cc(i+di, j+dj, i_phi)
       end do
    end select

  end subroutine sides_rb

end module m_ghost_cells
