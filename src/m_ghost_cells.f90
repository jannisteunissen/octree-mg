module m_ghost_cells
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: fill_ghost_cells_lvl

contains

  subroutine fill_ghost_cells_lvl(mg, lvl)
    use mpi
    type(mg_2d_t)        :: mg
    integer, intent(in)  :: lvl
    integer              :: i, id, n
    integer :: send_req(0:mg%n_cpu-1)
    integer :: recv_req(0:mg%n_cpu-1)
    integer              :: ierr

    if (lvl < 1) error stop "fill_ghost_cells_lvl: lvl < 1"
    if (lvl > mg%highest_lvl) error stop "fill_ghost_cells_lvl: lvl > highest_lvl"

    mg%gc(:)%i_send = 0
    mg%gc(:)%i_recv = 0
    mg%gc(:)%i_id = 0

    n = size(mg%lvls(lvl)%my_ids) * 4

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call set_and_buffer_ghost_cells(mg, id)
    end do

    ! Transfer data between processes

    do i = 0, mg%n_cpu - 1
       call sort_sendbuf(mg, mg%gc(i))
       ! print *, mg%my_rank, "receiving from", i, mg%gc(i)%idbuf(1:mg%gc(i)%i_id)
       call mpi_isend(mg%gc(i)%sendbuf, mg%gc(i)%i_send, MPI_DOUBLE, &
            i, 0, MPI_COMM_WORLD, send_req(i), ierr)
       call mpi_irecv(mg%gc(i)%recvbuf, mg%gc(i)%i_send, MPI_DOUBLE, &
            i, 0, MPI_COMM_WORLD, recv_req(i), ierr)
    end do

    call mpi_waitall(mg%n_cpu, recv_req, MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(mg%n_cpu, send_req, MPI_STATUSES_IGNORE, ierr)

    ! Set ghost cells to received data
    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call fill_buffered_ghost_cells(mg, id)
    end do

  end subroutine fill_ghost_cells_lvl

  subroutine sort_sendbuf(mg, gc)
    use m_mrgrnk
    type(mg_2d_t), intent(inout)  :: mg
    type(gc_buf_t), intent(inout) :: gc
    integer                       :: ix_sort(gc%i_id)
    real(dp)                      :: buf_cpy(gc%i_send)
    integer :: i, j, n

    call mrgrnk(gc%idbuf(1:gc%i_id), ix_sort)

    buf_cpy = gc%sendbuf(1:gc%i_send)

    do n = 1, gc%i_id
       i = (n-1) * mg%box_size
       j = (ix_sort(n)-1) * mg%box_size
       gc%sendbuf(i+1:i+mg%box_size) = buf_cpy(j+1:j+mg%box_size)
    end do
    gc%idbuf(1:gc%i_id) = gc%idbuf(ix_sort)

  end subroutine sort_sendbuf

  subroutine set_and_buffer_ghost_cells(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: nb, nb_id, nb_rank

    do nb = 1, num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank == mg%my_rank) then
             call copy_from_nb(mg, mg%boxes(id), mg%boxes(nb_id), nb)
          else
             ! call sendrecv_with_nb(mg, mg%boxes(id), nb_rank, nb)
             call buffer_for_nb(mg, mg%boxes(id), nb_id, nb_rank, nb)
          end if

          ! else if (nb_id == af_no_box) then
          !    ! Refinement boundary
          !    call subr_rb(boxes, id, nb, iv)
          ! else
          !    ! Physical boundary
          !    call subr_bc(boxes(id), nb, iv, bc_type)
          !    call bc_to_gc(boxes(id), nb, iv, bc_type)
       end if
    end do
  end subroutine set_and_buffer_ghost_cells

  subroutine fill_buffered_ghost_cells(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: nb, nb_id, nb_rank

    do nb = 1, num_neighbors
       nb_id = mg%boxes(id)%neighbors(nb)

       if (nb_id > no_box) then
          ! There is a neighbor
          nb_rank    = mg%boxes(nb_id)%rank

          if (nb_rank /= mg%my_rank) then
             call fill_buffered_nb(mg, mg%boxes(id), nb_rank, nb)
          end if
       end if
    end do
  end subroutine fill_buffered_ghost_cells

  subroutine copy_from_nb(mg, box, box_nb, nb)
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    type(box_2d_t), intent(in)    :: box_nb
    integer, intent(in)           :: nb
    real(dp)                      :: gc(mg%box_size)

    call box_gc_from_neighbor(box_nb, nb, mg%box_size, gc)
    call box_set_gc(box, nb, mg%box_size, gc)
  end subroutine copy_from_nb

  subroutine buffer_for_nb(mg, box, nb_id, nb_rank, nb)
    use mpi
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    integer, intent(in)           :: nb_id
    integer, intent(in)           :: nb_rank
    integer, intent(in)           :: nb
    integer                       :: i

    i = mg%gc(nb_rank)%i_send
    call box_gc_for_neighbor(box, nb, mg%box_size, &
         mg%gc(nb_rank)%sendbuf(i+1:i+mg%box_size))

    ! Later the buffer is sorted, using the fact that loops go from low to high
    ! box id, and we fill ghost cells according to the neighbor order
    i = mg%gc(nb_rank)%i_id
    mg%gc(nb_rank)%idbuf(i+1) = num_neighbors * nb_id + neighb_rev(nb)

    mg%gc(nb_rank)%i_send = mg%gc(nb_rank)%i_send + mg%box_size
    mg%gc(nb_rank)%i_id   = mg%gc(nb_rank)%i_id + 1
  end subroutine buffer_for_nb

  subroutine fill_buffered_nb(mg, box, nb_rank, nb)
    use mpi
    type(mg_2d_t), intent(inout)  :: mg
    type(box_2d_t), intent(inout) :: box
    integer, intent(in)           :: nb_rank
    integer, intent(in)           :: nb
    integer                       :: i

    i = mg%gc(nb_rank)%i_recv
    call box_set_gc(box, nb, mg%box_size, &
         mg%gc(nb_rank)%recvbuf(i+1:i+mg%box_size))
    mg%gc(nb_rank)%i_recv = mg%gc(nb_rank)%i_recv + mg%box_size

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
       gc =box%cc(1:nc, nc, i_phi)
    end select
  end subroutine box_gc_for_neighbor

  subroutine box_gc_from_neighbor(box_nb, nb, nc, gc)
    type(box_2d_t), intent(in) :: box_nb
    integer, intent(in)        :: nb, nc
    real(dp), intent(out)      :: gc(nc)

    select case (nb)
    case (neighb_lowx)
       gc = box_nb%cc(nc, 1:nc, i_phi)
    case (neighb_highx)
       gc = box_nb%cc(1, 1:nc, i_phi)
    case (neighb_lowy)
       gc = box_nb%cc(1:nc, nc, i_phi)
    case (neighb_highy)
       gc = box_nb%cc(1:nc, 1, i_phi)
    end select
  end subroutine box_gc_from_neighbor

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

  ! subroutine sendrecv_with_nb(mg, box, nb_rank, nb, recv_req, send_req)
  !   use mpi
  !   type(mg_2d_t), intent(inout)  :: mg
  !   type(box_2d_t), intent(inout) :: box
  !   integer, intent(in)           :: nb_rank
  !   integer, intent(in)           :: nb
  !   integer, intent(inout)        :: recv_req, send_req
  !   integer                       :: ierr
  !   real(dp)                      :: buf(mg%box_size)

  !   call mpi_irecv(box%gc(:, nb), mg%box_size, MPI_DOUBLE, nb_rank, &
  !        0, MPI_COMM_WORLD, recv_req, ierr)

  !   associate (l => mg%ix_sc_lo(:, nb), h => mg%ix_sc_hi(:, nb))
  !     buf = pack(box%cc(l(1):h(1), l(2):h(2), i_phi), .true.)
  !   end associate

  !   call mpi_ibsend(buf, mg%box_size, MPI_DOUBLE, &
  !           nb_rank, 0, MPI_COMM_WORLD, send_req, ierr)
  ! end subroutine sendrecv_with_nb

  ! subroutine create_mpi_types(mg)
  !   type(mg_2d_t), intent(in) :: mg
  !   integer :: block_size(2)

  !   block_size = mg%block_size

  !   call mpi_type_vector(1, nc, 
  ! end subroutine create_mpi_types

end module m_ghost_cells
