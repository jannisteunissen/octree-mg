module m_communication
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: comm_init
  public :: sort_and_transfer_buffers

contains

  subroutine comm_init(mg)
    use mpi
    type(mg_t), intent(inout) :: mg
    integer                      :: ierr

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, mg%my_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mg%n_cpu, ierr)
  end subroutine comm_init

  subroutine sort_and_transfer_buffers(mg, dsize)
    use mpi
    type(mg_t), intent(inout) :: mg
    integer, intent(in)          :: dsize
    integer                      :: i
    integer                      :: send_req(0:mg%n_cpu-1)
    integer                      :: recv_req(0:mg%n_cpu-1)
    integer                      :: ierr

    do i = 0, mg%n_cpu - 1
       call sort_sendbuf(mg%buf(i), dsize)
       call mpi_isend(mg%buf(i)%send, mg%buf(i)%i_send, MPI_DOUBLE, &
            i, 0, MPI_COMM_WORLD, send_req(i), ierr)
       call mpi_irecv(mg%buf(i)%recv, mg%buf(i)%i_recv, MPI_DOUBLE, &
            i, 0, MPI_COMM_WORLD, recv_req(i), ierr)
    end do

    call mpi_waitall(mg%n_cpu, recv_req, MPI_STATUSES_IGNORE, ierr)
    call mpi_waitall(mg%n_cpu, send_req, MPI_STATUSES_IGNORE, ierr)

  end subroutine sort_and_transfer_buffers

  !> Sort send buffers according to the idbuf array
  subroutine sort_sendbuf(gc, dsize)
    use m_mrgrnk
    type(buf_t), intent(inout) :: gc
    integer, intent(in)        :: dsize !< Size of send buffer elements
    integer                    :: ix_sort(gc%i_ix)
    real(dp)                   :: buf_cpy(gc%i_send)
    integer                    :: i, j, n

    call mrgrnk(gc%ix(1:gc%i_ix), ix_sort)

    buf_cpy = gc%send(1:gc%i_send)

    do n = 1, gc%i_ix
       i = (n-1) * dsize
       j = (ix_sort(n)-1) * dsize
       gc%send(i+1:i+dsize) = buf_cpy(j+1:j+dsize)
    end do
    gc%ix(1:gc%i_ix) = gc%ix(ix_sort)

  end subroutine sort_sendbuf

end module m_communication
