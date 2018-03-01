module m_prolong
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: prolong
  public :: prolong_buffer_size

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine prolong_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(out)         :: n_send(0:mg%n_cpu-1)
    integer, intent(out)         :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)         :: dsize
    integer                      :: lvl

    if (.not. allocated(mg%comm_restrict%n_send)) then
       error stop "Call restrict_buffer_size before prolong_buffer_size"
    end if

    allocate(mg%comm_prolong%n_send(0:mg%n_cpu-1, mg%highest_lvl))
    allocate(mg%comm_prolong%n_recv(0:mg%n_cpu-1, mg%highest_lvl))

    mg%comm_prolong%n_recv(:, mg%highest_lvl) = 0
    mg%comm_prolong%n_send(:, mg%highest_lvl) = 0

    do lvl = 1, mg%highest_lvl-1
       mg%comm_prolong%n_recv(:, lvl) = &
            mg%comm_restrict%n_send(:, lvl+1)
       mg%comm_prolong%n_send(:, lvl) = &
            mg%comm_restrict%n_recv(:, lvl+1)
    end do

    ! Send coarse grid points + ghost cells
    dsize = (mg%box_size/2 + 2)**2
    n_send = maxval(mg%comm_prolong%n_send, dim=2)
    n_recv = maxval(mg%comm_prolong%n_recv, dim=2)

  end subroutine prolong_buffer_size

  subroutine prolong(mg, lvl)
    use m_communication
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: lvl
    integer                      :: i, id, dsize

    if (lvl == mg%highest_lvl) error stop "cannot prolong highest level"
    if (lvl < 1) error stop "cannot prolong level < 1"

    dsize            = (mg%box_size/2 + 2)**2
    mg%buf(:)%i_send = 0
    mg%buf(:)%i_ix   = 0

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call prolong_set_buffer(mg, id)
    end do

    mg%buf(:)%i_recv = mg%comm_prolong%n_recv(:, lvl) * dsize
    call sort_and_transfer_buffers(mg, dsize)
    mg%buf(:)%i_recv = 0

    do i = 1, size(mg%lvls(lvl+1)%my_ids)
       id = mg%lvls(lvl+1)%my_ids(i)
       call prolong_onto(mg, id)
    end do
  end subroutine prolong

  subroutine prolong_set_buffer(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: i, hnc, dix(2)
    integer                      :: i_c, c_id, c_rank, dsize

    hnc   = mg%box_size/2
    dsize = (hnc+2)**2

    do i_c = 1, num_children
       c_id = mg%boxes(id)%children(i_c)
       if (c_id > no_box) then
          c_rank = mg%boxes(c_id)%rank
          if (c_rank /= mg%my_rank) then
             dix = get_child_offset(mg, c_id)

             i = mg%buf(c_rank)%i_send
             mg%buf(c_rank)%send(i+1:i+dsize) = &
                  pack(mg%boxes(id)%cc(dix(1):dix(1)+hnc+1, &
                  dix(2):dix(2)+hnc+1, i_phi), .true.)
             mg%buf(c_rank)%i_send = mg%buf(c_rank)%i_send + dsize

             i = mg%buf(c_rank)%i_ix
             mg%buf(c_rank)%ix(i+1) = c_id
             mg%buf(c_rank)%i_ix = mg%buf(c_rank)%i_ix + 1
          end if
       end if
    end do
  end subroutine prolong_set_buffer

  subroutine prolong_onto(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: i, j, hnc, p_id, p_rank, dix(2), dsize
    real(dp)                     :: tmp(0:mg%box_size/2+1, 0:mg%box_size/2+1)
    real(dp)                     :: f0, flx, fhx, fly, fhy

    hnc    = mg%box_size/2
    dsize  = (hnc+2)**2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank
    dix    = get_child_offset(mg, id)

    if (p_rank == mg%my_rank) then
       tmp = mg%boxes(p_id)%cc(dix(1):dix(1)+hnc+1, &
                  dix(2):dix(2)+hnc+1, i_phi)
    else
       i = mg%buf(p_rank)%i_recv
       tmp = reshape(mg%buf(p_rank)%recv(i+1:i+dsize), [hnc+2, hnc+2])
       mg%buf(p_rank)%i_recv = mg%buf(p_rank)%i_recv + dsize
    end if

    do j = 1, hnc
       do i = 1, hnc
          f0 = 0.5_dp * tmp(i, j)
          flx = 0.25_dp * tmp(i-1, j)
          fhx = 0.25_dp * tmp(i+1, j)
          fly = 0.25_dp * tmp(i, j-1)
          fhy = 0.25_dp * tmp(i, j+1)

          mg%boxes(id)%cc(2*i-1, 2*j-1, i_phi) = f0 + flx + fly
          mg%boxes(id)%cc(2*i  , 2*j-1, i_phi) = f0 + fhx + fly
          mg%boxes(id)%cc(2*i-1, 2*j,   i_phi) = f0 + flx + fhy
          mg%boxes(id)%cc(2*i  , 2*j,   i_phi) = f0 + fhx + fhy
       end do
    end do

  end subroutine prolong_onto

end module m_prolong
