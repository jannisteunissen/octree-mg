module m_restrict
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: restrict
  public :: restrict_buffer_size

contains

  !> Specify minimum buffer size (per process) for communication
  subroutine restrict_buffer_size(mg, n_send, n_recv, dsize)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(out)         :: n_send(0:mg%n_cpu-1)
    integer, intent(out)         :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)         :: dsize
    integer                      :: n_out(0:mg%n_cpu-1, mg%highest_lvl)
    integer                      :: n_in(0:mg%n_cpu-1, mg%highest_lvl)
    integer                      :: lvl, i, id, p_id, p_rank
    integer                      :: i_c, c_id, c_rank

    n_out(:, :) = 0
    n_in(:, :)  = 0

    do lvl = 2, mg%highest_lvl
       ! Number of messages to receive (at lvl-1)
       do i = 1, size(mg%lvls(lvl-1)%my_parents)
          id = mg%lvls(lvl-1)%my_parents(i)
          do i_c = 1, num_children
             c_id = mg%boxes(id)%children(i_c)

             if (c_id > no_box) then
                c_rank = mg%boxes(c_id)%rank
                if (c_rank /= mg%my_rank) then
                   n_in(c_rank, lvl) = n_in(c_rank, lvl) + 1
                end if
             end if
          end do
       end do

       ! Number of messages to send
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)

          p_id = mg%boxes(id)%parent
          p_rank = mg%boxes(p_id)%rank
          if (p_rank /= mg%my_rank) then
             n_out(p_rank, lvl) = n_out(p_rank, lvl) + 1
          end if

       end do
    end do

    allocate(mg%comm_restrict%n_send(0:mg%n_cpu-1, mg%highest_lvl))
    allocate(mg%comm_restrict%n_recv(0:mg%n_cpu-1, mg%highest_lvl))
    mg%comm_restrict%n_send = n_out
    mg%comm_restrict%n_recv = n_in

    dsize  = (mg%box_size/2)**2
    n_send = maxval(n_out, dim=2)
    n_recv = maxval(n_in, dim=2)
  end subroutine restrict_buffer_size

  subroutine restrict(mg, lvl)
    use m_communication
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: lvl
    integer                      :: i, id, dsize

    if (lvl <= 1) error stop "cannot restrict for lvl <= 1"

    dsize = (mg%box_size/2)**2

    mg%buf(:)%i_send = 0
    mg%buf(:)%i_ix   = 0

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call restrict_set_buffer(mg, id)
    end do

    mg%buf(:)%i_recv = mg%comm_restrict%n_recv(:, lvl) * dsize
    call sort_and_transfer_buffers(mg, dsize)
    mg%buf(:)%i_recv = 0

    do i = 1, size(mg%lvls(lvl-1)%my_parents)
       id = mg%lvls(lvl-1)%my_parents(i)
       call restrict_onto(mg, id)
    end do
  end subroutine restrict

  subroutine restrict_set_buffer(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    real(dp)                     :: tmp(mg%box_size/2, mg%box_size/2)
    integer                      :: i, j, n, hnc, p_id, p_rank

    hnc    = mg%box_size/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank /= mg%my_rank) then
       do j = 1, hnc
          do i = 1, hnc
             tmp(i, j) = 0.25_dp * &
                  sum(mg%boxes(id)%cc(2*i-1:2*i, 2*j-1:2*j, i_phi))
          end do
       end do

       ! Buffer
       n = hnc**2
       i = mg%buf(p_rank)%i_send
       mg%buf(p_rank)%send(i+1:i+n) = pack(tmp, .true.)
       mg%buf(p_rank)%i_send = mg%buf(p_rank)%i_send + n

       ! To later sort the send buffer according to parent order
       i = mg%buf(p_rank)%i_ix
       n = ix_to_ichild(mg%boxes(id)%ix)
       mg%buf(p_rank)%ix(i+1) = num_children * p_id + n
       mg%buf(p_rank)%i_ix = mg%buf(p_rank)%i_ix + 1
    end if
  end subroutine restrict_set_buffer

  subroutine restrict_onto(mg, id)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer                      :: i, j, hnc, dsize, i_c, c_id
    integer                      :: c_rank, dix(2)

    hnc   = mg%box_size/2
    dsize = hnc**2

    do i_c = 1, num_children
       c_id   = mg%boxes(id)%children(i_c)
       c_rank = mg%boxes(c_id)%rank
       dix    = get_child_offset(mg, c_id)

       if (c_rank == mg%my_rank) then
          do j = 1, hnc
             do i = 1, hnc
                mg%boxes(id)%cc(dix(1)+i, dix(2)+j, i_phi) = 0.25_dp * &
                     sum(mg%boxes(c_id)%cc(2*i-1:2*i, 2*j-1:2*j, i_phi))
             end do
          end do
       else
          i = mg%buf(c_rank)%i_recv
          mg%boxes(id)%cc(dix(1)+1:dix(1)+hnc, &
               dix(2)+1:dix(2)+hnc, i_phi) = &
               reshape(mg%buf(c_rank)%recv(i+1:i+dsize), [hnc, hnc])
          mg%buf(c_rank)%i_recv = mg%buf(c_rank)%i_recv + dsize
       end if
    end do

  end subroutine restrict_onto

end module m_restrict
