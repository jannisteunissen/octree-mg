#include "cpp_macros.h"
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
    type(mg_t), intent(inout) :: mg
    integer, intent(out)      :: n_send(0:mg%n_cpu-1)
    integer, intent(out)      :: n_recv(0:mg%n_cpu-1)
    integer, intent(out)      :: dsize
    integer                   :: lvl, min_lvl

    if (.not. allocated(mg%comm_restrict%n_send)) then
       error stop "Call restrict_buffer_size before prolong_buffer_size"
    end if

    min_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)
    allocate(mg%comm_prolong%n_send(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))
    allocate(mg%comm_prolong%n_recv(0:mg%n_cpu-1, &
         min_lvl:mg%highest_lvl))

    mg%comm_prolong%n_recv(:, mg%highest_lvl) = 0
    mg%comm_prolong%n_send(:, mg%highest_lvl) = 0

    do lvl = min_lvl, mg%highest_lvl-1
       mg%comm_prolong%n_recv(:, lvl) = &
            mg%comm_restrict%n_send(:, lvl+1)
       mg%comm_prolong%n_send(:, lvl) = &
            mg%comm_restrict%n_recv(:, lvl+1)
    end do

    ! Send coarse grid points + ghost cells
    dsize = (mg%box_size/2 + 2)**NDIM
    n_send = maxval(mg%comm_prolong%n_send, dim=2)
    n_recv = maxval(mg%comm_prolong%n_recv, dim=2)
  end subroutine prolong_buffer_size

  !> Prolong variable iv from lvl to variable iv_to at lvl+1
  subroutine prolong(mg, lvl, iv, iv_to, add)
    use m_communication
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: lvl   !< Level to prolong from
    integer, intent(in)       :: iv    !< Source variable
    integer, intent(in)       :: iv_to !< Target variable
    logical, intent(in)       :: add   !< If true, add to current values
    integer                   :: i, id, dsize, nc

    if (lvl == mg%highest_lvl) error stop "cannot prolong highest level"
    if (lvl < mg%lowest_lvl) error stop "cannot prolong below lowest level"

    if (lvl >= mg%first_normal_lvl-1) then
       dsize            = (mg%box_size/2 + 2)**NDIM
       mg%buf(:)%i_send = 0
       mg%buf(:)%i_ix   = 0

       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          call prolong_set_buffer(mg, id, iv)
       end do

       mg%buf(:)%i_recv = mg%comm_prolong%n_recv(:, lvl) * dsize
       call sort_and_transfer_buffers(mg, dsize)
       mg%buf(:)%i_recv = 0
    end if

    nc = mg%box_size_lvl(lvl+1)
    do i = 1, size(mg%lvls(lvl+1)%my_ids)
       id = mg%lvls(lvl+1)%my_ids(i)
       call prolong_onto(mg, id, nc, iv, iv_to, add)
    end do
  end subroutine prolong

  subroutine prolong_set_buffer(mg, id, iv)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: iv
    integer                      :: i, hnc, dix(NDIM)
    integer                      :: i_c, c_id, c_rank, dsize

    hnc   = mg%box_size/2
    dsize = (hnc+2)**NDIM

    do i_c = 1, num_children
       c_id = mg%boxes(id)%children(i_c)
       if (c_id > no_box) then
          c_rank = mg%boxes(c_id)%rank
          if (c_rank /= mg%my_rank) then
             dix = get_child_offset(mg, c_id)
             i   = mg%buf(c_rank)%i_send

#if NDIM == 2
             mg%buf(c_rank)%send(i+1:i+dsize) = &
                  pack(mg%boxes(id)%cc(dix(1):dix(1)+hnc+1, &
                  dix(2):dix(2)+hnc+1, iv), .true.)
#elif NDIM == 3
             mg%buf(c_rank)%send(i+1:i+dsize) = &
                  pack(mg%boxes(id)%cc(dix(1):dix(1)+hnc+1, &
                  dix(2):dix(2)+hnc+1, dix(3):dix(3)+hnc+1, iv), .true.)
#endif

             mg%buf(c_rank)%i_send  = mg%buf(c_rank)%i_send + dsize
             i                      = mg%buf(c_rank)%i_ix
             mg%buf(c_rank)%ix(i+1) = c_id
             mg%buf(c_rank)%i_ix    = mg%buf(c_rank)%i_ix + 1
          end if
       end if
    end do
  end subroutine prolong_set_buffer

  subroutine prolong_onto(mg, id, nc, iv, iv_to, add)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: iv
    integer, intent(in)       :: iv_to
    logical, intent(in)       :: add
    integer                   :: IJK, hnc, p_id, p_rank, dix(NDIM), dsize
#if NDIM == 2
    real(dp)                  :: f0, flx, fhx, fly, fhy
#elif NDIM == 3
    real(dp)                  :: f0, flx, fhx, fly, fhy, flz, fhz
#endif
    real(dp)                  :: tmp(DTIMES(0:nc/2+1))

    hnc    = nc/2
    p_id   = mg%boxes(id)%parent
    p_rank = mg%boxes(p_id)%rank

    if (p_rank == mg%my_rank) then
       dix    = get_child_offset(mg, id)
#if NDIM == 2
       tmp = mg%boxes(p_id)%cc(dix(1):dix(1)+hnc+1, &
            dix(2):dix(2)+hnc+1, iv)
#elif NDIM == 3
       tmp = mg%boxes(p_id)%cc(dix(1):dix(1)+hnc+1, &
            dix(2):dix(2)+hnc+1, dix(3):dix(3)+hnc+1, iv)
#endif
    else
       dsize  = (hnc+2)**NDIM
       i = mg%buf(p_rank)%i_recv
       tmp = reshape(mg%buf(p_rank)%recv(i+1:i+dsize), [DTIMES(hnc+2)])
       mg%buf(p_rank)%i_recv = mg%buf(p_rank)%i_recv + dsize
    end if

    if (.not. add) then
       mg%boxes(id)%cc(DTIMES(:), iv_to) = 0.0_dp
    end if

#if NDIM == 2
    do j = 1, hnc
       do i = 1, hnc
          f0  = 0.5_dp * tmp(i, j)
          flx = 0.25_dp * tmp(i-1, j)
          fhx = 0.25_dp * tmp(i+1, j)
          fly = 0.25_dp * tmp(i, j-1)
          fhy = 0.25_dp * tmp(i, j+1)

          mg%boxes(id)%cc(2*i-1, 2*j-1, iv_to) = f0 + flx + fly + &
               mg%boxes(id)%cc(2*i-1, 2*j-1, iv_to)
          mg%boxes(id)%cc(2*i  , 2*j-1, iv_to) = f0 + fhx + fly + &
               mg%boxes(id)%cc(2*i  , 2*j-1, iv_to)
          mg%boxes(id)%cc(2*i-1, 2*j,   iv_to) = f0 + flx + fhy + &
               mg%boxes(id)%cc(2*i-1, 2*j,   iv_to)
          mg%boxes(id)%cc(2*i  , 2*j,   iv_to) = f0 + fhx + fhy + &
               mg%boxes(id)%cc(2*i  , 2*j,   iv_to)
       end do
    end do
#elif NDIM == 3
    do k = 1, hnc
       do j = 1, hnc
          do i = 1, hnc
             f0  = 0.25_dp * tmp(i, j, k)
             flx = 0.25_dp * tmp(i-1, j, k)
             fhx = 0.25_dp * tmp(i+1, j, k)
             fly = 0.25_dp * tmp(i, j-1, k)
             fhy = 0.25_dp * tmp(i, j+1, k)
             flz = 0.25_dp * tmp(i, j, k-1)
             fhz = 0.25_dp * tmp(i, j, k+1)

             mg%boxes(id)%cc(2*i-1, 2*j-1, 2*k-1, iv_to) = &
                  mg%boxes(id)%cc(2*i-1, 2*j-1, 2*k-1, iv_to) + &
                  f0 + flx + fly + flz
             mg%boxes(id)%cc(2*i, 2*j-1, 2*k-1, iv_to)   = &
                  mg%boxes(id)%cc(2*i, 2*j-1, 2*k-1, iv_to) + &
                  f0 + fhx + fly + flz
             mg%boxes(id)%cc(2*i-1, 2*j, 2*k-1, iv_to)   = &
                  mg%boxes(id)%cc(2*i-1, 2*j, 2*k-1, iv_to) + &
                  f0 + flx + fhy + flz
             mg%boxes(id)%cc(2*i, 2*j, 2*k-1, iv_to)      = &
                  mg%boxes(id)%cc(2*i, 2*j, 2*k-1, iv_to) + &
                  f0 + fhx + fhy + flz

             mg%boxes(id)%cc(2*i-1, 2*j-1, 2*k, iv_to) = &
                  mg%boxes(id)%cc(2*i-1, 2*j-1, 2*k, iv_to) + &
                  f0 + flx + fly + fhz
             mg%boxes(id)%cc(2*i, 2*j-1, 2*k, iv_to)   = &
                  mg%boxes(id)%cc(2*i, 2*j-1, 2*k, iv_to) + &
                  f0 + fhx + fly + fhz
             mg%boxes(id)%cc(2*i-1, 2*j, 2*k, iv_to)   = &
                  mg%boxes(id)%cc(2*i-1, 2*j, 2*k, iv_to) + &
                  f0 + flx + fhy + fhz
             mg%boxes(id)%cc(2*i, 2*j, 2*k, iv_to)      = &
                  mg%boxes(id)%cc(2*i, 2*j, 2*k, iv_to) + &
                  f0 + fhx + fhy + fhz
          end do
       end do
    end do
#endif
  end subroutine prolong_onto

end module m_prolong
