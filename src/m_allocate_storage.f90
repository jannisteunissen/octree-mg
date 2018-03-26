#include "cpp_macros.h"
module m_allocate_storage
  use m_data_structures

  implicit none
  private

  public :: allocate_storage
  public :: deallocate_storage

contains

  subroutine deallocate_storage(mg)
    type(mg_t), intent(inout) :: mg

    if (.not. mg%is_allocated) &
         error stop "deallocate_storage: tree is not allocated"

    deallocate(mg%boxes)
    deallocate(mg%buf)

    deallocate(mg%comm_restrict%n_send)
    deallocate(mg%comm_restrict%n_recv)

    deallocate(mg%comm_prolong%n_send)
    deallocate(mg%comm_prolong%n_recv)

    deallocate(mg%comm_ghostcell%n_send)
    deallocate(mg%comm_ghostcell%n_recv)

    mg%is_allocated = .false.
    mg%n_boxes      = 0
  end subroutine deallocate_storage

  subroutine allocate_storage(mg)
    use m_ghost_cells, only: ghost_cell_buffer_size
    use m_restrict, only: restrict_buffer_size
    use m_prolong, only: prolong_buffer_size
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl, nc
    integer                   :: n_send(0:mg%n_cpu-1, 3)
    integer                   :: n_recv(0:mg%n_cpu-1, 3)
    integer                   :: dsize(3)
    integer                   :: n_in, n_out, n_id

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(DTIMES(0:nc+1), n_var))

          mg%boxes(id)%cc(DTIMES(:), :) = 0.0_dp
       end do
    end do

    allocate(mg%buf(0:mg%n_cpu-1))

    call ghost_cell_buffer_size(mg, n_send(:, 1), &
         n_recv(:, 1), dsize(1))
    call restrict_buffer_size(mg, n_send(:, 2), &
         n_recv(:, 2), dsize(2))
    call prolong_buffer_size(mg, n_send(:, 3), &
         n_recv(:, 3), dsize(3))

    do i = 0, mg%n_cpu-1
       n_out = maxval(n_send(i, :) * dsize(:))
       n_in = maxval(n_recv(i, :) * dsize(:))
       n_id = maxval(n_send(i, :))
       allocate(mg%buf(i)%send(n_out))
       allocate(mg%buf(i)%recv(n_in))
       allocate(mg%buf(i)%ix(n_id))
    end do

    mg%is_allocated = .true.
  end subroutine allocate_storage

end module m_allocate_storage
