module m_allocate_storage
  use m_data_structures

  implicit none
  private

  public :: allocate_storage

contains

  subroutine allocate_storage(mg)
    use m_ghost_cells, only: ghost_cell_buffer_size
    use m_restrict, only: restrict_buffer_size
    use m_prolong, only: prolong_buffer_size
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: i, id, lvl, nc
    integer                      :: n_send(0:mg%n_cpu-1, 3)
    integer                      :: n_recv(0:mg%n_cpu-1, 3)
    integer                      :: dsize(3)
    integer                      :: n_in, n_out, n_id

    nc = mg%box_size

    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(0:nc+1, 0:nc+1, n_var))

          mg%boxes(id)%cc(:, :, :) = 0.0_dp
       end do
    end do

    call ghost_cell_buffer_size(mg, n_send(:, 1), &
         n_recv(:, 1), dsize(1))
    call restrict_buffer_size(mg, n_send(:, 2), &
         n_recv(:, 2), dsize(2))
    call prolong_buffer_size(mg, n_send(:, 3), &
         n_recv(:, 3), dsize(3))

    allocate(mg%buf(0:mg%n_cpu-1))

    do i = 0, mg%n_cpu-1
       n_out = maxval(n_send(i, :) * dsize(:))
       n_in = maxval(n_recv(i, :) * dsize(:))
       n_id = maxval(n_send(i, :))
       allocate(mg%buf(i)%send(n_out))
       allocate(mg%buf(i)%recv(n_in))
       allocate(mg%buf(i)%ix(n_id))
    end do

  end subroutine allocate_storage

end module m_allocate_storage
