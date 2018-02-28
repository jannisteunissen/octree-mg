module m_allocate_storage
  use m_data_structures

  implicit none
  private

  public :: allocate_storage

contains

  subroutine allocate_storage(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: i, id, lvl, nc

    nc = mg%box_size

    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(0:nc+1, 0:nc+1, n_var))

          mg%boxes(id)%cc(:, :, :) = 0.0_dp
       end do
    end do

    call alloc_gc_buffers(mg)

  end subroutine allocate_storage

  subroutine alloc_gc_buffers(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: n_msg(0:mg%n_cpu-1, 1:mg%highest_lvl)
    integer                      :: i, id, lvl, nb, nb_id, nb_rank, n, m

    n_msg(:, :) = 0

    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)

          do nb = 1, num_neighbors
             nb_id = mg%boxes(id)%neighbors(nb)

             if (nb_id > no_box) then
                ! There is a neighbor
                nb_rank    = mg%boxes(nb_id)%rank
                if (nb_rank /= mg%my_rank) then
                   n_msg(nb_rank, lvl) = n_msg(nb_rank, lvl) + 1
                end  if
             end if
          end do
       end do
       ! print *, mg%my_rank, lvl, "msg", n_msg(:, lvl)
    end do

    allocate(mg%gc(0:mg%n_cpu-1))

    do i = 0, mg%n_cpu-1
       m = maxval(n_msg(i, :))
       n = m * mg%box_size
       allocate(mg%gc(i)%sendbuf(n))
       allocate(mg%gc(i)%recvbuf(n))
       allocate(mg%gc(i)%idbuf(m))
    end do
  end subroutine alloc_gc_buffers

end module m_allocate_storage
