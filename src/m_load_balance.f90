module m_load_balance
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: load_balance
  public :: load_balance_parents

contains

  subroutine load_balance(mg)
    type(mg_t), intent(inout) :: mg
    integer :: i, id, lvl, single_cpu_lvl
    integer :: work_left, my_work, i_cpu
    ! integer :: c_ids(4), c_ranks(4)

    single_cpu_lvl = max(mg%first_normal_lvl-1, mg%lowest_lvl)
    do lvl = mg%lowest_lvl, single_cpu_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = 0
       end do
    end do

    do lvl = single_cpu_lvl+1, mg%highest_lvl
       work_left = size(mg%lvls(lvl)%ids)
       my_work   = 0
       i_cpu     = 0

       do i = 1, size(mg%lvls(lvl)%ids)
          if ((mg%n_cpu - i_cpu - 1) * my_work >= work_left) then
             i_cpu   = i_cpu + 1
             my_work = 0
          end if

          my_work = my_work + 1
          work_left = work_left - 1

          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = i_cpu
          ! mg%boxes(id)%rank = modulo(i, mg%n_cpu)
       end do
    end do

    ! do lvl = mg%highest_lvl-1, 1, -1
    !    do i = 1, size(mg%lvls(lvl)%ids)
    !       id = mg%lvls(lvl)%ids(i)
    !       if (has_children(mg%boxes(id))) then
    !          c_ids = mg%boxes(id)%children
    !          c_ranks = mg%boxes(c_ids)%rank
    !          mg%boxes(id)%rank = most_popular(c_ranks)
    !       end if
    !    end do
    ! end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine load_balance

  subroutine load_balance_parents(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: i, id, lvl
    integer                   :: c_ids(num_children)
    integer                   :: c_ranks(num_children)

    do lvl = mg%highest_lvl-1, mg%lowest_lvl, -1
       do i = 1, size(mg%lvls(lvl)%parents)
          id = mg%lvls(lvl)%parents(i)

          c_ids = mg%boxes(id)%children
          c_ranks = mg%boxes(c_ids)%rank
          mg%boxes(id)%rank = most_popular(c_ranks)
       end do
    end do

    do lvl = mg%lowest_lvl, mg%highest_lvl
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine load_balance_parents

  pure integer function most_popular(list)
    integer, intent(in) :: list(:)
    integer             :: i, best_count, current_count

    best_count   = 0
    most_popular = -1

    do i = 1, size(list)
       current_count = count(list == list(i))

       if (current_count > best_count) then
          most_popular = list(i)
       end if
    end do

  end function most_popular

  subroutine update_lvl_info(mg, lvl)
    type(mg_t), intent(inout) :: mg
    type(lvl_t), intent(inout)   :: lvl

    lvl%my_ids = pack(lvl%ids, &
         mg%boxes(lvl%ids)%rank == mg%my_rank)
    lvl%my_leaves = pack(lvl%leaves, &
         mg%boxes(lvl%leaves)%rank == mg%my_rank)
    lvl%my_parents = pack(lvl%parents, &
         mg%boxes(lvl%parents)%rank == mg%my_rank)
    lvl%my_ref_bnds = pack(lvl%ref_bnds, &
         mg%boxes(lvl%ref_bnds)%rank == mg%my_rank)
  end subroutine update_lvl_info

end module m_load_balance
