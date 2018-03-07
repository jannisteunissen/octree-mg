module m_load_balance
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: load_balance

contains

  subroutine load_balance(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer :: i, id, lvl
    integer :: work_left, my_work, i_cpu

    do lvl = 1, mg%highest_lvl
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
          mg%boxes(id)%rank = modulo(i, mg%n_cpu)
       end do
       call update_lvl_info(mg, mg%lvls(lvl))
    end do

  end subroutine load_balance

  subroutine update_lvl_info(mg, lvl)
    type(mg_2d_t), intent(inout) :: mg
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
