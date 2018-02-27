module m_load_balance
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: load_balance

contains

  subroutine load_balance(mg, myrank)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in) :: myrank
    integer :: i, id, lvl

    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)
          mg%boxes(id)%rank = modulo(i, mg%n_cpu)
       end do
       call update_lvl_info(mg, mg%lvls(lvl), myrank)
    end do

  end subroutine load_balance

  subroutine update_lvl_info(mg, lvl, myrank)
    type(mg_2d_t), intent(inout) :: mg
    type(lvl_t), intent(inout)   :: lvl
    integer, intent(in)          :: myrank

    lvl%my_ids = pack(lvl%ids, mg%boxes(lvl%ids)%rank == myrank)
    lvl%my_leaves = pack(lvl%leaves, mg%boxes(lvl%leaves)%rank == myrank)
    lvl%my_parents = pack(lvl%parents, mg%boxes(lvl%parents)%rank == myrank)
  end subroutine update_lvl_info

end module m_load_balance
