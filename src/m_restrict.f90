module m_restrict

  implicit none
  private

  ! Public methods
  public :: 

contains

  subroutine restrict(mg, lvl)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in) :: lvl

    do i = 1, size(mg%lvls(lvl)%my_ids)
       id = mg%lvls(lvl)%my_ids(i)
       call restrict_or_buffer_box(mg, id)
    end do
  end subroutine restrict_to_lvl

  subroutine restrict_or_buffer_box(mg, id)
    real(dp) :: tmp(mg%box_size/2, mg%box_size/2)

    do j = 1, mg%box_size/2
       do i = 1, mg%box_size/2
          tmp(i, j) = sum(mg%boxes(id)%cc(2*i-1:2*i, 2*j-1:2*j, i_phi))
       end do
    end do

    p_id = mg%boxes(id)%parent
    if (mg%boxes(p_id)%rank == mg%my_rank) then
       dix = ...
       mg%boxes(p_id)%cc(dix(1)+1:dix(1):hnc, dix(1)+1:dix(1):hnc, i_phi)
    else
       ! Buffer
    end if
  end subroutine restrict_or_buffer_box


end module m_restrict
