module m_data_structures

  implicit none
  public

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: no_box = 0
  integer, parameter :: physical_boundary = -1

  integer, parameter :: nb_offset_2d(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])

  integer, parameter :: child_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  integer, parameter :: child_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])

  !> Lists of blocks per refinement level
  type lvl_t
     integer, allocatable :: leaves(:)
     integer, allocatable :: parents(:)
     integer, allocatable :: ids(:)
     integer, allocatable :: my_leaves(:)
     integer, allocatable :: my_parents(:)
     integer, allocatable :: my_ids(:)
  end type lvl_t

  !> Box data structure
  type box_2d_t
     integer  :: rank
     integer  :: id
     integer  :: lvl
     integer  :: ix(2)
     integer  :: parent
     integer  :: children(4)
     integer  :: neighbors(4)
     real(dp) :: r_min(2)

     ! integer(int64)        :: morton
     real(dp), allocatable :: cc(:, :, :)
  end type box_2d_t

  type mg_2d_t
     integer                     :: n_cpu
     integer                     :: block_size(2) = [-1, -1]
     integer                     :: highest_lvl   = 0
     integer                     :: n_boxes
     real(dp), allocatable       :: dr(:)
     type(box_2d_t), allocatable :: boxes(:)
     type(lvl_t), allocatable    :: lvls(:)
  end type mg_2d_t

end module m_data_structures
