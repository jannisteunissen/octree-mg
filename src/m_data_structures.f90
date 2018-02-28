module m_data_structures

  implicit none
  public

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: n_var = 3
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3

  integer, parameter :: no_box = 0
  integer, parameter :: physical_boundary = -1

  ! Numbering of children (same location as **corners**)
  integer, parameter :: num_children = 4
  integer, parameter :: child_lowx_lowy = 1
  integer, parameter :: child_highx_lowy = 2
  integer, parameter :: child_lowx_highy = 3
  integer, parameter :: child_highx_highy = 4

  ! Neighboring indices for each child
  integer, parameter :: child_neighbs(2, 4) = reshape([1,3,2,3,1,4,2,4], [2,4])
  ! Index offset for each child
  integer, parameter :: child_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  ! Reverse child index in each direction
  integer, parameter :: child_rev(4, 2) = reshape([2,1,4,3,3,4,1,2], [4,2])
  ! Children adjacent to a neighbor
  integer, parameter :: child_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])
  ! Which children have a low index per dimension
  logical, parameter :: child_low(4, 2) = reshape([.true., .false., .true., &
       .false., .true., .true., .false., .false.], [4,2])

  ! Neighbor topology information
  integer, parameter :: num_neighbors = 4
  integer, parameter :: neighb_lowx = 1
  integer, parameter :: neighb_highx = 2
  integer, parameter :: neighb_lowy = 3
  integer, parameter :: neighb_highy = 4

  ! Index offsets of neighbors
  integer, parameter :: neighb_dix(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])
  ! Which neighbors have a lower index
  logical, parameter :: neighb_low(4) = [.true., .false., .true., .false.]
  ! Opposite of nb_low, but now as 0,1 integers
  integer, parameter :: neighb_high_01(4) = [0, 1, 0, 1]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: neighb_high_pm(4) = [-1, 1, -1, 1]

  ! Reverse neighbors
  integer, parameter :: neighb_rev(4) = [2, 1, 4, 3]
  ! Direction (dimension) for a neighbor
  integer, parameter :: neighb_dim(4) = [1, 1, 2, 2]

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

  type gc_buf_t
     integer               :: i_send
     integer               :: i_recv
     integer               :: i_id
     integer, allocatable  :: idbuf(:)
     real(dp), allocatable :: sendbuf(:)
     real(dp), allocatable :: recvbuf(:)
  end type gc_buf_t

  type mg_2d_t
     integer                     :: n_cpu = -1
     integer                     :: my_rank = -1
     integer                     :: box_size
     integer                     :: highest_lvl
     integer                     :: n_boxes
     real(dp), allocatable       :: dr(:)
     type(box_2d_t), allocatable :: boxes(:)
     type(lvl_t), allocatable    :: lvls(:)
     type(gc_buf_t), allocatable :: gc(:)
  end type mg_2d_t

contains

  !> Compute the 'child index' for a box with spatial index ix. With 'child
  !> index' we mean the index in the children(:) array of its parent.
  integer function ix_to_ichild(ix)
    integer, intent(in) :: ix(2) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
    ix_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
  end function ix_to_ichild

end module m_data_structures
