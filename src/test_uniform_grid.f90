program test_one_level
  use mpi
  use m_data_structures
  use m_build_tree
  use m_load_balance
  use m_ghost_cells
  use m_allocate_storage
  use m_restrict
  use m_communication
  use m_prolong

  implicit none

  integer, parameter :: block_size = 16
  integer, parameter :: domain_size = 2048
  real(dp), parameter :: dr = 1.0_dp / block_size

  integer :: lvl, ierr
  type(mg_2d_t) :: mg

  call comm_init(mg)
  call build_uniform_tree(mg, block_size, domain_size, dr)

  if (mg%my_rank == 0) then
     do lvl = 1, mg%highest_lvl
        print *, lvl, ":", size(mg%lvls(lvl)%ids)
     end do
  end if

  call load_balance(mg)
  call allocate_storage(mg)
  call set_initial_conditions(mg)

  do lvl = 1, mg%highest_lvl
     call fill_ghost_cells_lvl(mg, lvl)
  end do
  call print_error(mg)

  call prolong(mg, mg%highest_lvl-1)
  call fill_ghost_cells_lvl(mg, mg%highest_lvl)
  call print_error(mg)

  ! call multigrid_vcycle(mg)
  call mpi_finalize(ierr)

contains

  subroutine set_initial_conditions(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: n, id, lvl, i, j
    real(dp)                     :: r(2), sol

    do lvl = 1, mg%highest_lvl
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do j = 1, mg%box_size
             do i = 1, mg%box_size
                r = mg%boxes(id)%r_min + &
                     [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
                sol = product(sin(2 * acos(-1.0_dp) * r))
                mg%boxes(id)%cc(i, j, i_phi) = sol
             end do
          end do
       end do
    end do
  end subroutine set_initial_conditions

  subroutine print_error(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: n, id, lvl, i, j
    real(dp)                     :: r(2), sol, val, err

    err = 0.0_dp

    do lvl = 1, mg%highest_lvl
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do j = 0, mg%box_size+1
             do i = 1, mg%box_size
                r = mg%boxes(id)%r_min + &
                     [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
                sol = product(sin(2 * acos(-1.0_dp) * r))
                val = mg%boxes(id)%cc(i, j, i_phi)
                err = max(err, abs(val-sol))
                ! print *, mg%my_rank, lvl, i, j, r, sol, val
             end do
          end do
       end do
    end do

    print *, mg%my_rank, "max err", err
  end subroutine print_error

end program test_one_level
