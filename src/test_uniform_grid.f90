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
  use m_multigrid

  implicit none

  integer, parameter :: block_size = 8
  integer, parameter :: domain_size = 1024
  real(dp), parameter :: dr = 1.0_dp / block_size
  real(dp), parameter :: pi = acos(-1.0_dp)

  integer :: lvl, n, ierr
  type(mg_2d_t) :: mg

  mg%boundary_cond => my_bc

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

  do n = 1, 10
     call mg_fas_vcycle(mg, .true.)
     call print_error(mg)
  end do
  call mpi_finalize(ierr)

contains

  subroutine set_initial_conditions(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: n, id, lvl, i, j
    real(dp)                     :: r(2), sol

    do lvl = 1, mg%highest_lvl
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do j = 0, mg%box_size+1
             do i = 0, mg%box_size+1
                r = mg%boxes(id)%r_min + &
                     [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
                sol = product(sin(2 * pi * r))
                mg%boxes(id)%cc(i, j, i_phi) = sol
             end do
          end do

          call box_lpl(mg, id, i_rhs)
          mg%boxes(id)%cc(:, :, i_phi) = 0
       end do
    end do
  end subroutine set_initial_conditions

  subroutine print_error(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer                      :: n, id, lvl, i, j
    real(dp)                     :: r(2), sol, val, err

    err = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do j = 0, mg%box_size+1
             do i = 1, mg%box_size
                r = mg%boxes(id)%r_min + &
                     [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
                sol = product(sin(2 * pi * r))
                val = mg%boxes(id)%cc(i, j, i_phi)
                err = max(err, abs(val-sol))
                ! err = max(err, abs(mg%boxes(id)%cc(i, j, i_res)))
                ! print *, lvl, id, i, j, r, sol, val, abs(sol-val)
             end do
          end do
       end do
    end do

    print *, mg%my_rank, "max err", err
  end subroutine print_error

  subroutine my_bc(mg, id, nb, bc_type)
    type(mg_2d_t), intent(inout) :: mg
    integer, intent(in)          :: id
    integer, intent(in)          :: nb      !< Direction
    integer, intent(out)         :: bc_type !< Type of b.c.

    call set_bc_dirichlet_zero(mg, id, nb, bc_type)
  end subroutine my_bc

end program
