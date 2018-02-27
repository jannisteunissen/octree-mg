program test_one_level
  use mpi
  use m_data_structures
  use m_build_tree
  use m_load_balance

  implicit none

  integer, parameter :: block_size = 4
  integer, parameter :: domain_size = 32
  real(dp), parameter :: dr = 1.0_dp / domain_size

  integer :: lvl, i, id, ierr
  integer :: n_cpu, myrank
  type(mg_2d_t) :: mg

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, n_cpu, ierr)

  mg%n_cpu = n_cpu
  call build_uniform_tree(mg, block_size, domain_size, dr)

  if (myrank == 0) then
     do lvl = 1, mg%highest_lvl
        print *, lvl, ":", size(mg%lvls(lvl)%ids)
        ! print *, "parents ", mg%lvls(lvl)%parents
        ! print *, "leaves ", mg%lvls(lvl)%leaves
     end do
  end if

  call load_balance(mg, myrank)
  call set_initial_conditions(mg)

  ! call multigrid_vcycle(mg)
  call mpi_finalize(ierr)

contains

  subroutine set_initial_conditions(mg)
    type(mg_2d_t), intent(inout) :: mg
    integer :: i, id, lvl

    do lvl = 1, mg%highest_lvl
       do i = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(i)
          allocate(mg%boxes(id)%cc(0:block_size+1, 0:block_size+1, 3))
       end do
    end do
  end subroutine set_initial_conditions

end program test_one_level
