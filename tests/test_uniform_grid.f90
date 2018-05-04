#include "../src/cpp_macros.h"
program test_one_level
  use mpi
  use m_octree_mg

  implicit none

  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr(NDIM), r_min(NDIM) = 0.0_dp
  logical             :: periodic(NDIM)  = .false.
  integer             :: n_finer         = 0
  real(dp), parameter :: pi              = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: lvl, n, ierr, n_args
  real(dp)            :: t0, t1
  integer             :: i_sol, i_eps
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args /= NDIM+1) then
     error stop "Usage: ./test_uniform_grid box_size nx ny [nz]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(1+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  dr =  1.0_dp / domain_size

  mg%n_extra_vars = 2
  i_eps = mg_num_vars + 1
  i_sol = mg_num_vars + 2

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_vlaplacian
  mg%smoother_type = smoother_gs

  call mg_set_methods(mg)

  ! if (mg%geometry_type == mg_cylindrical .and. NDIM == 3) then
  !    periodic(2) = .true.
  ! end if

  call mg_comm_init(mg)
  call mg_build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call mg_load_balance(mg)

  if (mg%my_rank == 0) then
     print *, "First normal level", mg%first_normal_lvl
     do lvl = mg%lowest_lvl, mg%highest_lvl
        print *, lvl, ":", size(mg%lvls(lvl)%ids), mg%box_size_lvl(lvl)
     end do
  end if

  call mg_allocate_storage(mg)

  call set_solution(mg)

  do lvl = mg%lowest_lvl, mg%highest_lvl
     call mg_fill_ghost_cells_lvl(mg, lvl, mg_iphi)
  end do

  call compute_rhs_and_reset(mg)

  do n = 1, num_neighbors
     mg%bc(n)%bc_type = bc_dirichlet
     mg%bc(n)%bc_value = 0.0_dp
  end do

  call print_error(mg)

  t0 = mpi_wtime()
  do n = 1, 10
     call mg_fas_fmg(mg, n > 1)
     call print_error(mg)
  end do
  t1 = mpi_wtime()

  if (mg%my_rank == 0) then
     print *, "time, time per iteration:", t1-t0, (t1-t0) / 10
     print *, "unknowns/microsec", 1e-6_dp * 10 * product(domain_size) / (t1-t0)
  end if
  call timers_show(mg)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

  subroutine set_solution(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK
    real(dp)                  :: r(NDIM), sol

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(0, nc+1)
#if NDIM == 2
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp] * mg%dr(:, lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, i_sol) = sol
             mg%boxes(id)%cc(i, j, i_eps) = 1.0_dp
#elif NDIM == 3

             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp, k-0.5_dp] * mg%dr(:, lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, k, i_sol) = sol
             mg%boxes(id)%cc(i, j, k, i_eps) = max(1.0_dp, 1.0_dp + r(1))
#endif
          end do; CLOSE_DO
       end do
    end do
  end subroutine set_solution

  subroutine compute_rhs_and_reset(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)

          mg%boxes(id)%cc(DTIMES(:), mg_iphi) = &
               mg%boxes(id)%cc(DTIMES(:), i_sol)
          call mg%box_op(mg, id, nc, mg_irhs)
          mg%boxes(id)%cc(DTIMES(:), mg_iphi) = 0.0_dp
       end do
    end do
  end subroutine compute_rhs_and_reset

  subroutine print_error(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, nc, id, lvl, IJK, ierr
    real(dp)                  :: sol, val, err, max_err

    err = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             sol = mg%boxes(id)%cc(IJK, i_sol)
             val = mg%boxes(id)%cc(IJK, mg_iphi)
             err = max(err, abs(val-sol))
             ! err = max(err, abs(mg%boxes(id)%cc(IJK, mg_ires)))
             ! print *, lvl, id, i, j, sol, val, abs(sol-val)
          end do; CLOSE_DO
       end do
    end do

    call mpi_reduce(err, max_err, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) print *, "max err", max_err
  end subroutine print_error

end program test_one_level
