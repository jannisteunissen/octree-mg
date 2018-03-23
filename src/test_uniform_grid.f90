#include "cpp_macros.h"
program test_one_level
  use mpi
  use m_octree_mg

  implicit none

  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr, r_min(NDIM) = 0.0_dp
  logical             :: periodic(NDIM) = .true.
  integer             :: n_finer = 0
  real(dp), parameter :: pi = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: lvl, n, ierr, n_args
  real(dp)            :: t0, t1
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

  dr                     =  1.0_dp / minval(domain_size)
  mg%boundary_cond       => my_bc
  mg%n_cycle_up          =  2
  mg%n_cycle_down        =  2
  mg%smoother_type       =  smoother_gsrb
  mg%residual_coarse_abs =  1e-10_dp
  mg%residual_coarse_rel =  1e-10_dp

  call comm_init(mg)
  call build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call load_balance(mg)

  if (mg%my_rank == 0) then
     print *, "First normal level", mg%first_normal_lvl
     do lvl = mg%lowest_lvl, mg%highest_lvl
        print *, lvl, ":", size(mg%lvls(lvl)%ids), mg%box_size_lvl(lvl)
     end do
  end if

  call allocate_storage(mg)
  call set_initial_conditions(mg)

  do n = 1, num_neighbors
     allocate(mg%bc(n)%d(mg%box_size))
     mg%bc(n)%d = 0.0_dp
     mg%bc(n)%bc_type = bc_dirichlet
  end do

  call print_error(mg)

  do lvl = mg%lowest_lvl, mg%highest_lvl
     call fill_ghost_cells_lvl(mg, lvl)
  end do

  t0 = mpi_wtime()
  do n = 1, 10
     ! call mg_fas_fmg(mg, n==10, n > 1)
     call mg_fas_vcycle(mg, .true.)
     call print_error(mg)
     ! print *, mg%boxes(1)%cc(1:box_size, 1:box_size, i_rhs)
     ! print *, mg%boxes(1)%cc(0:box_size+1, 2, i_phi)
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

  subroutine set_initial_conditions(mg)
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
                  [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, i_phi) = sol
#elif NDIM == 3

             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp, k-0.5_dp] * mg%dr(lvl)
             sol = product(sin(2 * pi * r))
             mg%boxes(id)%cc(i, j, k, i_phi) = sol
#endif
          end do; CLOSE_DO

          call box_lpl(mg, id, nc, i_rhs)
          mg%boxes(id)%cc(DTIMES(:), i_phi) = 0
       end do
    end do
  end subroutine set_initial_conditions

  subroutine print_error(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, nc, id, lvl, IJK, ierr
    real(dp)                  :: r(NDIM), sol, val, err, max_err

    err = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
#if NDIM == 2
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp] * mg%dr(lvl)
#elif NDIM == 3
             r = mg%boxes(id)%r_min + &
                  [i-0.5_dp, j-0.5_dp, k-0.5_dp] * mg%dr(lvl)
#endif
             sol = product(sin(2 * pi * r))
             val = mg%boxes(id)%cc(IJK, i_phi)
             err = max(err, abs(val-sol))
             ! err = max(err, abs(mg%boxes(id)%cc(IJK, i_res)))
             ! print *, lvl, id, i, j, r, sol, val, abs(sol-val)
          end do; CLOSE_DO
       end do
    end do

    call mpi_reduce(err, max_err, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) print *, "max err", max_err
  end subroutine print_error

  subroutine my_bc(mg, id, nc, nb, bc_type)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: nb      !< Direction
    integer, intent(out)      :: bc_type !< Type of b.c.

    call set_bc_dirichlet_zero(mg, id, nc, nb, bc_type)
  end subroutine my_bc

end program test_one_level
