#include "../src/cpp_macros.h"
program test_one_level
  use mpi
#ifndef SINGLE_MODULE
  use m_octree_mg
#elif NDIM == 1
  use m_octree_mg_1d
#elif NDIM == 2
  use m_octree_mg_2d
#elif NDIM == 3
  use m_octree_mg_3d
#endif

  implicit none

  integer             :: n_its          = 10
  integer             :: n_modes(NDIM)  = 5
  integer             :: box_size
  integer             :: domain_size(NDIM)
  real(dp)            :: dr(NDIM)
  real(dp)            :: r_min(NDIM)    = 0.0_dp
  logical             :: periodic(NDIM) = .false.
  logical             :: fmg_cycle      = .true.
  integer             :: n_finer        = 0
  real(dp), parameter :: pi             = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: n, ierr, n_args
  real(dp)            :: t0, t1, t2, n_unknowns
  integer             :: i_sol, i_eps
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args < NDIM+1 .or. n_args > NDIM+3) then
     error stop "Usage: ./test_uniform_grid box_size domain_size(NDIM)" &
          // " [n_its] [FMG?]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(1+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  if (n_args > NDIM+1) then
     call get_command_argument(NDIM+2, arg_string)
     read(arg_string, *) n_its
  end if

  if (n_args > NDIM+2) then
     call get_command_argument(NDIM+3, arg_string)
     read(arg_string, *) fmg_cycle
  end if

  dr =  1.0_dp / domain_size

  mg%n_extra_vars = 2
  i_eps = mg_num_vars + 1
  i_sol = mg_num_vars + 2

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_laplacian
  mg%smoother_type = mg_smoother_gs

  ! This sets Dirichlet zero boundary conditions
  ! do n = 1, mg_num_neighbors
  !    mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
  !    mg%bc(n, mg_iphi)%bc_value = 1.0_dp
  ! end do

  ! This specifies a function for the boundary conditions
  do n = 1, mg_num_neighbors
     mg%bc(n, mg_iphi)%boundary_cond => my_boundary_condition
  end do

  call mg_set_methods(mg)
  call mg_comm_init(mg)
  t0 = mpi_wtime()
  call mg_build_rectangle(mg, domain_size, box_size, dr, r_min, &
       periodic, n_finer)
  call mg_load_balance(mg)
  t1 = mpi_wtime()
  call mg_allocate_storage(mg)
  t2 = mpi_wtime()

  if (mg%my_rank == 0) then
     print *, "mesh construction (s) ", t1-t0
     print *, "allocate storage (s)  ", t2-t1
  end if

  call set_solution(mg)
  call compute_rhs_and_reset(mg)

  call print_error(mg, 0)

  t0 = mpi_wtime()
  do n = 1, n_its
     if (fmg_cycle) then
        call mg_fas_fmg(mg, n > 1)
     else
        call mg_fas_vcycle(mg)
     end if
     call print_error(mg, n)
  end do
  t1 = mpi_wtime()

  if (mg%my_rank == 0) then
     if (fmg_cycle) then
        print *, "cycle type        FMG"
     else
        print *, "cycle type        V-cycle"
     end if
     print *, "n_cpu            ", mg%n_cpu
     print *, "problem_size     ", domain_size
     print *, "box_size         ", box_size
     print *, "n_iterations     ", n_its
     print *, "time/iteration   ", (t1-t0) / n_its
     print *, "total_time(s)    ", (t1-t0)
     n_unknowns = real(mg%n_boxes, dp) * mg%box_size**NDIM
     print *, "unknowns/microsec", 1e-6_dp * n_its * &
          n_unknowns / (t1-t0)
     print *, ""
  end if
  call mg_timers_show(mg)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

  real(dp) function solution(r)
    real(dp), intent(in) :: r(NDIM)
    solution = product(sin(2 * pi * n_modes * r))
  end function solution

  subroutine set_solution(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK
    real(dp)                  :: r(NDIM)

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(0, nc+1)
             r = mg%boxes(id)%r_min + ([IJK] - 0.5_dp) * mg%dr(:, lvl)
             mg%boxes(id)%cc(IJK, i_sol) = solution(r)
             mg%boxes(id)%cc(IJK, i_eps) = max(1.0_dp, 1.0_dp + r(1))
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

  subroutine print_error(mg, it)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: it
    integer                   :: n, nc, id, lvl, IJK, ierr
    real(dp)                  :: sol, val, err, max_err, res, max_res

    err = 0.0_dp
    res = 0.0_dp

    do lvl = mg%highest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(1, nc)
             sol = mg%boxes(id)%cc(IJK, i_sol)
             val = mg%boxes(id)%cc(IJK, mg_iphi)
             err = max(err, abs(val-sol))
             res = max(res, abs(mg%boxes(id)%cc(IJK, mg_ires)))
             ! print *, lvl, id, i, j, sol, val, abs(sol-val)
          end do; CLOSE_DO
       end do
    end do

    call mpi_reduce(err, max_err, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    call mpi_reduce(res, max_res, 1, MPI_DOUBLE, MPI_MAX, 0, &
         mpi_comm_world, ierr)
    if (mg%my_rank == 0) write(*, "(I4,A,E12.4,A,E12.4)") it, &
         "  max solution error", max_err, &
         "  max residual", max_res
  end subroutine print_error

  subroutine my_boundary_condition(box, nc, iv, nb, bc_type, bc)
    type(mg_box_t), intent(in) :: box
    integer, intent(in)        :: nc
    integer, intent(in)        :: iv         !< Index of variable
    integer, intent(in)        :: nb         !< Direction
    integer, intent(out)       :: bc_type    !< Type of b.c.
#if NDIM == 2
    real(dp), intent(out)      :: bc(nc)     !< Boundary values
    real(dp)                   :: x(nc, 2)
    integer                    :: i
#elif NDIM == 3
    real(dp), intent(out)      :: bc(nc, nc) !< Boundary values
    real(dp)                   :: x(nc, nc, 3)
    integer                    :: i, j
#endif

    call mg_get_face_coords(box, nb, nc, x)
    bc_type = mg_bc_dirichlet

#if NDIM == 2
    do i = 1, nc
       bc(i) = solution(x(i, :))
    end do
#elif NDIM == 3
    do j = 1, nc
       do i = 1, nc
          bc(i, j) = solution(x(i, j, :))
       end do
    end do
#endif
  end subroutine my_boundary_condition

end program test_one_level
