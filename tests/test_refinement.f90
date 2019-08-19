#include "../src/cpp_macros.h"
program test_refinement
  use mpi
#ifndef SINGLE_MODULE
  use m_octree_mg
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
  logical             :: periodic(NDIM) = .false.
  logical             :: fmg_cycle      = .true.
  integer             :: lvl, n_levels
  real(dp), parameter :: pi             = acos(-1.0_dp)
  character(len=40)   :: arg_string
  integer             :: n, ierr, n_args
  real(dp)            :: t0, t1, n_unknowns
  integer             :: i_sol
  type(mg_t)          :: mg

  n_args = command_argument_count()
  if (n_args < NDIM+2 .or. n_args > NDIM+4) then
     error stop "Usage: ./test_refinement n_levels box_size nx ny [nz] [n_its] [FMG?]"
  end if

  call get_command_argument(1, arg_string)
  read(arg_string, *) n_levels

  call get_command_argument(2, arg_string)
  read(arg_string, *) box_size

  do n = 1, NDIM
     call get_command_argument(2+n, arg_string)
     read(arg_string, *) domain_size(n)
  end do

  if (n_args > NDIM+2) then
     call get_command_argument(NDIM+3, arg_string)
     read(arg_string, *) n_its
  end if

  if (n_args > NDIM+3) then
     call get_command_argument(NDIM+4, arg_string)
     read(arg_string, *) fmg_cycle
  end if

  dr =  1.0_dp / domain_size

  mg%n_extra_vars = 1
  i_sol = mg_num_vars + 1

  mg%geometry_type = mg_cartesian
  mg%operator_type = mg_laplacian
  mg%smoother_type = mg_smoother_gs

  do n = 1, mg_num_neighbors
     mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
     mg%bc(n, mg_iphi)%bc_value = 0.0_dp
  end do

  call mg_set_methods(mg)
  call mg_comm_init(mg)

  call build_amr_tree(mg, n_levels, domain_size, box_size, dr, periodic)

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
     do lvl = 1, mg%highest_lvl
        write(*, '(A,I0,A,I0,A,I0,A)') " lvl_", lvl, ": ", &
             size(mg%lvls(lvl)%ids), " boxes, ", &
             size(mg%lvls(lvl)%leaves), " leaves"
     end do

     if (fmg_cycle) then
        print *, "cycle type        FMG"
     else
        print *, "cycle type        V-cycle"
     end if
     print *, "n_cpu            ", mg%n_cpu
     print *, "coarse grid      ", domain_size
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

  subroutine set_solution(mg)
    type(mg_t), intent(inout) :: mg
    integer                   :: n, id, lvl, nc, IJK
    real(dp)                  :: r(NDIM), sol

    do lvl = mg%lowest_lvl, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_ids)
          id = mg%lvls(lvl)%my_ids(n)
          do KJI_DO(0, nc+1)
             r = mg%boxes(id)%r_min + ([IJK] - 0.5_dp) * mg%dr(:, lvl)
             sol = product(sin(2 * pi * n_modes * r))
             mg%boxes(id)%cc(IJK, i_sol) = sol
          end do; CLOSE_DO
       end do
    end do

    ! Set ghost cells for the solution. Because they depend on coarse-grid
    ! values, first restrict the solution.
    call mg_restrict(mg, i_sol)
    call mg_fill_ghost_cells(mg, i_sol)
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
    real(dp)                  :: sol, val, err, max_err

    err = 0.0_dp

    do lvl = 1, mg%highest_lvl
       nc = mg%box_size_lvl(lvl)
       do n = 1, size(mg%lvls(lvl)%my_leaves)
          id = mg%lvls(lvl)%my_leaves(n)
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
    if (mg%my_rank == 0) print *, it, "max err", max_err
  end subroutine print_error

  subroutine build_amr_tree(mg, n_amr_levels, lvl1_size, box_size, dr, periodic)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: n_amr_levels
    integer, intent(in)       :: lvl1_size(NDIM)
    integer, intent(in)       :: box_size
    real(dp), intent(in)      :: dr(NDIM)
    logical, intent(in)       :: periodic(NDIM)
    integer                   :: lvl, i, id
    integer                   :: n_finer
    real(dp)                  :: r_min(NDIM), domain_len(NDIM)
    real(dp)                  :: r0(NDIM), r1(NDIM), box_center(NDIM)

    ! Estimate for the number of boxes on refined levels (can be large)
    n_finer    = n_amr_levels * product(lvl1_size / box_size) + 1000
    r_min      = 0.0_dp
    domain_len = lvl1_size * dr

    call mg_build_rectangle(mg, lvl1_size, box_size, dr, &
         r_min, periodic, n_finer)

    ! Add refinement
    do lvl = 1, n_amr_levels-1
       do i = 1, size(mg%lvls(lvl)%ids)
          id = mg%lvls(lvl)%ids(i)

          ! Refine center region
          r0 = 0.5_dp * domain_len - domain_len * 0.5**(lvl+1)
          r1 = 0.5_dp * domain_len + domain_len * 0.5**(lvl+1)

          box_center = mg%boxes(id)%r_min + 0.5_dp * box_size * mg%boxes(id)%dr

          if (all(box_center >= r0 .and. box_center <= r1)) then
             call mg_add_children(mg, id)
          end if
       end do

       call mg_set_leaves_parents(mg%boxes, mg%lvls(lvl))
       call mg_set_next_level_ids(mg, lvl)
       call mg_set_neighbors_lvl(mg, lvl+1)
    end do

    call mg_set_leaves_parents(mg%boxes, mg%lvls(n_amr_levels))

    mg%highest_lvl = n_amr_levels

    ! Store boxes with refinement boundaries (from the coarse side)
    do lvl = 1, mg%highest_lvl
       call mg_set_refinement_boundaries(mg%boxes, mg%lvls(lvl))
    end do

    ! Assign boxes to MPI processes
    call mg_load_balance(mg)

    ! Allocate storage for boxes owned by this process
    call mg_allocate_storage(mg)

  end subroutine build_amr_tree

end program test_refinement
