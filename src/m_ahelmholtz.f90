#include "cpp_macros.h"
!> Module which contains multigrid procedures for a Helmholtz operator of the
!> form: div(D grad(phi)) - lambda*phi = f, where D has a smooth spatial
!> variation and a component in each spatial direction
module m_ahelmholtz
  use m_data_structures

  implicit none
  private

  !> The lambda used for the Helmholtz equation (should be >= 0)
  real(dp), public, protected :: ahelmholtz_lambda = 0.0_dp

  public :: ahelmholtz_set_methods
  public :: ahelmholtz_set_lambda

contains

  subroutine ahelmholtz_set_methods(mg)
    type(mg_t), intent(inout) :: mg

    if (mg%n_extra_vars == 0 .and. mg%is_allocated) then
       error stop "ahelmholtz_set_methods: mg%n_extra_vars == 0"
    else
      mg%n_extra_vars = max(NDIM, mg%n_extra_vars)
    end if

    ! Use Neumann zero boundary conditions for the variable coefficient, since
    ! it is needed in ghost cells.
    mg%bc(:, mg_iveps1)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps1)%bc_value = 0.0_dp

#if NDIM > 1
    mg%bc(:, mg_iveps2)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps2)%bc_value = 0.0_dp
#endif

#if NDIM > 2
    mg%bc(:, mg_iveps3)%bc_type = mg_bc_neumann
    mg%bc(:, mg_iveps3)%bc_value = 0.0_dp
#endif

    select case (mg%geometry_type)
    case (mg_cartesian)
       mg%box_op => box_ahelmh

       select case (mg%smoother_type)
       case (mg_smoother_gs, mg_smoother_gsrb)
          mg%box_smoother => box_gs_ahelmh
       case default
          error stop "ahelmholtz_set_methods: unsupported smoother type"
       end select
    case default
       error stop "ahelmholtz_set_methods: unsupported geometry"
    end select

  end subroutine ahelmholtz_set_methods

  subroutine ahelmholtz_set_lambda(lambda)
    real(dp), intent(in) :: lambda

    if (lambda < 0) &
         error stop "ahelmholtz_set_lambda: lambda < 0 not allowed"

    ahelmholtz_lambda = lambda
  end subroutine ahelmholtz_set_lambda

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine box_gs_ahelmh(mg, id, nc, redblack_cntr)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: redblack_cntr !< Iteration counter
    integer                   :: IJK, i0, di
    logical                   :: redblack
    real(dp)                  :: idr2(2*NDIM), u(2*NDIM)
    real(dp)                  :: a0(2*NDIM), a(2*NDIM), c(2*NDIM)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*NDIM:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    i0  = 1

    redblack = (mg%smoother_type == mg_smoother_gsrb)
    if (redblack) then
       di = 2
    else
       di = 1
    end if

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if NDIM == 1
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, i_eps1 => mg_iveps1)
      if (redblack) i0 = 2 - iand(redblack_cntr, 1)

      do i = i0, nc, di
         a0(1:2) = cc(i, i_eps1)
         u(1:2) = cc(i-1:i+1:2, n)
         a(1:2) = cc(i-1:i+1:2, i_eps1)
         c(:)   = 2 * a0(:) * a(:) / (a0(:) + a(:)) * idr2

         cc(i, n) = &
              (sum(c(:) * u(:)) - cc(i, mg_irhs)) / &
              (sum(c(:)) + ahelmholtz_lambda)
      end do
    end associate
#elif NDIM == 2
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps1 => mg_iveps1, &
         i_eps2 => mg_iveps2)

      do j = 1, nc
         if (redblack) &
              i0 = 2 - iand(ieor(redblack_cntr, j), 1)

         do i = i0, nc, di
            a0(1:2)     = cc(i, j, i_eps1)
            a0(3:4)     = cc(i, j, i_eps2)
            u(1:2) = cc(i-1:i+1:2, j, n)
            a(1:2) = cc(i-1:i+1:2, j, i_eps1)
            u(3:4) = cc(i, j-1:j+1:2, n)
            a(3:4) = cc(i, j-1:j+1:2, i_eps2)
            c(:)   = 2 * a0(:) * a(:) / (a0(:) + a(:)) * idr2

            cc(i, j, n) = &
                 (sum(c(:) * u(:)) - cc(i, j, mg_irhs)) / &
                 (sum(c(:)) + ahelmholtz_lambda)
         end do
      end do
    end associate
#elif NDIM == 3
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps1 => mg_iveps1, &
         i_eps2 => mg_iveps2, &
         i_eps3 => mg_iveps3)

      do k = 1, nc
         do j = 1, nc
            if (redblack) &
                 i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, di
               a0(1:2)     = cc(i, j, k, i_eps1)
               a0(3:4)     = cc(i, j, k, i_eps2)
               a0(4:5)     = cc(i, j, k, i_eps3)
               u(1:2) = cc(i-1:i+1:2, j, k, n)
               a(1:2) = cc(i-1:i+1:2, j, k, i_eps1)
               u(3:4) = cc(i, j-1:j+1:2, k, n)
               a(3:4) = cc(i, j-1:j+1:2, k, i_eps2)
               u(5:6) = cc(i, j, k-1:k+1:2, n)
               a(5:6) = cc(i, j, k-1:k+1:2, i_eps3)
               c(:)   = 2 * a0(:) * a(:) / (a0(:) + a(:)) * idr2

               cc(i, j, k, n) = (sum(c(:) * u(:)) - &
                    cc(i, j, k, mg_irhs)) / &
                    (sum(c(:)) + ahelmholtz_lambda)
            end do
         end do
      end do
    end associate
#endif
  end subroutine box_gs_ahelmh

  !> Perform Helmholtz operator on a box
  subroutine box_ahelmh(mg, id, nc, i_out)
    type(mg_t), intent(inout) :: mg
    integer, intent(in)       :: id
    integer, intent(in)       :: nc
    integer, intent(in)       :: i_out !< Index of variable to store Helmholtz in
    integer                   :: IJK
    real(dp)                  :: idr2(2*NDIM), a0(2*NDIM), u0, u(2*NDIM), a(2*NDIM)

    ! Duplicate 1/dr^2 array to multiply neighbor values
    idr2(1:2*NDIM:2) = 1/mg%dr(:, mg%boxes(id)%lvl)**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)

#if NDIM == 1
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, i_eps1 => mg_iveps1)
      do i = 1, nc
         a0(1:2)     = cc(i, i_eps1)
         a(1:2) = cc(i-1:i+1:2, i_eps1)
         u0     = cc(i, n)
         u(1:2) = cc(i-1:i+1:2, n)

         cc(i, i_out) = sum(2 * idr2 * &
              a0(:)*a(:)/(a0(:) + a(:)) * (u(:) - u0)) - &
              ahelmholtz_lambda * u0
      end do
    end associate
#elif NDIM == 2
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps1 => mg_iveps1, &
         i_eps2 => mg_iveps2)
      do j = 1, nc
         do i = 1, nc
            a0(1:2)     = cc(i, j, i_eps1)
            a0(3:4)     = cc(i, j, i_eps2)
            a(1:2) = cc(i-1:i+1:2, j, i_eps1)
            a(3:4) = cc(i, j-1:j+1:2, i_eps2)
            u0     = cc(i, j, n)
            u(1:2) = cc(i-1:i+1:2, j, n)
            u(3:4) = cc(i, j-1:j+1:2, n)

            cc(i, j, i_out) = sum(2 * idr2 * &
                 a0(:)*a(:)/(a0(:) + a(:)) * (u(:) - u0)) - &
                 ahelmholtz_lambda * u0
         end do
      end do
    end associate
#elif NDIM == 3
    associate (cc => mg%boxes(id)%cc, n => mg_iphi, &
         i_eps1 => mg_iveps1, &
         i_eps2 => mg_iveps2, &
         i_eps3 => mg_iveps3)
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc
               u0 = cc(i, j, k, n)
               a0(1:2) = cc(i, j, k, i_eps1)
               a0(3:4) = cc(i, j, k, i_eps2)
               a0(5:6) = cc(i, j, k, i_eps3)
               u(1:2) = cc(i-1:i+1:2, j, k, n)
               u(3:4) = cc(i, j-1:j+1:2, k, n)
               u(5:6) = cc(i, j, k-1:k+1:2, n)
               a(1:2) = cc(i-1:i+1:2, j, k, i_eps1)
               a(3:4) = cc(i, j-1:j+1:2, k, i_eps2)
               a(5:6) = cc(i, j, k-1:k+1:2, i_eps3)

               cc(i, j, k, i_out) = sum(2 * idr2 * &
                    a0(:)*a(:)/(a0(:) + a(:)) * (u(:) - u0)) - &
                    ahelmholtz_lambda * u0
            end do
         end do
      end do
    end associate
#endif
  end subroutine box_ahelmh

end module m_ahelmholtz
