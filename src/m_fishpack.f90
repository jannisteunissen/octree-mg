module m_fishpack
  use m_data_structures

  implicit none
  private

  ! Public methods
  public :: fishpack_2d

  interface
     subroutine hstcrt(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd, &
          elmbda,f,idimf,pertrb,ierror,w)
       import
       real(dp), intent(in)    :: a, b
       integer, intent(in)     :: m, mbdcnd
       real(dp), intent(in)    :: bda(n), bdb(n)
       real(dp), intent(in)    :: c, d
       integer, intent(in)     :: n, nbdcnd
       real(dp), intent(in)    :: bdc(m), bdd(m)
       real(dp), intent(in)    :: elmbda
       real(dp), intent(inout) :: F(m, n)
       integer, intent(in)     :: idimf
       real(dp), intent(inout) :: w(*)
       real(dp), intent(inout) :: pertrb
       integer, intent(inout)  :: ierror
     end subroutine hstcrt
  end interface

contains

  subroutine fishpack_2d(nx, rhs, bndcnd, rmin, rmax)
    integer, intent(in)     :: nx(2)
    real(dp), intent(inout) :: rhs(nx(1), nx(2))
    type(bc_t), intent(in)  :: bndcnd(4)
    real(dp), intent(in)    :: rmin(2)
    real(dp), intent(in)    :: rmax(2)

    real(dp), allocatable :: w(:)
    real(dp)              :: pertrb
    integer               :: ierror
    integer               :: nwork
    integer               :: mbdcnd, nbdcnd

    if (bndcnd(neighb_lowx)%bc_type == bc_periodic) then
       mbdcnd = 0
    else if (bndcnd(neighb_lowx)%bc_type == bc_dirichlet) then
       if (bndcnd(neighb_highx)%bc_type == bc_dirichlet) then
          mbdcnd = 1
       else
          mbdcnd = 2
       end if
    else
       if (bndcnd(neighb_highx)%bc_type == bc_dirichlet) then
          mbdcnd = 4
       else
          mbdcnd = 3
       end if
    end if

    if (bndcnd(neighb_lowy)%bc_type == bc_periodic) then
       nbdcnd = 0
    else if (bndcnd(neighb_lowy)%bc_type == bc_dirichlet) then
       if (bndcnd(neighb_highy)%bc_type == bc_dirichlet) then
          nbdcnd = 1
       else
          nbdcnd = 2
       end if
    else
       if (bndcnd(neighb_highy)%bc_type == bc_dirichlet) then
          nbdcnd = 4
       else
          nbdcnd = 3
       end if
    end if

    nwork = 13 * nx(1) + 4 * nx(2) + &
         nx(1) * int(log(real(nx(2), dp))/log(2.0_dp))
    allocate(w(nwork))

    call hstcrt(rmin(1), rmax(1), nx(1), mbdcnd, bndcnd(1)%d, &
         bndcnd(2)%d, rmin(2), rmax(2), nx(2), nbdcnd, bndcnd(3)%d, &
         bndcnd(4)%d, 0.0_dp, rhs, nx(1), pertrb, ierror, w)

    if (ierror /= 0) then
       print *, "hstcrt returned ierror = ", ierror
       error stop
    end if
  end subroutine fishpack_2d

end module m_fishpack
