module m_ghost_cells

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: fill_ghost_cells_lvl

contains

  subroutine fill_ghost_cells_lvl(tree, lvl)
    type(tree_2d_t)     :: tree
    integer, intent(in) :: lvl
    integer             :: i, id

    do i = 1, size(tree(lvl)%my_ids)
       id = tree(lvl)%my_ids(i)
       call fill_ghost_cells_box(tree%boxes, id)
    end do
  end subroutine fill_ghost_cells_lvl

  subroutine fill_ghost_cells_box(boxes, id)

    integer, intent(in) :: id
    integer :: nb

    do nb = 1, num_neighbors
       nb_id = boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          ! There is a neighbor
          nb_rank = boxes(nb_id)%rank
          if (nb_rank == my_rank) then
             call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
          else
             call sendrecv_from_nb(boxes(id))
          end if
          nb_dim     = a$D_neighb_dim(nb)
          lo(:)      = 1
          hi(:)      = boxes(id)%n_cell
          lo(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell + 1)
          hi(nb_dim) = lo(nb_dim)
          dnb        = a$D_neighb_offset([nb])
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
       else if (nb_id == af_no_box) then
          ! Refinement boundary
          call subr_rb(boxes, id, nb, iv)
       else
          ! Physical boundary
          call subr_bc(boxes(id), nb, iv, bc_type)
          call bc_to_gc(boxes(id), nb, iv, bc_type)
       end if
    end do
  end subroutine fill_ghost_cells_box


end module m_ghost_cells
