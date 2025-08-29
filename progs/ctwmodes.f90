! version 0.7.1 (007b): Minor update on 007 for publishing the code.
!   - Version numbering scheme has been changed.
!   - The new format of the Ne file now requires the number
!     datapoints (km+1) at the top. (Before, it required km.)
!
! version 007: Cleaned-up version of 006c.
!   - The module to define the finite-difference scheme is included
!     in this source file.
!   - The format of Grid.txt is improved (and incompatible
!     with previous version). See read_grid().
!
! version 006c: Based on 006. New finite-difference scheme
!    to use "diagonal" slopes and "boundary points".
!    The definition of the matrices is separated out to a module
!    diagonal_slopes_sym.f90 to organize the transition to the new
!    finite-difference scheme.
!    The original "steppy" scheme is found in 006b.
!
! version 006: Same as 005a but output includes im and km.
!   Former name was Eigenmode2d-006.f90 .
!
! version 005a: Same as 005 but the eigenvalues are 1/c.
!
! version 005:
!
! version 004: Place N(z) at the edges of gridcells
!   to put it inside the derivative as (F_z/N^2)_z .
!
! version 003: Introduces the other wall to support eastern
!   bounday Kelvin waves.
!
!---Solve an eigenvalue problem in a x-z plane---
!
!        N2e[0]   N2e[1]         N2e[im-1]   N2e[im]
!   F[0]--+--F[1]--+--F[2]-- . . . +---F[im]---+---F[im+1]
!
! The matrix form is
!   A x = (f/c)B x
! Both A are B symmetric;
!   B is positive semi-definite with a large kernel (rank(B) = km);
!   A is negative semi-definite when rigid-lid (barotropic mode is the kernel)
!      and negative definite when free-surface.
! We use the general solver ?GGEV of lapack
!   https://www.netlib.org/lapack/lug/node35.html
!   https://www.netlib.org/lapack/explore-html/da/da1/group__ggev.html
! or switch to the symmetric-definite solver ?SYGV
!   https://www.netlib.org/lapack/lug/node48.html
!   https://www.netlib.org/lapack/explore-html/d1/d39/group__hegv.html
! when free-surface.
!---------------------------------------
! This source file consists of
!   - module define_matrices_diag (formerly "diagonal_slopes_sym")
!   - program main
!
! 1. The main program determines the basic parameters,
!      kocn, dx, dz, N2e, f0, gravit
!   where kocn(i) is the k index of the deepest ocean grid cell for each i.
!
! 2. The main program then calls the subroutines from the module
!   to calculate the matrices A and B:
!      call calc_ocn(iocn, mask, &
!                    kocn)
!      call calc_num(num, &
!                    iocn, kocn)
!      M = maxval(num)             !-> # of active grid cells.
!      call eqs_all(A, B, l, &
!                   iocn, kocn, num, dx, dz, N2e, f0, M, gravit=gravit)
!
! 3. The main program then solves the matrix eigenvalue problem:
!     Input: A, B, f0
!     Output: VR(1:M, n) ... eigenvector for mode n.
!             alphaR(n) ... c, charateristic speed for mode n.
!             alphaI(n) ... imag part of eigenvalue for mode n.
!         where n=1,...,M.
! ==========================================================
!
module define_matrices_diag
!
! Same as diagonal_slopes.f90 but produces symmetrix A and B.
!
! Set the boundary conditions along the slope and bottom.
! We place special points right on the bottom-slope boundary.
!
! Where the slope is a concave step, the center of the box falls
! on the diagonal and that is the boundary point:
!   +
!   |╲  :
!   | Ⓒ :
!   |  ╲:
!   +−−−+
! The F value is stored at gridpoint (i,k).
!
! Where the "slope" is vertical, the special boundary is to the left
! of the center of the box:
!           |
!           +
!           |
! c(i-1,k)  ◯  c(i,k)
!           |
!           +
!           |
! The F value is stored at gridpoint (i-1,k).
!
! Where the "slope" is horizonal, the boundary point is below
! the center of the box:
!
!      c(i,k)
!
! -+-- ◯ --+-
!
!      c(i,k+1)
!
! The F value is stored at gridpoint (i,k+1).
!
! An isolated convex vertex
!   −−−+−−−+
!      :   |
!      :   |
!      +   +
!      :   |
!      :   |
! will be "shaved off", i.e., converted to a diagonal slope,
! the first kind discussed above.
!
! special characters
! ┖───
!   ╳
! ╱
!

  implicit NONE

  private
  public:: calc_ocn, calc_num, eqs_all

contains

! Helper subroutine
!   Just to check if (gridpoint num) == (equation num) for the "central"
!   gridpoint
subroutine eqnum(nc,l,i,k)
  use iso_fortran_env, only: ERROR_UNIT
  implicit NONE
  integer, intent(in):: nc, l, i, k
  if (nc /= l) then
    write(ERROR_UNIT,*) "nc, l, i, k = ", nc, l, i, k
  end if
end subroutine eqnum


! Define all the equations.
! Output:
!   A(:,:), B(:,:)
!   l: total equation number
!
! Input:
!   gravit: gravity of Earth; if gravit==NaN, rigid-lid is assumed.
!
! The free-surface factors, sfOne and sfTwo,
!   are defined and used only here.
!
subroutine eqs_all(A, B, l, iocn, kocn, num, dx, dz, N2e, f0, M, gravit)
  use ieee_arithmetic, only: ieee_is_nan
  use ctwmodes_pars, only: double
  implicit NONE
  real(double), intent(out),allocatable:: A(:,:), B(:,:)
  integer,      intent(out)            :: l ! equation number
  real(double), intent(in) :: dx(0:), dz(0:)
  integer,      intent(in) :: iocn(0:), kocn(0:), num(0:,0:)
  real(double), intent(in) :: N2e(0:), f0
  integer, intent(in):: M
  real(double), intent(in):: gravit

  ! (sfOne,sfTwo) = (1/2, 1/2 f^2/g) if freesurface and k = 1.
  ! (sfOne,sfTwo) = (1,   0)   otherwise
  real(double), allocatable:: sfOne(:), sfTwo(:)
  real(double):: fsq
  integer:: km

  fsq = f0*f0

  km = size(dz) - 2
  allocate(sfOne(1:km), sfTwo(1:km))
  sfOne(:) = 1
  sfTwo(:) = 0
  if (.not. ieee_is_nan(gravit)) then
    sfOne(1) = 1.d0/2
    sfTwo(1) = (fsq/gravit)/2
  end if

  allocate( A(M,M) )
  allocate( B(M,M) )
  A(:,:) = 0.0d0
  B(:,:) = 0.0d0
  l = 0 ! will be incremented by the following
  call eqs_interior(     A,    l, iocn, kocn, num, dx, dz, N2e, fsq, &
      sfOne, sfTwo)
  call eqs_surface_right(A,    l, iocn, kocn, num, dx, dz, N2e, fsq, &
      sfOne, sfTwo)
  call eqs_slope_bottom( A, B, l, iocn, kocn, num, dx, dz, N2e, fsq, &
      sfOne, sfTwo)
end subroutine eqs_all


! Gnerates the equations for the interior;
! the inerior points adjacent to the slope-bottom boundary is
! But, diagonal-slope boundary cells are NOT included in the "interior";
! they are boundary points.
! See eqs_slope_bottom().
!
! On input, l is the last equation number.
! l is incremented for each equation this subroutine "generates".
!
! Multiply dx(i) * dz(k) to make A symmetric
!
! free surface or rigid lid depending on sfOne and sfTwo
subroutine eqs_interior(A, l, iocn, kocn, num, dx, dz, N2e, fsq, &
    sfOne, sfTwo)
  use iso_fortran_env, only: ERROR_UNIT
  use ctwmodes_pars, only: double
  implicit NONE
  real(double), intent(inout):: A(:,:)
  integer,      intent(inout):: l ! equation number
  real(double), intent(in)   :: dx(0:), dz(0:)
  integer,      intent(in)   :: iocn(0:), kocn(0:), num(0:,0:)
  real(double), intent(in)   :: N2e(0:), fsq
  real(double), intent(in)   :: sfOne(:), sfTwo(:)
  real(double):: dxl, dxr, dzu, dzd
  integer:: i, k, im, km, nc, nl, nr, nu, nd
  im = size(dx) - 2
  km = size(dz) - 2
  do k = 1, km
    do i = iocn(k), im
      if (i==iocn(k) .and. k==kocn(i)) cycle ! diagonal is not inerior
      l = l + 1
      nc = num(i,k)
      nl = num(i-1,k)
      nr = num(i+1,k)
      nu = num(i,k-1)
      nd = num(i,k+1)
      if (i==iocn(k)) then ! wall
        dxl = dx(i)/2
      else
        dxl = (dx(i-1) + dx(i))/2
      endif
      dxr = (dx(i) + dx(i+1))/2
      if (k==kocn(i)) then ! bottom
        dzd = dz(k)/2
      else
        dzd = (dz(k) + dz(k+1))/2
      endif
      dzu = (dz(k-1) + dz(k))/2
      A(l,nc) = - dz(k) / dxl &
                - dz(k) / dxr &
                - dx(i) / dzu * fsq / N2e(k-1) * sfOne(k)&
                - dx(i) / dzd * fsq / N2e(k) &
                - dx(i) / 2 * sfTwo(k) ! free surface
      A(l,nl) = dz(k) / dxl
      A(l,nr) = dz(k) / dxr
      A(l,nu) = dx(i) / dzu * fsq / N2e(k-1) * sfOne(k)&
                - dx(i) / 2 * sfTwo(k) ! free surface
      A(l,nd) = dx(i) / dzd * fsq / N2e(k)
      call eqnum(nc,l,i,k)
    end do
  end do
end subroutine eqs_interior


! Gnerates the equations for the surface and right boundaries.
! Rigid-lid only.
! p_x = 0 at x = Lx only.
!
! On input, l is the last equation number.
! l is incremented for each equation this subroutine "generates".
!
! Multiply dx(i) or dz(k) to make A symmetric
!
! free surface or rigid lid depending on sfOne and sfTwo
subroutine eqs_surface_right(A, l, iocn, kocn, num, dx, dz, N2e, fsq, &
    sfOne, sfTwo)
  use iso_fortran_env, only: ERROR_UNIT
  use ctwmodes_pars, only: double
  implicit NONE
  real(double), intent(inout):: A(:,:)
  integer,      intent(inout):: l ! equation number
  real(double), intent(in)   :: dx(0:), dz(0:)
  integer,      intent(in)   :: iocn(0:), kocn(0:), num(0:,0:)
  real(double), intent(in)   :: N2e(0:), fsq
  real(double), intent(in)   :: sfOne(:), sfTwo(:)
  real(double):: dxc, dzc
  integer:: im, km, i, k, nc, nd, nl
  im = size(dx,1) - 2
  km = size(dz,1) - 2

! Boundary Condition at the Surface (z = 0)
! free surface or rigid lid depending on sfOne and sfTwo
  do i = iocn(1), im
    l = l + 1
    nc = num(i,0)
    nd = num(i,1)
    dzc = (dz(0) + dz(1))/2
    A(l,nc) = -dx(i) / dzc * fsq / N2e(0) * sfOne(1) &
              -dx(i) / 2 * sfTwo(1)
    A(l,nd) =  dx(i) / dzc * fsq / N2e(0) * sfOne(1) &
              -dx(i) / 2 * sfTwo(1)
    call eqnum(nc,l,i,k)
  end do

! Boundary condition offshore (x = L):
!  p_x = 0
  do k = 1, kocn(im)
    l = l + 1
    nc = num(im+1,k)
    nl = num(im,  k)
    dxc = (dx(im) + dx(im+1))/2
    A(l,nc) = -dz(k) / dxc
    A(l,nl) =  dz(k) / dxc
    call eqnum(nc,l,i,k)
  end do

end subroutine eqs_surface_right


! Gnerates the equations for the slope-bottom boundary condition.
! On input, l is the last equation number.
! l is incremented for each equation this subroutine "generates".
!
! Multiply dx(i) or dz(k) to make A and B symmetric
subroutine eqs_slope_bottom(A, B, l, iocn, kocn, num, dx, dz, N2e, fsq, &
    sfOne, sfTwo)
  use iso_fortran_env, only: ERROR_UNIT
  use ctwmodes_pars, only: double
  implicit NONE
  real(double), intent(inout):: A(:,:), B(:,:)
  integer,      intent(inout):: l ! equation number
  real(double), intent(in)   :: dx(0:), dz(0:)
  integer,      intent(in)   :: iocn(0:), kocn(0:), num(0:,0:)
  real(double), intent(in)   :: N2e(0:), fsq
  real(double), intent(in)   :: sfOne(:), sfTwo(:)
  character(*), parameter:: fmt&
      = "(1X,'l nc =',I8,I8,'; i k =',I5,I5,A)"
  real(double):: dxr, dzu, hx
  integer:: i, k, no, nr, nu, nc, im, km
  im = size(dx) - 2
  km = size(dz) - 2
  ! Scan the cell-center gridpoints adjacent to the slope-bottom.
  k = 1
  i = iocn(k)
  do while (i <= im)
    l = l + 1 ! next equation
    no = num(i,k) !! ocean point
    nr = num(i+1,k) !! to its right
    nu = num(i,k-1) !! to its above
    if (i == iocn(k) .and. k == kocn(i)) then ! next to bottom and wall
      nc = no ! bndry point is ocean point
      hx = dz(k)/dx(i) ! slope
      dxr = (dx(i) + dx(i+1))/2 ! right side
      dzu = (dz(k) + dz(k-1))/2 ! above
      !         - dz(k) * (- 1/dxr - 1/dzu * fsq/N2e(k-1) / hx
      A(l,nc) = - dz(k) / dxr &
                - dx(i) / dzu * fsq / N2e(k-1) * sfOne(k)&
                - dx(i) / 2 * sfTwo(k) ! free surface
      A(l,nr) = dz(k) / dxr
      !         dz(k)/dzu * fsq/N2e(k-1) / hx
      A(l,nu) = dx(i) / dzu * fsq / N2e(k-1) * sfOne(k)&
                - dx(i) / 2 * sfTwo(k) ! free surface
      B(l,nc) = dz(k)
      write(*,fmt) l,nc, i,k,": slope"
      !write(*,*) "dx^2/dz^2, N^2/f^2 = ", (dxr/dzu)**2, N2e(k-1)/(fsq)
    else if (i == iocn(k)) then ! vertical wall
      nc = num(i-1,k) ! bndry point is to the left
      dxr = dx(i)/2 ! right of the bndry point nc
      A(l,nc) = - dz(k)/dxr
      A(l,no) =   dz(k)/dxr
      B(l,nc) = dz(k)
      write(*,fmt) l,nc, i,k,": vertical"
    else if (k == kocn(i)) then ! flat bottom
      nc = num(i,k+1) ! bndry point is downward
      dzu = dz(k)/2 ! above the bottom point
      A(l,nc) = - dx(i)/dzu * fsq/N2e(k)
      A(l,no) =   dx(i)/dzu * fsq/N2e(k)
      write(*,fmt) l,nc, i,k,": horizontal"
    else
      write(ERROR_UNIT,*)&
          "Something is wrong: Point must be adjacent to bndry."
      error stop 9
    end if
    call eqnum(nc,l,i,k)

    ! move to next boundary point, following the bottom
    if (k == kocn(i)) then ! already at the bottom.
      i = i + 1 ! move right
      if (k < kocn(i)) then ! if not at the bottom,
        k = k + 1 ! move down
      end if
    else
      k = k + 1 ! move down
    end if
  end do
end subroutine eqs_slope_bottom


! ### OBSOLETE ###
! ### This subroutine calculates iocn and kocn from mask
! ### but we now calculate mask from kocn with calc_mask()
! ### and then calculate iocn from mask with calc_iocn().
!
subroutine calc_ocn_sub(iocn, kocn, mask)
  use iso_fortran_env, only: ERROR_UNIT
  implicit NONE
  integer, intent(out) :: iocn(0:), kocn(0:)
  integer, intent(in) :: mask(0:, 0:)
  integer:: im, km, i, k
  !###>>>
  write(ERROR_UNIT,*) "Don't use this. Calculate mask from kocn."
  error stop 1
  !<<<###
  im = size(mask,1) - 2
  km = size(mask,2) - 2
  iocn(:) = - huge(1)
  kocn(:) = - huge(1)
  do k = 1, km
    i = 1
    do while (i <= im)
      if (mask(i,k) == 1) exit
      i = i + 1
    end do
    if (i > im) then
      write(ERROR_UNIT,*) "iocn not found at k =", k
      error stop 3
    end if
    iocn(k) = i
  end do
  iocn(0) = iocn(1)

  do i = 1, im
    k = km
    do while (k >= 1)
      if (mask(i,k) == 1) exit
      k = k - 1
    end do
    if (k < 1) then
      write(ERROR_UNIT,*) "kocn not found at i =", i
      error stop 5
    end if
    kocn(i) = k
  end do
  kocn(im+1) = kocn(im)
end subroutine calc_ocn_sub


! "flag" is 1 in the lobes; "mask" is 1 only in the ocean.
! "flag" is real; "mask" is integer.
! Otherwise they are the same. In particular
!
!    mask(1:im,1:km) == flag(1:im,1:km)
!
subroutine calc_mask(mask, kocn)
  implicit NONE
  integer, intent(out):: mask(0:, 0:)
  integer, intent(out):: kocn(0:)
  integer:: i, k, im
  im = size(mask,1) - 2
  mask = 0
  do i = 1, im
    do k = 1, kocn(i)
      mask(i,k) = 1
    end do
  end do
end subroutine calc_mask

! calculate iocn from mask
subroutine calc_iocn(iocn, mask)
  use iso_fortran_env, only: ERROR_UNIT
  implicit NONE
  integer, intent(out):: iocn(0:)
  integer, intent(in) :: mask(0:, 0:)
  integer:: im, km, i, k

  im = size(mask,1) - 2
  km = size(mask,2) - 2
  iocn(:) = - huge(1)

  do k = 1, km
    i = 1
    do while (i <= im)
      if (mask(i,k) == 1) exit
      i = i + 1
    end do
    if (i > im) then
      write(ERROR_UNIT,*) "iocn not found at k =", k
      error stop 3
    end if
    iocn(k) = i
    !!>>>>>>>>>
    write(*,*) "calc_iocn: i, k, iocn(k) = ", i, k, iocn(k)
    !!<<<<<<<<
  end do
  iocn(0) = iocn(1)
end subroutine calc_iocn

! Follow the slope-bottom to look for isolated convex land cell.
! If we find one, we convert it to a "corner" cell by moving
!   in there and changing mask(i,k), kocn(i), and iocn(k)
!   consistently. Then, we move on to follow the slope-bottom cell.
subroutine shave_convex(iocn, kocn, mask, nshaved)
  integer, intent(inout):: iocn(0:), kocn(0:)
  integer, intent(inout):: mask(0:, 0:)
  integer, intent(out)  :: nshaved
  integer:: k, i, im
  im = size(mask,1) - 2
  nshaved = 0
  k = 1
  i = iocn(k)
  do while (i <= im)
    ! move, following the bottom
    if (k == kocn(i)) then ! already at the bottom.
      i = i + 1 ! move right
      if (k < kocn(i)) then ! not bottom
        k = k + 1 ! move down
        ! If the land cell to the left is isolated convex
        if (i-1 /= iocn(k-1) .and. k < kocn(i)) then
          i = i - 1 ! move back leftward under ground
          mask(i,k) = 1 ! land -> ocean
          kocn(i)   = k ! push bottom downward
          iocn(k)   = i ! push wall leftward
          write(*,"(1X,A,I4,',',I4,A)") "Shaved: convex virtex (i,k) = ("&
              , i,k, ")"
          nshaved = nshaved + 1
        end if
      end if
    else
      k = k + 1 ! move down
    end if
  end do
end subroutine shave_convex


! Calculate mask and iocn from kocn.
!
! But, if land includes isolated protrusions (convex verticies),
!   they will be shaved off and kocn, iocn, and mask will be revised.
subroutine calc_ocn(iocn, mask, kocn)
  use iso_fortran_env, only: ERROR_UNIT
  implicit NONE
  integer, intent(out)  :: iocn(0:)
  integer, intent(out)  :: mask(0:, 0:)
  integer, intent(inout):: kocn(0:)
  integer:: n
  ! First calculate iocn and mask from kocn
  ! and fix them as necessary.
  call calc_mask(mask=mask, kocn=kocn)
  call calc_iocn(iocn=iocn, mask=mask)
  !!>>>>>
  !!write(*,*) "calc_ocn before shaving: iocn=",iocn
  !!<<<<<
  call shave_convex(iocn=iocn, kocn=kocn, mask=mask, nshaved=n)
  !!>>>>>
  !!write(*,*) "calc_ocn after shaving: iocn=",iocn
  !!<<<<<
  if (n /= 0) then
    write(*,"(*(1X,G0))") "nshaved =", n, "convex steps have been shaved."
    call shave_convex(iocn=iocn, kocn=kocn, mask=mask, nshaved=n)
    if (n /= 0) then
      write(ERROR_UNIT,*) "***Second shaving: n = ", n
      error stop 7
    end if
  else
    write(*,*) "No shaving was necessary."
  end if
end subroutine calc_ocn


! Gridpoint numbers.
! num(i,k) > 0 for active points
! num(i,k) < 0 for inactive (unused) points.
!
! It's convinenent if the equations are defined
!   in the same order as nc = n(i,k).
! Compare with the order in which the equations are generated
!   in eqs_interior, eqs_surface_right, and eqs_slope_bottom.
subroutine calc_num(num, iocn, kocn) !, mask)
  use iso_fortran_env, only: ERROR_UNIT
  implicit NONE
  integer, intent(out):: num(0:, 0:)
  integer, intent(in) :: iocn(0:), kocn(0:)
  integer, allocatable:: flag(:,:) ! active gridpoints
  integer:: im, km, l, i, k, ntot
  im = size(num,1) - 2
  km = size(num,2) - 2

  allocate(flag(0:im+1, 0:km+1))
  flag = 0
  num(:,:) = -huge(1)
  l = 0

! Interior
  do k = 1, km
    do i = iocn(k), im
      if (i==iocn(k) .and. k==kocn(i)) cycle ! diagonal is not inerior
      l = l + 1
      num(i,k) = l
      flag(i,k) = 1
    end do
  end do

! Boundary Condition at the Surface (z = 0)
  do i = iocn(1), im
    l = l + 1
    num(i,0) = l
    flag(i,0) = 1
  end do

! Boundary condition offshore (x = L):
  do k = 1, kocn(im)
    l = l + 1
    num(im+1,k) = l
    flag(im+1,k) = 1
  end do

! Boundary condition along the slope-bottom.
! Under-the-ground gridpoints are active only for
!  "vertical" or "horizontal" sections.
  k = 1
  i = iocn(k)
  do while (i <= im)
    l = l + 1
    if (i == iocn(k) .and. k == kocn(i)) then ! next to bottom and wall
      num(i,k) = l
      flag(i,k) = 1
      write(*,"(1X,A,I4,I4,I6)") "calc_num: diag slope: i k l =", i, k, l
    else if (i == iocn(k)) then ! vertical
      num(i-1,k) = l
      flag(i-1,k) = 1
      write(*,"(1X,A,I4,I4,I6)") "calc_num: vert wall : i k l =", i, k, l
    else if (k == kocn(i)) then ! horizontal
      num(i,k+1) = l
      flag(i,k+1) = 1
      write(*,"(1X,A,I4,I4,I6)") "calc_num: flat bot  : i k l =", i, k, l
    else
      write(ERROR_UNIT,*)&
          "Something is wrong: Point must be adjacent to bndry."
      error stop 11
    end if

    ! move to next boundary point, following the bottom
    if (k == kocn(i)) then ! already at the bottom.
      i = i + 1 ! move right
      if (k < kocn(i)) then ! if not at the bottom,
        k = k + 1 ! move down
      end if
    else
      k = k + 1 ! move down
    end if
  end do

! Verify sizes
  ntot = sum(flag)
  if(l /= ntot) then
    write(ERROR_UNIT,*) 'num: Matrix Size Mismatch! l = ', l, ", N = ", ntot
    error stop 5
  else
    write(*,"(*(1X,G0))") "Number of eqs. l =", l
  end if

! Verify: num should be defined ( > 0) for all the active gridpoints.
  do k = 0, km+1
    do i = 0, im+1
      if ((flag(i,k) == 0 .and. num(i,k) >= 0) &
          .or. (flag(i,k) == 1 .and. num(i,k) <= 0) ) then
        write(ERROR_UNIT,*) "i, k, flag, num =", i, k, flag(i,k), num(i,k)
        error stop 2
      end if
    end do
  end do

end subroutine calc_num

end module define_matrices_diag

! =====================================================================
program main
  use ctwmodes_pars, only: im, km, f0&
      , double, single, outtype, freesurface&
      , grid_file, Ne_file, verbose
  use define_matrices_diag
  use, intrinsic:: ieee_arithmetic
  use, intrinsic:: iso_fortran_env, only: ERROR_UNIT
  implicit none

  !-- data files: set them to empty to use default data
  !character(*), parameter:: grid_file = "Grid.txt"
  !character(*), parameter:: Ne_file   = "Ne.txt"

  !-- Surface boundary condition
  real(double), parameter:: gravit_default = 9.8d0 !ignored if not freesurface

  ! If true, a lot of messages and diagnostics will be printed.
  !logical, parameter:: verbose = .true.

  integer:: i, k, l, num(0:im+1,0:km+1)
  integer:: M
  real(double):: dx(0:im+1)
  integer:: kocn(0:im+1), iocn(0:km+1)
  real(double):: N2e(0:km) ! N^2(z) defined at the edges of gridcells.
  real(double):: dz(0:km+1)
  integer:: mask(0:im+1, 0:km+1)

  real(double), allocatable:: A(:,:) ! lhs matrix
  real(double), allocatable:: B(:,:) ! rhs matrix
  real(double), allocatable:: alphaR(:) ! real part of eigenvalues
  real(double), allocatable:: alphaI(:) ! imag part of eigenvalues

  real(double), allocatable:: VR(:,:) ! right eigenvectors (output)
  real(double), allocatable:: tmpVec(:) ! for debug

  !-- Boundary type at x → ∞ (x = Lx)
  ! Use one of the three options below.
!  character(*), parameter:: rightboundary = "Px=0"
!  character(*), parameter:: rightboundary = "P=0"
!  character(*), parameter:: rightboundary = "wall"

  character(*), parameter:: ofile = "Eigenmodes-z.fort.bin"

  real(double):: gravit
  real(double):: tmp ! temporary variable
  logical:: Asym, Bsym, ltmp ! temporary variables
  integer:: itmp ! temporary variable
  integer:: uni

  write(*,*) "Ver 0.7.1: Eigenvalues will be c ."
!  write(*,*) "rightboundary = "//trim(rightboundary)
  write(*,*) "f0 = ", f0

! Grid and Topography
!   Set within subroutines or read from a file.
  call set_grid(dx, dz, kocn, infile=grid_file)

! Buoyancy Frequency
!   Set within subroutines or read from a file.
  call set_N2e(N2e, dz=dz, infile=Ne_file)

! calculate other geometry parameters.
  call calc_ocn(iocn, mask, kocn)
  call calc_num(num, iocn, kocn)

  if (verbose) then
    write(*,*) "num(0,   0)    = ", num(0,   0)
    write(*,*) "num(im+1,0)    = ", num(im+1,0)
    write(*,*) "num(im+1,km+1) = ", num(im+1,km+1)
    write(*,*) "num(0,   km+1) = ", num(0,   km+1)
  end if

  if (verbose) then
    do k = 0, km+1
      write(*,"(1X,A,I4,I12,F21.15)") "k iocn(k) dz(k) =", k, iocn(k), dz(k)
    end do
    do i = 0, im+1
      write(*,"(1X,A,I4,I12,F23.13)") "i kocn(i) dx(i) =", i, kocn(i), dx(i)
    end do
  end if

  M = maxval(num) ! max gridpoint number == the number of gridpoints.
  if (verbose) write(*,"(*(1X,G0))") "im, km, M, im*km - M =",&
      im, km, M, im*km - M

!=== Build the A and B matrices ===
! l is the equation number.
! For the gridpoint (i,k), nc is its gridpoint number and
!   nu, nd, nl, and nr are the gridpoint numbers of the upward, downward,
!   leftward, and rightward neighbors.
  if (freesurface) then
    gravit = gravit_default
    write(*,"(*(1X,G0))") "free surface: gravit =", gravit
  else
    gravit = ieee_value(1.d0, IEEE_QUIET_NAN)
    write(*,"(*(1X,G0))") "rigid lid"
  end if

  call eqs_all(A, B, l, iocn, kocn, num, dx, dz, N2e, f0, M, gravit=gravit)
  if (l == M) then
    write(*,"(*(1X,G0))") "l =", l, "equations have been defined."
  else
    write(ERROR_UNIT,"(*(1X,G0))") "eqs-gridpoints mismatch: l,M =", l, M
    error stop 2
  end if

  !--- Check whether the barotropic mode gives AF = 0
  ! Check whether BF[l] /= 0 only at the wall/slope.
  allocate(tmpVec(M))
  if (.not. freesurface) then
    write(*,*) "AF for barotropic mode:"
    tmp = 1.0d-14
    tmpVec(:) = 1
    tmpVec = matmul(A,tmpVec)
    ltmp = .true.
    do i = 1, size(tmpVec)
      if (abs(tmpVec(i)) > tmp) then
        write(*,*) "l, (AF)[l] = ", i, tmpVec(i)
        ltmp = .false.
      end if
    end do
    if (ltmp) write(*,*) "abs(AF[l]) <", tmp, "for all l"
  end if
  if (verbose) then
    write(*,*) "BF when F = (1,...,1):"
    tmpVec(:) = 1
    tmpVec = matmul(B,tmpVec)
    itmp = 0
    do i = 1, size(tmpVec)
      if (abs(tmpVec(i)) > 1.d-21) then
        write(*,"(2X,A,I8,1X,G0)") "B: l (FB)[l] =", i, tmpVec(i)
        itmp = itmp + 1
      end if
    end do
    if (itmp /= km) then
      write(ERROR_UNIT, *) &
          "There are", itmp, "nonzero BF values whereas km =", km
      error stop 77
    else
      write(*,"(1X,*(1X,G0))") &
          "There are exactly km = ", km, " nonzero vals in BF."
    end if
  end if
  deallocate(tmpVec)

!--- Solve the matrix eigenvalue problem.

  !-- Internal check: A and B should be symmetric ---
  Asym = is_symmetric(A)
  Bsym = is_symmetric(B)
  if (.not.(Asym .and. Bsym)) then
    write(ERROR_UNIT,*) "Asym, Bsym =", Asym, Bsym
    error stop 23
  end if

  if (freesurface) then
    write(*,*) "symmetric definite solver is used."
    call solve_sym(VR, alphaR, alphaI, A, B, f0)
  else
    write(*,*) "genenral solver is used."
    call solve_gen(VR, alphaR, alphaI, A, B, f0)
  end if

!---Output the results---
  open(newunit=uni, file=ofile, form="UNFORMATTED", access="SEQUENTIAL", &
      status="REPLACE", position="REWIND")
  write(uni) M, im, km
  write(uni) num
  write(uni) kocn
  write(uni) iocn
  write(uni) outtype
  write(uni) real(alphaR, kind=outtype)
  write(uni) real(alphaI, kind=outtype)
  write(uni) real(VR,     kind=outtype)
  write(uni) real(dx,     kind=outtype)
  write(uni) real(dz,     kind=outtype)
  write(uni) real(f0,     kind=outtype)
  write(uni) real(N2e,    kind=outtype)
  close(uni)

  write(*,"(1X,A)",advance="NO") "Output to "//trim(ofile)
  if (outtype == double) then
    write(*,"(A)") " in double precision"
  else if (outtype == single) then
    write(*,"(A)") " in single precision"
  else
    write(*,"(A)") " in unkown type"
  end if

  write(*,*) "Program finished successfully."
!  stop ! <- Will print an unnecessary message depending on the compiler.

contains
  ! General solver. Solves
  !    A x = (f/c) B x
  ! and returns c's in alphaR and corresponding eigenvectors in VR.
  ! The imaginary part alphaI should be all zero.
  !  (Otherwise, there must be some error somewhere. To signal such error
  !   is the only purpose of retaining alphaI.)
  subroutine solve_gen(VR, alphaR, alphaI, A, B, f0)
    implicit NONE
    real(double), intent(out), allocatable:: VR(:,:), alphaR(:), alphaI(:)
    real(double), intent(in):: A(:,:), B(:,:), f0

    character(1), parameter  :: JOBVL='N' ! Don't compute left eigenvectors
    character(1), parameter  :: JOBVR='V' ! Compute right eigenvectors.

    ! (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N
    ! are the eigenvalues. (output)
    real(double), allocatable:: BETA(:)
    real(double), allocatable:: VL(:,:) ! left eigenvectors (unused)
    integer                  :: LWORK   ! size of WORK(:) (input)
    real(double), allocatable:: WORK(:) ! work space
    integer                  :: INFO ! Error information if not 0 (output)

    integer:: M
    M = size(A,1)

    allocate( alphaR(M) )
    allocate( alphaI(M) )
    allocate( BETA(M) )
    !allocate( VL(M,M) ) ! not used.
    allocate( VR(M,M) )

    LWORK = 8 * M
    allocate( WORK(LWORK) )

    !-- Solve ---
    call DGGEV(JOBVL, JOBVR, M, A, M, B, M,&
        alphaR, alphaI, BETA, VL, M, VR, M, &
        WORK, LWORK, INFO)
    if (INFO /= 0) then; write(ERROR_UNIT,*) "INFO =", INFO; end if

    !-- Determine eigenvalues ---
    ! Input: alphaR(:), alphaI(:), beta(:).
    ! Output: alphaR = c, alphaI = 0.
    ! Explanation:
    !   alpha/beta is the eigenvalue = f/c.
    !   beta == 0 means that f/c is infinity.
    !   By the way that A and B are constructed, all eigenvalues should be real.
    !   So, we set
    !      alphaR = f0 * beta / alphaR = c
    !   and make sure that alphaI = 0.
    ! The special values INF and NAN are provided by the
    ! standard intrinsic module IEEE_ARITHMETIC .
    do l = 1, M
      if (alphaI(l) /= 0) then
        write(ERROR_UNIT,*) "l, alphaI(l) =", l, alphaI(l)
      end if
      if (alphaR(l) /= 0) then
        alphaR(l) = BETA(l)/alphaR(l) * f0 ! -> c
      else
        if (BETA(l) /= 0) then
          alphaR(l) = ieee_value(tmp, IEEE_POSITIVE_INF)
        else
          alphaR(l) = ieee_value(tmp, IEEE_QUIET_NAN)
        end if
      end if
    end do

    deallocate(BETA, WORK)
  end subroutine solve_gen

! Symmetric definite solver. Solves
!   (-B) x = (c/f) (-A) x
! because A is negative deifinte.
! B and A should be both symmetric.
! On exit, B contains the eigenvectors and A <- -A
  subroutine solve_sym(VR, cee, alphaI, A, B, f0)
    implicit NONE
    real(double), intent(out), allocatable:: VR(:,:), cee(:), alphaI(:)
    real(double), intent(inout):: A(:,:), B(:,:)
    real(double), intent(in):: f0

    integer, parameter:: ITYPE = 1 ! A x = lambda B x
    character(*), parameter::&
        JOBZ = 'V', & !Compute eigenvalues and eigenvectors.
        UPLO = 'U' ! Upper triangles of A and B are stored
    integer:: M, LDA, LDB, LWORK, INFO
    real(double), allocatable:: WORK(:)
    real(double):: tmpwork(1)

    M = size(A,1)

    allocate(cee(M))
    allocate(alphaI(M))
    allocate(VR(M,M))

    LDA = M ! leading dimension of A
    LDB = M ! leading dimension of B

    A = -A ! has to be positive definite
    B = -B

    !LWORK = (NB+2)*M ! optimal??

    ! Query dsygv to determine an optimal work-space size.
    ! Optimal value will be returned as WORK(1).
    ! It should be an integral value, but I use ceiling()
    ! just in case and to silence the compiler warning.
    LWORK = -1 ! query.
    call dsygv(ITYPE, JOBZ, UPLO, &
      M, B, LDB, A, LDA, &
      cee, tmpwork, LWORK, INFO )
    LWORK = ceiling(tmpwork(1)) ! Optimal value.
    if (verbose) write(*,*) "solve_sym: LWORK =", LWORK

    ! Actually solve
    allocate(WORK(LWORK))
    call dsygv(ITYPE, JOBZ, UPLO, &
        M, B, LDB, A, LDA, &
        cee, WORK, LWORK, INFO )
    if (INFO /= 0) then
      write(ERROR_UNIT,*) "INFO =", INFO
      error stop 15
    end if

    cee(:) = cee(:) * f0 ! eigenvalue = c/f
    alphaI = 0 ! dummy
    VR = B ! solution returned as B.

    deallocate(WORK)
  end subroutine solve_sym

!--- geometry

  ! Call read_grid if infile /= " ";
  ! or call set_dx, set_dz, and set_kocn.
  subroutine set_grid(dx,dz,kocn,infile)
    implicit NONE
    real(double), intent(out):: dx(0:), dz(0:)
    integer,      intent(out):: kocn(0:)
    character(*), intent(in) :: infile
    if (infile == " ") then
      write(*,*) "set grid internally."
      call set_dx(dx)
      call set_dz(dz)
      call set_kocn(kocn)
    else
      call read_grid(dx, dz, kocn, infile)
    end if
  end subroutine set_grid

  ! Set dx(:).
  ! This is just an example. Rewrite the subroutine as necessary.
  subroutine set_dx(dx)
    implicit NONE
    real(double), intent(out):: dx(0:)
    integer:: i, ii
    if (im < km) then
      write(ERROR_UNIT,*) "im < km: im, km =", im, km
      error stop 9
    end if
    ii = km - 2 - 1
    do i = 1, ii
      dx(i) = 1.5d+3
    end do
    do i = ii+1, im
      dx(i) = min( 1.1d0 * dx(i-1), 3.0d+3 )
    end do
    dx(0) = dx(1)
    dx(im+1) = dx(im)
  end subroutine set_dx

  ! Set dz(:).
  ! This is just an example. Rewrite the subroutine as necessary.
  subroutine set_dz(dz)
    implicit NONE
    real(double), intent(out):: dz(0:)
    integer:: k
    do k = 1, km
      dz(k) = 40.0d0
    end do
    dz(0) = dz(1)
    dz(km+1) = dz(km)
  end subroutine set_dz

  ! Set kocn(:)
  !  which is the bottom profile: the maximum k in the ocean at each i.
  subroutine set_kocn(kocn)
    implicit NONE
    integer, intent(out):: kocn(0:)
    integer:: i, kbot
    write(*,*) "Initially set kocn = undefined."
    kocn = - huge(1)
    if (im < km) error stop 3 ! There cannot be offshore region.
    do i = 1, km ! one-by-one step
      kbot = i ! one-by-one step
      kocn(i) = kbot
    end do
  end subroutine set_kocn

  ! Read grid and topography from file.
  ! *** The new grid file for ver.007 does NOT include lobes.
  !   dx(0), dx(im+1), dz(0), and dz(km+1) are therefore generated
  !   in this subroutine.
  subroutine read_grid(dx, dz, kocn, infile)
    implicit NONE
    real(double), intent(out):: dx(0:), dz(0:)
    integer,      intent(out):: kocn(0:)
    character(*), intent(in) :: infile
    !character(*), parameter:: infile = "Grid.txt"
    integer:: im_tmp, km_tmp, im, km, i_tmp, k_tmp
    !!character(255):: tmp_debug ! only for debugging

    !!>>> debug >>>
    !!call getcwd(tmp_debug)
    !!write(ERROR_UNIT, *) "pwd = ", trim(tmp_debug)
    !!<<< debug <<<

    if (size(kocn,1) /= size(dx,1)) then
      write(ERROR_UNIT,*) "size mismatch: size(kocn,1), size(dx,1) =",&
          size(kocn,1), size(dx,1)
      error stop 11
    end if
    im = size(dx,1) - 2
    km = size(dz,1) - 2

    kocn(:) = - huge(1) ! undef

    write(*,*) "reading grid from "//trim(infile)
    open(newunit=uni, file=infile, form="FORMATTED", &
        access="SEQUENTIAL", status="OLD", position="REWIND")

    ! Read header
    read(uni, *) im_tmp, km_tmp
    if (im_tmp /= im .or. km_tmp /= km) then
      write(ERROR_UNIT,*) "size mismatch: im, im_tmp =", im, im_tmp
      write(ERROR_UNIT,*) "               km, km_tmp =", km, km_tmp
      error stop 13
    end if

    ! Read dx(1:im), kocn(1:im)
    do i = 1, im
      read(uni, *) i_tmp, dx(i), kocn(i)
      if (i_tmp /= i) then
        write(ERROR_UNIT,*) "i_tmp, i =", i, i_tmp
        error stop 15
      end if
    end do

    ! Read dz(1:km)
    do k = 1, km
      read(uni, *) k_tmp, dz(k)
      if (k_tmp /= k) then
        write(ERROR_UNIT,*) "k_tmp, k =", k, k_tmp
        error stop 17
      end if
    end do

    close(uni)

    ! set dx and dz for the lobe regions.
    dx(0)    = dx(1)
    dx(im+1) = dx(im)
    dz(0)    = dz(1)
    dz(km+1) = dz(km)
  end subroutine read_grid

!--- N^2 ---

  ! Set N2(:), which is N^2.
  ! Call set_N2e if infile /= " ";
  ! or set N2e in here.
  subroutine set_N2e(N2e, dz, infile)
    implicit NONE
    real(double), intent(out):: N2e(:) !0:km
    real(double), intent(in) :: dz(0:) !0:km+1, cell center
    character(*), intent(in) :: infile
    real(double), parameter:: bvf = 0.0016d0 !c_1=ND/π =~ 2.8 m/s.
    if (infile == " ") then
      N2e(:) = bvf*bvf
      write(*,*) "const N^2 =", N2e(1)
    else
      call read_N2e(N2e, dz, infile)
    end if
  end subroutine set_N2e

  ! N^2 at cell edges k = 0, . . . , km
  ! get N(z) from a file
  subroutine read_N2e(N2e, dz, infile)
    implicit NONE
    !character(*), parameter:: infile = "Ne.txt"
    real(double), intent(out):: N2e(0:) !0:km, cell edges
    real(double), intent(in) :: dz(0:) !0:km+1, cell center
    character(*), intent(in) :: infile
    real(double), allocatable:: Ne(:)
    real(double):: zz, zz_file, tmp
    integer:: uni, kmp ! km+1
    open(newunit=uni,file=infile,action="READ",status="OLD", &
        form="FORMATTED", access="SEQUENTIAL", position="REWIND")
    read(uni, *) kmp
    if (kmp /= km+1) then
      write(ERROR_UNIT,"(*(1X,G0))") trim(infile) &
          ,"starts with", kmp, "but has to start with km+1 =", km+1
      error stop 5
    end if
    allocate(Ne(0:km))
    zz = 0
    do k = 0,km
      read(uni,*) zz_file, Ne(k)
      ! write(*,*) k, zz
      if (zz_file /= zz) then
        tmp = max(abs(zz_file), abs(zz))
        tmp = abs(zz_file - zz)/tmp ! -> relative error.
        if (tmp > 1.d-14) then
          write(ERROR_UNIT,*) "zz_file, zz =", zz_file, zz
          error stop 7
        end if
      end if
      if (k < km) then; zz = zz - dz(k+1); end if
    end do
    close(uni)
    N2e(:) = Ne(:)**2
    deallocate(Ne)
    write(*,*) "Ne(z) is read from "//trim(infile) &
        //" and N2e(z) is calculated."
  end subroutine read_N2e


!-- Utilities to check coding ---

! Check if the matrix is symmetric or not.
function is_symmetric(mat,eps) result(sym)
  implicit NONE
  real(double), intent(in):: mat(:,:)
  real(double), intent(in), optional::eps
  logical:: sym ! return value
  integer:: l,n
  real(double):: eps1, da
  if (present(eps)) then
    eps1 = eps
  else
    eps1 = 0 ! strict. To relax, use 1.d-14.
  end if
  if (size(mat,1) /= size(mat,2)) then
    write(ERROR_UNIT,*) "mat isn't square: size(mat) =", size(mat)
    error stop 13
  end if
  sym = .true.
  do n=2,size(mat,2)
    do l=1,n-1 ! upper triangle
      da = abs(mat(l,n) - mat(n,l))
      if (da /= 0) then
        da = da / max(abs(mat(l,n)), abs(mat(n,l)))
        if (da > eps1) then
          write(ERROR_UNIT,"(1x,A,*(1x,G0))") &
              "l, n, da, mat(l,n), mat(n,l) = ", l, n, da, mat(l,n), mat(n,l)
          sym = .false.
        else
          write(ERROR_UNIT,"(1x,A,*(1x,G0))") &
              "l, n, mat(l,n), da  = ", l, n, mat(l,n), da
        end if
      end if
    end do
  end do
  !if (sym) write(*,*) "mat is symmetric."
end function is_symmetric


end program main
