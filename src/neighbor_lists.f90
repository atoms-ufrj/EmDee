!   This file is part of EmDee.
!
!    EmDee is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    EmDee is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with EmDee. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu@eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module neighbor_lists

use EmDeeData

implicit none

integer, parameter, private :: nb(3,nbcells) = reshape( [ &
 0,  0,  1,    0,  1,  0,    1,  0,  0,   -1,  0,  1,   -1,  1,  0,    0, -1,  1,    0,  1,  1, &
 1,  0,  1,    1,  1,  0,   -1, -1,  1,   -1,  1,  1,    1, -1,  1,    1,  1,  1,    0,  0,  2, &
 0,  2,  0,    2,  0,  0,   -2,  0,  1,   -2,  1,  0,   -1,  0,  2,   -1,  2,  0,    0, -2,  1, &
 0, -1,  2,    0,  1,  2,    0,  2,  1,    1,  0,  2,    1,  2,  0,    2,  0,  1,    2,  1,  0, &
-2, -1,  1,   -2,  1,  1,   -1, -2,  1,   -1, -1,  2,   -1,  1,  2,   -1,  2,  1,    1, -2,  1, &
 1, -1,  2,    1,  1,  2,    1,  2,  1,    2, -1,  1,    2,  1,  1,   -2,  0,  2,   -2,  2,  0, &
 0, -2,  2,    0,  2,  2,    2,  0,  2,    2,  2,  0,   -2, -2,  1,   -2, -1,  2,   -2,  1,  2, &
-2,  2,  1,   -1, -2,  2,   -1,  2,  2,    1, -2,  2,    1,  2,  2,    2, -2,  1,    2, -1,  2, &
 2,  1,  2,    2,  2,  1,   -2, -2,  2,   -2,  2,  2,    2, -2,  2,    2,  2,  2], shape(nb) )

contains

!===================================================================================================

  real(rb) function maximum_approach_sq( N, delta )
    integer,  intent(in) :: N
    real(rb), intent(in) :: delta(3,N)

    integer  :: i
    real(rb) :: maximum, next, deltaSq

    maximum = sum(delta(:,1)**2)
    next = maximum
    do i = 2, N
      deltaSq = sum(delta(:,i)**2)
      if (deltaSq > maximum) then
        next = maximum
        maximum = deltaSq
      end if
    end do
    maximum_approach_sq = maximum + 2*sqrt(maximum*next) + next

  end function maximum_approach_sq

!===================================================================================================

  subroutine distribute_atoms( me, M, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: M
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer :: MM, cells_per_thread
    integer :: maxNatoms(me%nthreads), threadNatoms(me%nthreads), next(me%natoms)
    logical :: make_cells

    MM = M*M
    make_cells = M /= me%mcells
    if (make_cells) then
      me%mcells = M
      me%ncells = M*MM
      if (me%ncells > me%maxcells) then
        deallocate( me%cell, me%cellAtom%first, me%cellAtom%last, me%atomsInCell )
        allocate( me%cell(me%ncells), me%cellAtom%first(me%ncells), me%cellAtom%last(me%ncells) )
        allocate( me%atomsInCell(me%ncells) )
        call me % threadCell % allocate( 0, me%nthreads )
        me%maxcells = me%ncells
      end if
      cells_per_thread = (me%ncells + me%nthreads - 1)/me%nthreads
    end if

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call distribute( thread, maxNatoms(thread) )
    end block
    !$omp end parallel
    me%maxatoms = maxval(maxNatoms)
    me%maxpairs = (me%maxatoms*((2*nbcells + 1)*me%maxatoms - 1))/2

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine distribute( thread, maxNatoms )
        integer, intent(in)  :: thread
        integer, intent(out) :: maxNatoms

        integer :: i, j, k, icell, ix, iy, iz, first, last
        integer :: icoord(3)
        integer, allocatable :: head(:)

        if (make_cells) then
          first = (thread - 1)*cells_per_thread + 1
          last = min( thread*cells_per_thread, me%ncells )
          do icell = first, last
            k = icell - 1
            iz = k/MM
            j = k - iz*MM
            iy = j/M
            ix = j - iy*M
            me%cell(icell)%neighbor = 1 + pbc(ix+nb(1,:)) + pbc(iy+nb(2,:))*M + pbc(iz+nb(3,:))*MM
          end do
          me%threadCell%first(thread) = first
          me%threadCell%last(thread) = last
        else
          first = me%threadCell%first(thread)
          last = me%threadCell%last(thread)
        end if

        do i = (thread - 1)*me%threadAtoms + 1, min( thread*me%threadAtoms, me%natoms )
          icoord = int(M*(Rs(:,i) - floor(Rs(:,i))),ib)
          me%atomCell(i) = 1 + icoord(1) + M*icoord(2) + MM*icoord(3)
        end do
        !$omp barrier

        allocate( head(first:last) )
        head = 0
        me%atomsInCell(first:last) = 0
        do i = 1, me%natoms
          icell = me%atomCell(i)
          if ((icell >= first).and.(icell <= last)) then
            next(i) = head(icell)
            head(icell) = i
            me%atomsInCell(icell) = me%atomsInCell(icell) + 1
          end if
        end do
        threadNatoms(thread) = sum(me%atomsInCell(first:last))
        !$omp barrier

        maxNatoms = 0
        k = sum(threadNatoms(1:thread-1))
        do icell = first, last
          me%cellAtom%first(icell) = k + 1
          i = head(icell)
          do while (i /= 0)
            k = k + 1
            me%cellAtom%item(k) = i
            i = next(i)
          end do
          me%cellAtom%last(icell) = k
          if (me%atomsInCell(icell) > maxNatoms) maxNatoms = me%atomsInCell(icell)
        end do
      end subroutine distribute
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental integer function pbc( x )
        integer, intent(in) :: x
        if (x < 0) then
          pbc = x + M
        else if (x >= M) then
          pbc = x - M
        else
          pbc = x
        end if
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine distribute_atoms

!===================================================================================================

  subroutine handle_neighbor_lists( me, builds, time, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(inout) :: builds
    real(rb),    intent(inout) :: time
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer :: M

    time = time - omp_get_wtime()
    if (maximum_approach_sq( me%natoms, me%R - me%R0 ) > me%skinSq) then
      M = floor(ndiv*me%Lbox/me%xRc)
      call distribute_atoms( me, max(M,2*ndiv+1), Rs )
      me%R0 = me%R
      builds = builds + 1
      !$omp parallel num_threads(me%nthreads)
      call build_neighbor_lists( me, omp_get_thread_num() + 1, Rs )
      !$omp end parallel
    endif
    time = time + omp_get_wtime()

  end subroutine handle_neighbor_lists

!===================================================================================================

  subroutine build_neighbor_lists( me, thread, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer  :: i, j, k, m, icell, npairs, itype, ibody, ipairs, middle
    integer  :: nlocal, ntotal, first
    real(rb) :: xRc2, xInRc2, r2
    logical  :: include(me%natoms)
    integer  :: atom((nbcells+1)*me%maxatoms)

    integer,  allocatable :: xlist(:)
    real(rb), allocatable :: Ratom(:,:), rsq(:)

    xRc2 = me%xRcSq*me%invL2
    xInRc2 = me%xExRcSq*me%invL2

    include = .true.
    npairs = 0
    associate ( neighbor => me%neighbor(thread), c1 => me%cellAtom%first, cN => me%cellAtom%last )
      do icell = me%threadCell%first(thread), me%threadCell%last(thread)

        associate ( neigh => me%cell(icell)%neighbor )
          nlocal = me%atomsInCell(icell)
          ntotal = nlocal + sum(me%atomsInCell(neigh))
          atom(1:nlocal) = me%cellAtom%item(c1(icell):cN(icell))
          atom(nlocal+1:ntotal) = [(me%cellAtom%item(c1(neigh(m)):cN(neigh(m))),m=1,nbcells)]
        end associate

        if (neighbor%nitems < npairs + nlocal*ntotal) then
          call neighbor % resize( npairs + nlocal*ntotal + extra )
        end if

        Ratom = Rs(:,atom(1:ntotal))
        do k = 1, nlocal
          i = atom(k)
          first = npairs + 1
          neighbor%first(i) = first
          associate ( jlist => atom(k+1:ntotal),      &
                      item => neighbor%item(first:),  &
                      value => neighbor%value(first:) )
            ipairs = 0
            middle = 0
            itype = me%atomType(i)
            ibody = me%atomBody(i)
            xlist = me%excluded%item(me%excluded%first(i):me%excluded%last(i))
            include(xlist) = .false.
            rsq = [(sum(pbc(Ratom(:,k) - Ratom(:,m))**2),m=k+1,ntotal)]
            do m = 1, size(jlist)
              r2 = rsq(m)
              if (r2 < xRc2) then
                j = jlist(m)
                if (include(j)) then
                  if (me%atomBody(j) /= ibody) then
                    if (me%interact(me%atomType(j),itype)) then
                      call insert_neighbor( item, value, ipairs, j, r2 )
                      if (r2 < xInRc2) middle = middle + 1
                    end if
                  end if
                end if
              end if
            end do
          end associate
          include(xlist) = .true.

          neighbor%middle(i) = npairs + middle
          npairs = npairs + ipairs
          neighbor%last(i) = npairs

        end do

      end do
      neighbor%count = npairs
    end associate

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental function pbc( x )
        real(rb), intent(in) :: x
        real(rb)              :: pbc
        pbc = x - anint(x)
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine insert_neighbor( item, value, npairs, j, r2 )
        integer,     intent(inout) :: item(:)
        real(rb),    intent(inout) :: value(:)
        integer,     intent(inout) :: npairs
        integer,     intent(in)    :: j
        real(rb),    intent(in)    :: r2
        integer :: k
        do k = npairs, 1, -1
          if (value(k) < r2) exit
          value(k+1) = value(k)
          item(k+1) = item(k)
        end do
        value(k+1) = r2
        item(k+1) = j
        npairs = npairs + 1
      end subroutine insert_neighbor
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine build_neighbor_lists

!===================================================================================================

end module neighbor_lists
