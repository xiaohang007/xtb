! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module xtb_io_reader_top
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_strings, only : parse
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_type_molecule
   use xtb_pbc_tools
   use xtb_readin, only : getline => strip_line
   implicit none
   private

   public :: readMoleculeTopo


   logical, parameter :: debug = .false.


contains


subroutine readMoleculeTopo(mol, unit, status, iomsg)
   class(TMolecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer  :: i,ii, n,nn, iat
   character(len=symbolLength), allocatable :: sym(:)
   integer,allocatable ::     nbf(:,:)
   real(wp) :: x, y, z
   real(wp) :: conv

   character(len=:),allocatable :: message
   character(len=:),allocatable :: line
   character(len=symbolLength) :: chdum
   character(len=symbolLength) :: args(1000000)   ! 原子数上限为256，估计后面要改高一些
   integer  :: err
   integer  :: xy

   status = .false.

   conv = aatoau

   read(unit,*,iostat=err) n
   if (err.ne.0) then
      iomsg = "Could not read number of atoms, check format!"
      return
   endif

   if (n.lt.1) then
      iomsg = "Found no atoms, cannot work without atoms!"
      return
   endif

   allocate(sym(n))
   allocate(nbf(20, n))

   ! read sym
   call getline(unit,line,err)
   call parse(line,' ',args,nn)
   if (nn.ne.n) then
   iomsg = "Wrong number of atomic symbols"
   return
   endif
   do i = 1, nn
       sym(i) = trim(args(i))
   enddo

   do ii = 1, 20
      call getline(unit,line,err)
      if (nn.ne.n) then
      iomsg = "Wrong number of atomic symbols"
      return
      endif
      read(line,*,iostat=err) nbf(ii,:)
      if (err.ne.0) then
         iomsg = "Could not parse coordinates from Xmol file"
         return
      endif
   enddo

   call init2(mol, sym, nbf)
   status = .true.

end subroutine readMoleculeTopo


end module xtb_io_reader_top
