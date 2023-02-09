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

module xtb_prog_main
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stderr
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert
   use xtb_mctc_param
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment, init
   use xtb_prog_argparser
   use xtb_solv_state
   use xtb_setparam
   use xtb_sphereparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_fixparam
   use xtb_constrain_param, only : read_userdata
   use xtb_shake, only: init_shake
   use xtb_gfnff_shake, only: gff_init_shake => init_shake
   use xtb_embedding, only : init_pcem
   use xtb_io_reader, only : readMolecule
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_filetypes, only : fileType, getFileType, generateFileMetaInfo, &
      & generateFileName
   use xtb_readin
   use xtb_printout
   use xtb_setmod
   use xtb_propertyoutput
   use xtb_io_writer_turbomole, only : writeResultsTurbomole
   use xtb_io_writer_orca, only : writeResultsOrca
   use xtb_io_writer_gaussian, only : writeResultsGaussianExternal
   use xtb_restart
   use xtb_readparam
   use xtb_scc_core, only : iniqshell
   use xtb_aespot, only : get_radcn
   use xtb_iniq, only : iniqcn
   use xtb_eeq
   use xtb_disp_ncoord, only : ncoord_gfn, dncoord_erf, dncoord_d3, ncoord_erf, &
      & ncoord_d3
   use xtb_basis
   use xtb_axis, only : axis3
   use xtb_hessian, only : numhess
   use xtb_dynamic, only : md
   use xtb_modef, only : modefollow
   use xtb_mdoptim, only : mdopt
   use xtb_screening, only : screen
   use xtb_xtb_calculator
   use xtb_gfnff_calculator
   use xtb_paramset
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_main_setup
   use xtb_main_defaults, only : initDefaults
   use xtb_main_json, only : main_json, write_json_gfnff_lists
   use xtb_geoopt
   use xtb_metadynamic
   use xtb_biaspath
   use xtb_coffee
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_param, only : gff_print
   use xtb_gfnff_topology, only : TPrintTopo
   use xtb_gfnff_convert, only : struc_convert
   use xtb_scan
   use xtb_kopt
   use xtb_oniom, only : oniom_input
   implicit none
   private

   public :: xtbMain


contains


subroutine xtbMain(env, argParser)

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_main"

   type(TEnvironment), intent(inout) :: env

   type(TArgParser), intent(inout) :: argParser

!! ========================================================================
!  use some wrapper types to bundle information together
   type(TMolecule) :: mol
   type(scc_results) :: res
   class(TCalculator), allocatable :: calc
   type(freq_results) :: fres
   type(TRestart) :: chk
   type(chrg_parameter) :: chrgeq
   type(oniom_input) :: oniom
!  store important names and stuff like that in FORTRAN strings
   character(len=:),allocatable :: fname    ! geometry input file
   character(len=:),allocatable :: xcontrol ! instruction file
   character(len=:),allocatable :: xrc      ! global instruction file
   character(len=:),allocatable :: fnv      ! parameter file
   character(len=:),allocatable :: tmpname  ! temporary string
   character(len=:),allocatable :: cdum     ! temporary string
   character(len=:),allocatable :: extension, basename, directory
   integer :: ftype

!! ========================================================================
!  default names for important files in xtb
   character(len=*),parameter :: p_fname_rc = '.xtbrc'
   character(len=*),parameter :: p_fname_param_gfn0  = 'param_gfn0-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfn1  = 'param_gfn1-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfn2  = 'param_gfn2-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfnff = '.param_gfnff.xtb'
   character(len=*),parameter :: p_fname_param_ipea  = 'param_ipea-xtb.txt'

   integer :: gsolvstate
   integer :: i,j,k,l,idum,xy
   integer :: ich,ictrl,iprop ! file handle
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: cn  (:)
   real(wp),allocatable :: sat (:)
   real(wp),allocatable :: g   (:,:)
   real(wp) :: vec3(3)
   type(TxTBParameter) :: globpar
   real(wp),allocatable :: dcn (:,:,:)
   real(wp),allocatable :: dq  (:,:,:)
   real(wp),allocatable :: dumdumdum  (:,:,:)
   real(wp),allocatable :: q  (:)
   real(wp),allocatable :: ql  (:)
   real(wp),allocatable :: qr  (:)
   character (len=2550) :: cwd
!! ------------------------------------------------------------------------
   integer,external :: ncore

!! ------------------------------------------------------------------------
   logical :: struc_conversion_done = .false.

!! ========================================================================
!  debugging variables for numerical gradient
   logical, parameter    :: gen_param = .false.
   logical, parameter    :: debug = .false.
   type(TRestart) :: wf0
   real(wp),allocatable  :: coord(:,:),numg(:,:),gdum(:,:)
   real(wp) :: sdum(3,3)
   real(wp),parameter    :: step = 0.00001_wp, step2 = 0.5_wp/step
   real(wp) :: er,el
   logical  :: coffee ! if debugging gets really though, get a coffee

!! ------------------------------------------------------------------------
!  undocumented and unexplainable variables go here
   integer  :: nFiles, iFile
   integer  :: rohf,err
   real(wp) :: dum5,egap,etot,ipeashift
   real(wp) :: zero,t0,t1,w0,w1,acc,etot2,g298
   real(wp) :: one,two
   real(wp) :: ea,ip
   real(wp) :: vomega,vfukui
   real(wp),allocatable :: f_plus(:), f_minus(:)
   type(TRestart) :: wf_p, wf_m
   parameter (zero=0.0_wp)
   parameter (one =1.0_wp)
   parameter (two =2.0_wp)
   logical :: ex,okbas
   logical :: epr,diff,murks
   logical :: exist
   logical :: lgrad,restart
   logical :: copycontrol
   logical :: newreader
   logical :: strict
   logical :: exitRun
   logical :: cold_fusion

!  OMP stuff
   integer :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
   integer :: nproc

   type(TPrintTopo) :: printTopo ! gfnff topology printout list

   xenv%home = env%xtbhome
   xenv%path = env%xtbpath


   ! ------------------------------------------------------------------------
   !> read the command line arguments
   call parseArguments(env, argParser, xcontrol, fnv, acc, lgrad, &
      & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom)

   nFiles = argParser%countFiles()
   select case(nFiles)
   case(0)
      if (.not.coffee) then
         if(printTopo%warning) call env%error("Eventually the input file was given to wrtopo as an argument.",source)
         call env%error("No input file given, so there is nothing to do", source)
      else
         fname = 'coffee'
      end if
   case(1:)
      do iFile = 1, nFiles-1
         call argParser%nextFile(fname)
         call env%warning("Input file '"//fname//"' will be ignored", source)
      end do
      call argParser%nextFile(fname)
   end select

   if (.not.allocated(xcontrol)) then
      if (copycontrol) then
         xcontrol = 'xtb.inp'
      else
         xcontrol = fname
      end if
   end if

   call env%checkpoint("Command line argument parsing failed")


   ! ------------------------------------------------------------------------
   !> read the detailed input file
   call rdcontrol(xcontrol, env, copy_file=copycontrol)

   call env%checkpoint("Reading '"//xcontrol//"' failed")


   ! ------------------------------------------------------------------------
   !> read dot-Files before reading the rc and after reading the xcontrol
   !> Total molecular charge
   call open_file(ich,'.CHRG','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.CHRG is empty!', source)
      else
         call set_chrg(env,cdum)
         call close_file(ich)
      end if
   end if

   call env%checkpoint("Reading charge from file failed")

   !> Number of unpaired electrons
   call open_file(ich,'.UHF','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.UHF is empty!', source)
      else
         call set_spin(env,cdum)
         call close_file(ich)
      end if
   endif
   
   !> efield read: gfnff only
   call open_file(ich,'.EFIELD','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.EFIELD is empty!', source)
      else
         call set_efield(env,cdum)
         call close_file(ich)
      end if
   endif

   call env%checkpoint("Reading multiplicity from file failed")


   ! ------------------------------------------------------------------------
   !> read the xtbrc if you can find it (use rdpath directly instead of xfind)
   call rdpath(env%xtbpath, p_fname_rc, xrc, exist)
   if (exist) then
      call rdcontrol(xrc, env, copy_file=.false.)

      call env%checkpoint("Reading '"//xrc//"' failed")
   endif


   ! ------------------------------------------------------------------------
   !> FIXME: some settings that are still not automatic
   !> Make sure GFN0-xTB uses the correct exttyp
   if(set%gfn_method == 0)  call set_exttyp('eht')
   rohf = 1 ! HS default
   egap = 0.0_wp
   ipeashift = 0.0_wp


   ! ========================================================================
   !> no user interaction up to now, time to show off!
   !> print the xtb banner with version number and compilation date
   !> making a fancy version of this is hard, x is difficult in ASCII art
   !> print current time


   ! ------------------------------------------------------------------------
   !> get molecular structure
   if (coffee) then ! it's coffee time
      fname = 'caffeine'
      call get_coffee(mol)
      call generateFileMetaInfo(fname, directory, basename, extension)
   else                                                              
      call generateFileMetaInfo(fname, directory, basename, extension)
      ftype = getFileType(basename, extension)
      call open_file(ich, fname, 'r')
      call readMolecule(env, mol, ich, ftype)
      
        

!   call getCWD(cwd)
!   open (12, file = trim(cwd)//'/nbf.txt', status = 'new', action="write")
!   write(12,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!nbf'')')
!   write(12,*) mol%n
!   write(12,*) mol%sym(1)
!   do xy = 1, size(mol%nbf, dim=1)
!   write(12,"(*(g0,1X))") mol%nbf(xy,:) ! 一次打印一行
!   end do


      call close_file(ich)
      if (mol%struc%two_dimensional) then
         call env%warning("Two dimensional input structure detected", source)
      end if

      ! Special CT input file case
      if (mol%chrg /= 0.0_wp) then
         if (set%ichrg == 0) then
            set%ichrg = nint(mol%chrg)
         else
            call env%warning("Charge in sdf/mol input was overwritten", source)
         end if
      end if

      call env%checkpoint("reading geometry input '"//fname//"' failed")
   end if

   ! ------------------------------------------------------------------------
   !> get some memory
   allocate(cn(mol%n),sat(mol%n),g(3,mol%n), source = 0.0_wp)
   atmass = atomic_mass(mol%at) * autoamu ! from splitparam.f90
   set%periodic = mol%npbc > 0

   do i=1,mol%n
      mol%z(i) = mol%at(i) - ncore( mol%at(i) )
      ! lanthanides without f are treated as La
      if(mol%at(i).gt.57.and.mol%at(i).lt.72) mol%z(i)=3
   enddo

   !> initialize time step for MD if requested autocomplete
   if (set%tstep_md < 0.0_wp) then
      set%tstep_md = (minval(atmass)/(atomic_mass(1)*autoamu))**(1.0_wp/3.0_wp)
   endif

   mol%chrg = real(set%ichrg, wp)
   mol%uhf = set%nalphabeta
   call initrand





   ! ------------------------------------------------------------------------
   !> write copy of detailed input
   if (copycontrol) then
      call open_set(ictrl,xcontrol)
      call write_set(ictrl)
      call close_set(ictrl)
   endif

   ! ------------------------------------------------------------------------
   !> if you have requested a define we stop here...
   if (set%define) then
      if (set%verbose) call main_geometry(env%unit,mol)
      call eval_define(set%veryverbose)
   endif
   call env%show('Please study the warnings concerning your input carefully')
   call raise('F', 'Please study the warnings concerning your input carefully')

   ! ========================================================================
   !> From here we switch to the method setup
   !> enable error on warnings
   if (strict) call mctc_strict
   env%strict = strict



   ! ------------------------------------------------------------------------
   !> Obtain the parameter data
   call newCalculator(env, mol, calc, fnv, restart, acc, oniom)
   call env%checkpoint("Could not setup single-point calculator")

   call initDefaults(env, calc, mol, gsolvstate)
   call env%checkpoint("Could not setup defaults")






   ! ========================================================================
   !> PRINTOUT SECTION
   if (allocated(set%property_file)) then
      call open_file(iprop,set%property_file,'w')
      if (iprop.eq.-1) then
         iprop = env%unit
         deallocate(set%property_file)
      else
         write(env%unit,'(/,a)') "Property printout bound to '"//set%property_file//"'"
         if (allocated(cdum)) deallocate(cdum)
         call get_command(length=l)
         allocate( character(len=l) :: cdum )
         call get_command(cdum)
         write(iprop,'("command:  ''",a,"''")') cdum
         call rdvar('HOSTNAME',cdum,err)
         if (err.eq.0) &
            write(iprop,'("hostname: ''",a,"''")') cdum
         write(iprop,'("date:     ",a)') prtimestring('S')
      endif
   else
      iprop = env%unit
   endif



   if(printTopo%any()) then
     select type(calc)
       type is(TGFFCalculator)
         call write_json_gfnff_lists(mol%n,calc%topo,chk%nlist,printTopo)
     end select
   endif



   ! ------------------------------------------------------------------------
   !  make some post processing afterward, show some timings and stuff

   call terminate(0)


end subroutine xtbMain


!> Parse command line arguments and forward them to settings
subroutine parseArguments(env, args, inputFile, paramFile, accuracy, lgrad, &
      & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom)
   use xtb_mctc_global, only : persistentEnv

   !> Name of error producer
   character(len=*), parameter :: source = "prog_main_parseArguments"

   !> Calculation environment
   type(TEnvironment) :: env

   !> Command line argument parser
   type(TArgParser) :: args

   !> Detailed input file name
   character(len=:),allocatable,intent(out) :: inputFile

   !> Parameter file name
   character(len=:),allocatable,intent(out) :: paramFile

   !> Accuracy number for numerical thresholds
   real(wp), intent(out) :: accuracy

   !> Reference state for solvation free energies
   integer, intent(out) :: gsolvstate

   !> Restart calculation
   logical, intent(out) :: restart

   !> Handle warnings as errors
   logical, intent(out) :: strict

   !> Debugging with a lot of caffeine
   logical, intent(out) :: coffee

   !> topology printout list
   type(TPrintTopo), intent(out) :: printTopo

   !> Print the gradient to file
   logical, intent(out) :: lgrad

   !> Copy the detailed input file
   logical, intent(out) :: copycontrol

   !> Input for ONIOM model
   type(oniom_input), intent(out) :: oniom

!$ integer :: omp_get_num_threads, nproc
   integer :: nFlags
   integer :: idum, ndum
   real(wp) :: ddum
   character(len=:), allocatable :: flag, sec
   logical :: exist

   set%gfn_method = 2
   coffee = .false.
   strict = .false.
   restart = .true.
   copycontrol = .false.
   lgrad = .false.
   accuracy = 1.0_wp
   gsolvstate = solutionState%gsolv

   nFlags = args%countFlags()
   call args%nextFlag(flag)
   do while(allocated(flag))
      if (len(flag) > 2 .and. flag(1:1) == '-' .and. flag(1:2) /= '--') then
         call env%warning("the use of '"//flag//"' is discouraged, "// &
            & "please use '-"//flag//"' next time", source)
         flag = '-'//flag
      end if
      select case(flag)
      case default
         call env%warning("Unknown option '"//flag//"' provided", source)

      case('-h', '--help')
         call help(env%unit)
         call terminate(0)

      case('--citation')
         call citation(env%unit)
         call terminate(0)

      case('--license')
         call disclamer(env%unit)
         call terminate(0)

      case('--version')
         call xtb_header(env%unit)
         call terminate(0)

      case('-v','--verbose')
         set%verbose = .true.

      case('-V','--very-verbose')
         set%verbose = .true.
         set%veryverbose = .true.

      case(     '--define')
         call set_define

      case('-P','--parallel')
   !$    if (.false.) then
            call env%warning('Program compiled without threading support', source)
   !$    endif
         ! Always remove next argument to keep argument parsing consistent
         call args%nextArg(sec)
   !$    if (allocated(sec)) then
   !$    if (getValue(env,sec,idum)) then
   !$       nproc = omp_get_num_threads()
   !$       call omp_set_num_threads(idum)
#ifdef WITH_MKL
   !$       call mkl_set_num_threads(idum)
#endif
   !$    endif
   !$    endif

      case('--gfnff')
         call set_exttyp('ff')
      
      case('--gff')
         call set_exttyp('ff')

      case('--wrtopo')
         call args%nextArg(sec)
         if (allocated(sec)) then
           call setWRtopo(sec,printTopo)
           if(printTopo%warning) call env%error("A wrtopo argument has been misspelled.",source)
         else
           call env%error("The wrtopo keyword is missing an argument.",source)
         endif
      end select
      call args%nextFlag(flag)
   end do

end subroutine parseArguments

function read_whole_file(fname) result(list)
   character(len=*), intent(in) :: fname
   character(len=:), allocatable :: list
   integer :: io, stat
   character(len=:), allocatable :: line
   open(newunit=io, file=fname, iostat=stat)
   call getline(io, list, stat)
   do while(stat == 0)
      call getline(io, line, stat)
      if (stat == 0) list = list // "," // line
   end do
   close(io, iostat=stat)
end function read_whole_file

! set booleans for requested topology list printout
subroutine setWRtopo(sec,printTopo)
   ! command line argument
   character(len=*), intent(in) :: sec
   ! type holds booleans of to be printed topology lists
   type(TPrintTopo), intent(inout) :: printTopo
   ! seperator for lists is ","
   character, parameter :: sep = ","
   ! current and old position of seperator
   integer :: curr_pos, old_pos
   integer :: lenSec, i

   curr_pos = 0
   old_pos = 0
   lenSec = len(sec)
   do i=1, lenSec
     curr_pos = scan(sec(curr_pos+1:lenSec),sep)+old_pos
     if(curr_pos.ne.old_pos) then
       call selectList(sec(old_pos+1:curr_pos-1),printTopo)
     else
       call selectList(sec(old_pos+1:lenSec),printTopo)
       exit
     endif
     old_pos=curr_pos
   enddo

end subroutine setWRtopo

subroutine selectList(secSplit, printTopo)
   ! part of command line argument
   character(len=*), intent(in) :: secSplit
   ! holds booleans of to be printed topology lists
   type(TPrintTopo), intent(inout) :: printTopo

   select case(secSplit)
   case("nb")
     printTopo%nb = .true.
   case("bpair")
     printTopo%bpair = .true.
   case("alist")
     printTopo%alist = .true.
   case("blist")
     printTopo%blist = .true.
   case("tlist")
     printTopo%tlist = .true.
   case("vtors")
     printTopo%vtors = .true.
   case("vbond")
     printTopo%vbond = .true.
   case("vangl")
     printTopo%vangl = .true.
   case("hbbond")
      printTopo%hbbond = .true.
   case("eeq")
      printTopo%eeq = .true.
   case default
     printTopo%warning = .true.
   end select
end subroutine selectList

end module xtb_prog_main
