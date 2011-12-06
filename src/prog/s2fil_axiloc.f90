!------------------------------------------------------------------------------
! s2fil_axiloc
!
!! Detect azimuthally symmetric (i.e. axisymmetric) objects of
!! variying size in a sky.
!!
!! Usage: s2fil_axiloc
!!   - [-help]: Display usage information.
!!   - [-inp filename_inp]: Name of sky file (that we attempt to detect 
!!     objects in).
!!   - [-filter_data filename_filter_data]: Name of filter data file.
!!   - [-nside nside]: Healpix resolution to consider.
!!   - [-lmax lmax]: Maximum harmonic l to consider.
!!   - [-out_prefix filename_out_prefix]: Prefix of output files.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   December 2012 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2fil_axiloc

  use s2_types_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_inp
  character(len=S2_STRING_LEN) :: filename_filter_data
  character(len=S2_STRING_LEN) :: filename_out_prefix

  character(len=1), parameter :: COMMENT_CHAR = '#'
  integer :: fileid = 21

  integer :: nside, lmax
  type(s2_sky) :: sky

  ! Read command line options.
  call parse_options()


  !----------------------------------------------------------------------------
  ! Read data
  !----------------------------------------------------------------------------

  ! Read input sky.
  sky = s2_sky_init(filename_inp, S2_SKY_FILE_TYPE_MAP)

  ! Read filter data file.



  !----------------------------------------------------------------------------
  ! Apply matched filters
  !----------------------------------------------------------------------------





  !----------------------------------------------------------------------------
  ! Threshold filtered maps
  !----------------------------------------------------------------------------



  ! ...


  !----------------------------------------------------------------------------
  ! Free memory
  !----------------------------------------------------------------------------

  call s2_sky_free(sky)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - February 2005
    !
    ! Revisions:
    !   February 2005 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: n, i
      character(len=S2_STRING_LEN) :: opt
      character(len=S2_STRING_LEN) :: arg
      
      n = nArguments()
     
      do i = 1,n,2
        
        call getArgument(i,opt)
     
        if (i == n .and. trim(opt) /= '-help') then
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a,a)') 'Usage: s2fil_axiloc ', &
              '[-inp filename_inp]'
            write(*,'(a,a)') '                    ', &
              '[-filter_data filename_filter_data]'
            write(*,'(a,a)') '                    ', &
              '[-nside nside]'
            write(*,'(a,a)') '                    ', &
              '[-lmax lmax]'
            write(*,'(a,a)') '                    ', &
              '[-out out_prefix]'
            stop

          case ('-inp')
            filename_inp = trim(arg) 

          case ('-filter_data')
            filename_filter_data = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax

         case ('-out')
            filename_out_prefix = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2fil_axiloc
