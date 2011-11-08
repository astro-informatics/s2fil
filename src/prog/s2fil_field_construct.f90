!------------------------------------------------------------------------------
! s2fil_field_construct
!
!! Construct filtered field from optimal filter and data sky.
!!
!! Usage: s2fil_field_construct
!!   - [-help]: Display usage information.
!!   - [-filter filename_filter]: Name of optimal filter to read and to be 
!!     used to filter the data sky map.  Note that filter file header also
!!     contains names of associated child files that also contain filter 
!!     data to be read.  
!!   - [-data filename_sky]: Name of data sky fits file to read sky data from.
!!   - [-filetype sky_file_type_str]: String specifying type of data sky
!!     fits file, either 'map' or 'sky'.
!!   - [-field filename_field]: Name of field file to write.  Note field 
!!     data is also written to assocaited child files, the names of which
!!      are stored in the field file header.
!!   - [-ngamma n_gamma]: Number of orientations to consider when filtering.
!!   - [-write_filter write_filter]: Logical to specify whether to write 
!!     an optimal filter file as an associated child file of the s2fil field.
!!     This may not be necessary since the optimal filter will most likely be
!!     read from a file of this type and so to avoid duplication another
!!     need not be written.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2fil_field_construct

  use s2_types_mod
  use s2_sky_mod
  use s2fil_error_mod
  use s2fil_field_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_filter
  character(len=S2_STRING_LEN) :: filename_sky
  character(len=S2_STRING_LEN) :: filename_field
  character(len=S2_STRING_LEN) :: sky_file_type_str
  character(len=S2_STRING_LEN), parameter :: &
    FILE_TYPE_MAP = 'map', &
    FILE_TYPE_SKY = 'sky'
  integer :: sky_file_type = S2_SKY_FILE_TYPE_MAP
  integer :: n_gamma = 1
  type(s2fil_field) :: field
  logical :: write_filter = .false.

  ! Read command line options.
  call parse_options()

  ! Set data sky file type.
  select case(trim(sky_file_type_str))

     case(trim(FILE_TYPE_MAP))
        sky_file_type = S2_SKY_FILE_TYPE_MAP

     case(trim(FILE_TYPE_SKY))
        sky_file_type = S2_SKY_FILE_TYPE_SKY

     case default
        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
          's2fil_field_construct', comment_add='Invalid file type')

  end select

  ! Initialise filter field.
  field = s2fil_field_init(filename_filter, filename_sky, &
    sky_file_type, n_gamma)

  ! Compute filtered field coefficients.
  call s2fil_field_compute(field)

  ! Save field.  
  call s2fil_field_io_fits_write(filename_field, field, write_filter)

  ! Free memory.
  call s2fil_field_free(field)


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
            write(*,'(a,a)') 'Usage: s2fil_field_construct ', &
              '[-filter filename_filter]'
            write(*,'(a,a)') '                             ', &
              '[-data filename_sky]'
            write(*,'(a,a)') '                             ', &
              '[-filetype sky_file_type_str]'
            write(*,'(a,a)') '                             ', &
              '[-field filename_field]'
            write(*,'(a,a)') '                             ', &
              '[-ngamma n_gamma]'
            write(*,'(a,a)') '                             ', &
              '[-write_filter write_filter]'
            stop

          case ('-filter')
            filename_filter = trim(arg)

          case ('-data')
            filename_sky = trim(arg) 

          case ('-filetype')
            sky_file_type_str = trim(arg) 

          case ('-field')
            filename_field = trim(arg)

         case ('-ngamma')
            read(arg,*) n_gamma

         case ('-write_filter')
            read(arg,*) write_filter

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2fil_field_construct
