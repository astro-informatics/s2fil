!------------------------------------------------------------------------------
! s2fil_draw_dots
!
!! Draw dots on sky map at specified positions.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of map HEAPix fits file to read.
!!   - [-ext ext (optional)]: Optional extension of HEALpix map fits file to 
!!     read map from.
!!   - [-out filename_out]: Name of map HEALPix fits file to write.  Map 
!!     contains dots drawn at specified positions.
!!   - [-dots filename_dots]: Name of text file containing dot alpha and
!!     beta positions.
!!   - [-large large (optional)]: Logical to indicate whether to draw large
!!     dots.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2fil_draw_dots

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2fil_field_mod, only: s2fil_field_io_txt_dots_read

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_out, filename_dots
  type(s2_sky) :: sky
  integer :: ext = 1
  real(s2_sp), allocatable :: alpha(:), beta(:)
  logical :: large = .true.

  ! Parse input parameters.
  call parse_options()

  ! Read alpha and beta values from file.
  call s2fil_field_io_txt_dots_read(filename_dots, alpha, beta)

  ! Initialse sky with map read in from map fits file.
  sky = s2_sky_init(filename_in, S2_SKY_FILE_TYPE_MAP, ext)

  ! Draw dots.
  call s2_sky_draw_dot(sky, alpha, beta, large)

  ! Save output file.
  call s2_sky_write_map_file(sky, filename_out)

  ! Free memory.
  call s2_sky_free(sky)
  deallocate(alpha, beta)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen 
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
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2_sky2map [-inp filename_in]'
            write(*,'(a)') '                   [-ext ext (optional)]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-dots filename_dots]'
            write(*,'(a)') '                   [-large large (optional)]'

            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-ext')
            read(arg,*) ext

          case ('-out')
            filename_out = trim(arg)

          case ('-dots')
            filename_dots = trim(arg)

          case ('-large')
            read(arg,*) large

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program s2fil_draw_dots
