!------------------------------------------------------------------------------
! s2fil_draw_dots_only
!
!! Draw dots at specified positions on constant value map.
!!
!! Usage: 
!!   - [-help]: Display usage information.
!!   - [-nside nside]: Healpix nside of map to produce.
!!   - [-out filename_out]: Name of map HEALPix fits file to write.  Map 
!!     contains dots drawn at specified positions.
!!   - [-dots filename_dots]: Name of text file containing dot alpha and
!!     beta positions.
!!   - [-large large (optional)]: Logical to indicate whether to draw large
!!     dots.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2005
!
! Revisions:
!   November 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2fil_draw_dots_only

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2fil_field_mod, only: s2fil_field_io_txt_dots_read

  implicit none

  interface 
     function func(theta, phi, param) result(val)
       use s2_types_mod
       real(s2_sp), intent(in) :: theta, phi
       real(s2_sp), intent(in), optional :: param(:)
       real(s2_sp) :: val
     end function func
  end interface
  
  character(len=S2_STRING_LEN) :: filename_out, filename_dots
  integer :: nside
  type(s2_sky) :: sky
  real(s2_sp), allocatable :: alpha(:), beta(:)
  logical :: large = .true.

  ! Parse input parameters.
  call parse_options()

  ! Read alpha and beta values from file.
  call s2fil_field_io_txt_dots_read(filename_dots, alpha, beta)

  ! Initialise constant sky.
  sky = s2_sky_init(func, nside, S2_SKY_RING)

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
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: s2fil_draw_dots_only [-nside nside]'
            write(*,'(a)') '                            [-out filename_out]'
            write(*,'(a)') '                            [-dots filename_dots]'
            write(*,'(a)') '                            [-large large (optional)]'

            stop
          
          case ('-nside')
            read(arg,*) nside

          case ('-out')
            filename_out = trim(arg)

          case ('-dots')
            filename_dots = trim(arg)

          case ('-large')
            read(arg,*) large

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options




end program s2fil_draw_dots_only


function func(theta, phi, param) result(val)
  use s2_types_mod
  real(s2_sp), intent(in) :: theta, phi
  real(s2_sp), intent(in), optional :: param(:)
  real(s2_sp) :: val
  
  val = -1.6375e30
  
end function func
