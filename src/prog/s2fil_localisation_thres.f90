!------------------------------------------------------------------------------
! s2fil_localisation_thres
!
!! Find localised regions in filtered field using thresholding strategy.  
!! Locations and values found are written to the standard output.
!!
!! Usage: s2fil_localisation_thres
!!   - [-help]: Display usage information.
!!   - [-field filename_field]: Name of precomputed field file to read.
!!   - [-filter filename_filter]: Name of filter file used to compute field.
!!     Read and stored in field object.  Must match field attributes, although
!!     not directly used in this program.
!!   - [-nsigma nsigma]: Number of sigma level to perform thresholding at.
!!   - [-thres filename_thres (optional)]: If present save cswt file of
!!     thresholded coefficients.
!!   - [-connected filename_connected (optional)]: If present save cswt file
!!     of connected components.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program s2fil_localisation_thres

  use s2_types_mod
  use s2_vect_mod
  use s2fil_types_mod
  use s2fil_error_mod
  use s2fil_field_mod

  implicit none

  integer, parameter :: MIN_REG_SIZE = 50
  integer, parameter :: S2FIL_RECLASS_SIZ = 1
  integer, parameter :: S2FIL_RECLASS_VAL = 2
  integer, parameter :: S2FIL_RECLASS_HEU = 3
  integer :: reclass_strategy = S2FIL_RECLASS_SIZ

  character(len=S2_STRING_LEN) :: filename_field, filename_filter
  character(len=S2_STRING_LEN) :: filename_out
  character(len=S2_STRING_LEN) :: filename_thres, filename_connected
  logical :: save_thres = .false., save_connected = .false.

  type(s2fil_field) :: field
  real(s2_sp), allocatable :: max_val(:,:)
  integer, allocatable :: max_loc(:,:,:)
  integer, allocatable :: max_siz(:,:)
  integer, allocatable :: n_regions(:)
  real(s2_sp) :: nsigma = 3.0e0
  integer :: i_reg, i_gamma0, irp, igp, fail=0
  integer :: n_alpha, n_beta, n_gamma
  logical :: no_gamma_search = .false.

  logical, allocatable :: keep(:,:)
  real(s2_sp) :: ANG_EXTEND_SIZE = 1.0e0
  real(s2_sp) :: dil(1:S2FIL_SCALE_DIM2_SIZE)
  real(s2_sp) :: alpha, beta, alpha_p, beta_p
  type(s2_vect) :: v, vp
  real(s2_sp) :: dist, dist_near, dot
  integer :: n_keep, i_keep

  character(len=1), parameter :: COMMENT_CHAR = '#'
  integer :: fileid = 21

  ! Read command line options.
  call parse_options()

  ! Read field from files.
  field = s2fil_field_init(filename_field, filename_filter)

  ! Set local variables to reuse.
  n_alpha = s2fil_field_get_n_alpha(field)
  n_beta = s2fil_field_get_n_beta(field)
  n_gamma = s2fil_field_get_n_gamma(field)

  ! Find local maxima.
  if(save_thres .and. save_connected) then
     call s2fil_field_loc_thres(field, nsigma, n_regions, max_val, max_loc, &
       max_siz, filename_thres, filename_connected)
  elseif(save_thres) then
     call s2fil_field_loc_thres(field, nsigma, n_regions, max_val, max_loc, &
       max_siz, filename_thres)
  elseif(save_connected) then
     call s2fil_field_loc_thres(field, nsigma, n_regions, max_val, max_loc, &
       max_siz, filename_connected=filename_connected)
  else
     call s2fil_field_loc_thres(field, nsigma, n_regions, max_val, max_loc, &
       max_siz)
  end if

  ! Refine candidate regions over orientations.

  ! Allocate keep matrix.
  allocate(keep(0:size(max_val,1)-1,1:size(max_val,2)), stat=fail)
  if(fail /= 0) then
     call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
          's2fil_localisation_thres')
  end if

  ! Set proximity ranges for alpha and beta (in radians).
  call s2fil_field_get_scale_val(field, 1, dil)
  dist_near = ANG_EXTEND_SIZE * 4.0e0 * atan(maxval(dil) / sqrt(2.0e0))
    ! Was bug here - bracket in wrong place.
    !  4 * atan(0.1/sqrt(2)) = 0.2824  (correct)
    !  4 * atan(0.1)/sqrt(2) = 0.2819  (before, but made minimal numerical difference)

  ! Consider each candidate, decide whether to keep or throw.
  do i_gamma0 = 0,s2fil_field_get_n_gamma(field)-1
     do i_reg = 1,n_regions(i_gamma0)

        ! If too small discard
        if(max_siz(i_gamma0, i_reg) < MIN_REG_SIZE) then

           keep(i_gamma0, i_reg) = .false.           

        else

           alpha = max_loc(i_gamma0,i_reg,1) / &
                real(n_alpha, s2_sp) * 2.0e0 * pi
           beta = max_loc(i_gamma0,i_reg,2) / &
                real(n_beta, s2_sp) * pi
           v = s2_vect_init(1e0, beta, alpha)
           call s2_vect_convert(v, S2_VECT_TYPE_CART)

           ! Initialise keep status to true, set to false if decide to throw.
           keep(i_gamma0, i_reg) = .true.

           ! Check all candidates found on other orientations.
           do igp = 0,s2fil_field_get_n_gamma(field)-1 
              if(igp == i_gamma0) cycle
              do irp = 1,n_regions(igp)
                 
                 ! Decide whether to keep or throw.

                 ! Compute angle between locations.                 
                 alpha_p = max_loc(igp,irp,1) / &
                      real(n_alpha, s2_sp) * 2.0e0 * pi
                 beta_p = max_loc(igp,irp,2) / &
                      real(n_beta, s2_sp) * pi
                 vp = s2_vect_init(1e0, beta_p, alpha_p)
                 dot = s2_vect_dot(v, vp)
                 dist = acos(dot)
                 call s2_vect_free(vp) 

                 if(dist < dist_near) then

                    ! Region nearby, so could potentially throw.

                    select case(reclass_strategy)

                       ! Choose based on size of region.
                       case (S2FIL_RECLASS_SIZ)

                          if( max_siz(igp, irp) > max_siz(i_gamma0, i_Reg) ) then
                             keep(i_gamma0, i_reg) = .false.
                             ! Loop done, but see out for simplicity of code.
                             ! (Otherwise require goto statements.)
                          end if

                       ! Choose based on max value of region.
                       case (S2FIL_RECLASS_VAL)

                          if( max_val(igp, irp) > max_val(i_gamma0, i_Reg) ) then
                             keep(i_gamma0, i_reg) = .false.
                             ! Loop done, but see out for simplicity of code.
                             ! (Otherwise require goto statements.)
                          end if
                         
                       ! Choose based on heuristic that weights size and max
                       ! value of region.
                       case (S2FIL_RECLASS_HEU)
stop 'heuristic reclassification strategy not done yet'
!** todo
                             ! Loop done, but see out for simplicity of code.
                             ! (Otherwise require goto statements.)

                       case default
                          call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
                               's2fil_localisation_thres', &
                               comment_add='Invalid reclassification strategy')
                 
                    end select


                 end if

              end do
           end do

           call s2_vect_free(v) 

        end if

     end do
  end do

  ! If want to consider only 0 orientation then set keep accordingly.
  if (no_gamma_search) then
     keep(0:size(max_val,1)-1,1:size(max_val,2)) = .false.
     keep(0,1:size(max_val,2)) = .true.
  end if

  ! Write local maxima data to the standard output.
  n_keep = 0
  do i_gamma0 = 0,s2fil_field_get_n_gamma(field)-1
     do i_reg = 1,n_regions(i_gamma0)

        write(*,'(a,l16)') 'keep: ', keep(i_gamma0, i_reg)
        write(*,'(a,i12)') 'i_gamma0: ', i_gamma0
        write(*,'(a,i12)') 'i_reg:    ', i_reg 
        write(*,'(a,f14.4)') 'max_val:', max_val(i_gamma0,i_reg)
        write(*,'(a,i9,a,i4)') 'max_loc:', max_loc(i_gamma0,i_reg,1), ',', max_loc(i_gamma0,i_reg,2) 
        write(*,'(a,i14)') 'max_siz:', max_siz(i_gamma0,i_reg)
        write(*,'(a,f16.4)') 'alpha:',  max_loc(i_gamma0,i_reg,1) / &
             real(n_alpha, s2_sp) * 2.0e0 * pi
        write(*,'(a,f16.4)') 'beta: ',  max_loc(i_gamma0,i_reg,2) / &
             real(n_beta, s2_sp) * pi
        write(*,*)
        
        if (keep(i_gamma0, i_reg)) then
           n_keep = n_keep + 1
        end if

     end do
  end do

  ! Write localisation data to file.  
  open(unit=fileid, file=filename_out, status='replace', action='write')
  write(fileid,'(a,a)') COMMENT_CHAR, ' Localised source positions'
  write(fileid,'(a,a)') COMMENT_CHAR, ' Written by s2fil-0.1'
  write(fileid,'(a)') COMMENT_CHAR
  write(fileid,'(a,i7)') 'n_sources= ', n_keep

  i_keep = 0
  do i_gamma0 = 0,s2fil_field_get_n_gamma(field)-1
     do i_reg = 1,n_regions(i_gamma0)

        if(keep(i_gamma0,i_reg)) then

           i_keep = i_keep + 1
           write(fileid, '(a)') COMMENT_CHAR
           write(fileid,'(a,f7.4)') 'amplitude= ', max_val(i_gamma0,i_reg)
           write(fileid,'(a,f11.4)') 'alpha= ', max_loc(i_gamma0,i_reg,1) / &
                real(n_alpha, s2_sp) * 2.0e0 * pi
           write(fileid,'(a,f11.4)') 'beta=  ', max_loc(i_gamma0,i_reg,2) / &
                real(n_beta, s2_sp) * pi
           write(fileid,'(a,f11.4)') 'gamma= ', &
                i_gamma0 / real(n_gamma, s2_sp) * 2e0 * pi

        end if

     end do
  end do

  close(fileid)

  ! Write output text file in dots format. :: DON'T WRITE DOTS ANYMORE, WRITE FULL INFO INSTEAD
!  call s2fil_field_io_txt_dots_write(filename_out, &
!    max_loc(0,:,1) / real(s2fil_field_get_n_alpha(field), s2_sp)*2.0e0*pi, &
!    max_loc(0,:,2) / real(s2fil_field_get_n_beta(field), s2_sp) * pi)

  ! Free memory used.
  call s2fil_field_free(field)
  deallocate(max_val, max_loc, max_siz, n_regions, keep)


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
            write(*,'(a,a)') 'Usage: s2fil_localisation_thres ', &
              '[-field filename_field]'
            write(*,'(a,a)') '                                ', &
              '[-filter filename_filter]'
            write(*,'(a,a)') '                                ', &
              '[-nsigma nsigma]'
            write(*,'(a,a)') '                                ', &
              '[-no_gamma_search no_gamma_search]'
            write(*,'(a,a)') '                                ', &
              '[-thres filename_thres (optional)]'
            write(*,'(a,a)') '                                ', &
              '[-connected filename_connected (optional)]'
            stop

          case ('-field')
            filename_field = trim(arg) 

          case ('-filter')
            filename_filter = trim(arg)

         case ('-nsigma')
            read(arg,*) nsigma

         case ('-no_gamma_search')
            read(arg,*) no_gamma_search

          case ('-thres')
            filename_thres = trim(arg)
            save_thres = .true.

          case ('-connected')
            filename_connected = trim(arg)
            save_connected = .true.

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2fil_localisation_thres
