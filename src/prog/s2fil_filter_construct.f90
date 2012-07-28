!------------------------------------------------------------------------------
! s2fil_filter_construct
!
!! Compute optimal filter from template and background process.
!!
!! Notes:
!!   - Complete for current version of s2fil_filter_mod but will no doubt 
!!     change as s2fil_filter_mod is extended to additional functionality.
!!     (Extensions planned: save full filter structure; functionality for SAF;
!!     ability to construct multiple dilations; dilatable filter; others)
!!     -- these additions have now been implemented.
!!
!! Usage: s2fil_filter_construct
!!   - [-help]: Display usage information.
!!   - [-out filename_filter]: Name of filter file to save.
!!   - [-nside nside]: Healpix nside resolution of construction.
!!   - [-filename_dil filename_dil]: Filename containing dilations to construct
!!     filter for. 
!!   - [-bkgnd_cmb filename_background_cmb]:
!!   - [-noise_var noise_var]: Variance of noise to add.
!!   - [-tmpl tmpl_type]: String defining template function type.
!!   - [-beam_fwhm beam_fwhm (arcmin)]: FWHM in arcmin of beam to apply.  If
!!     not present then no beam is applied.
!!   - [-filter_heu filter_heu]: Logical to specify whether heuristic 
!!     truncation of optimal filter alms is to be performed (recommended).
!!   - [-filter_type filter_type_string]: String defining filter type (either mf 
!!     or saf).
!!   - [-scale_type scale_type_string]: String defining scale type (either tmpl 
!!     or filter).
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.2 - April 2005
!
! Revisions:
!   February 2005 - Written by Jason McEwen 
!   April 2005 - Jason McEwen: Functionality added for SAF, multiple 
!     dilations, dilatable filter, saving of full filter structure.
!------------------------------------------------------------------------------

program s2fil_filter_construct

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2_vect_mod, only: s2_vect_arcmin_to_rad
  use s2fil_filter_mod
  use s2fil_error_mod
  use comb_tmplmap_mod   ! Templates should be defined in a s2 module in future,
                      ! but for now used ones defined when constructing simulations
!  use cswt_tmpl_mod   ! Templates should be defined in a sph module in future
  use cswt_tr_mod, only: cswt_tr_io_txt_dilation_read

  implicit none

#ifdef WMAP1
  integer, parameter :: NUM_COMMENT_LINES_BKGND_FILE = 29 !WMAP1
#endif
#ifdef WMAP3
  integer, parameter :: NUM_COMMENT_LINES_BKGND_FILE = 45   ! Value for WMAP3
#endif
#ifdef WMAP7
  integer, parameter :: NUM_COMMENT_LINES_BKGND_FILE = 0   ! Value for WMAP7
#endif

  character(len=S2_STRING_LEN) :: &
    filter_type_string = S2FIL_FILTER_TYPE_STR_MF
  character(len=S2_STRING_LEN) :: &
    scale_type_string = S2FIL_FILTER_SCALE_TYPE_STR_TMPL
  integer :: filter_type = S2FIL_FILTER_TYPE_MF
  integer :: scale_type = S2FIL_FILTER_SCALE_TYPE_TMPL
  
  character(len=S2_STRING_LEN) :: filename_background_cmb
  logical :: cmb_present = .false.
  logical :: noise_present = .false., beam_present = .false.
  real(s2_sp) :: noise_var = 0.0e0, beam_fwhm = 0.0e0
  character(len=S2_STRING_LEN) :: filename_dil
  character(len=S2_STRING_LEN) :: filename_filter
  character(len=S2_STRING_LEN), parameter :: &
    TMPL_TYPE_GAUSSIAN = 'gaussian', &
    TMPL_TYPE_MEXHAT = 'mexhat', &
    TMPL_TYPE_MORLET = 'morlet', &
    TMPL_TYPE_BUTTERFLY = 'butterfly', &
    TMPL_TYPE_BUBBLE = 'bubble', &
    TMPL_TYPE_TEXTURE = 'texture'
  character(len=S2_STRING_LEN) :: tmpl_type
  character(len=S2_STRING_LEN) :: tmpl_param_file
  logical :: param_file_present = .false.
  character(len=1), parameter :: COMMENT_CHAR = '#'
  character(len=S2_STRING_LEN) :: line, line2
  integer :: fileid_tmpl_param = 31
  real(s2_sp), allocatable :: tmpl_params(:)
  integer :: fail = 0
  integer :: n_tmpl_params, iparam

  integer :: nside, lmax, mmax, pix_scheme, background_read_lmin
  logical :: norm_pres_dil, background_read_scale
  logical :: filter_heu = .true.
  real(s2_sp), allocatable :: scale(:,:)
  type(s2_sky) :: tmpl
  type(s2_pl) :: background, background_cmb, background_noise
  type(s2_pl) :: beam
  type(s2fil_filter) :: filter


  ! --------------------------------------
  ! Set parameters
  ! --------------------------------------

  nside = 128
  pix_scheme = S2_SKY_RING
  norm_pres_dil = .false.
  filename_background_cmb = 'data_in/wmap_lcdm_pl_model_yr1_v1.txt'
  cmb_present = .true.
  noise_var = 0.0e0
  noise_present = .false.
  beam_fwhm = 0.0e0
  beam_present = .false.
  filename_filter = 'filter.fil'
  background_read_lmin = 2
  background_read_scale = .true.
  filter_heu = .true.

  tmpl_type = TMPL_TYPE_MEXHAT

  filter_type = S2FIL_FILTER_TYPE_MF
  scale_type = S2FIL_FILTER_SCALE_TYPE_TMPL

  ! Read command line options.
  call parse_options()

  lmax = 2*nside
  mmax = lmax

  ! Read dilations from file.
  call cswt_tr_io_txt_dilation_read(filename_dil, scale)


  ! --------------------------------------
  ! Initialise template
  ! --------------------------------------

  ! If template parameter file present, then read template parameters
  ! from file.
  if (param_file_present) then

     ! Open file
     open(fileid_tmpl_param, file=tmpl_param_file, &
          form='formatted', status='old')

     ! Ignore leading comment lines.
     line = COMMENT_CHAR
     do while(line(1:1) == COMMENT_CHAR)
        read(fileid_tmpl_param,'(a)') line
     end do

     ! Read number of template paramaters from last line read (that is
     ! actually not a comment line).
     read(line, *) line2, n_tmpl_params

     ! Allocate space for parameters.
     allocate(tmpl_params(1:n_tmpl_params), stat=fail)
     if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_filter_construct')
     end if

     ! Read data for sources found.
     do iparam = 1,n_tmpl_params
        read(fileid_tmpl_param,*) tmpl_params(iparam)
     end do

     ! Close file
     close(fileid_tmpl_param)

     ! Define template with parameters read from file.
     select case(trim(tmpl_type))

     case(trim(TMPL_TYPE_GAUSSIAN))

        tmpl = s2_sky_init(comb_tmplmap_gaussian, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_MEXHAT))

        tmpl = s2_sky_init(comb_tmplmap_mexhat, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_MORLET))

        tmpl = s2_sky_init(comb_tmplmap_morlet, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_BUTTERFLY))

        tmpl = s2_sky_init(comb_tmplmap_butterfly, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_BUBBLE))

        tmpl = s2_sky_init(comb_tmplmap_bubble, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_TEXTURE))

        tmpl = s2_sky_init(comb_tmplmap_texture, nside, pix_scheme, &
             lmax, mmax, param=tmpl_params, &
             fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case default

        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
             's2fil_filter_construct', comment_add='Invalid template type')

     end select

  else

     ! Define template with default parameters.
     select case(trim(tmpl_type))

     case(trim(TMPL_TYPE_GAUSSIAN))

        tmpl = s2_sky_init(comb_tmplmap_gaussian, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_MEXHAT))

        tmpl = s2_sky_init(comb_tmplmap_mexhat, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_MORLET))

        tmpl = s2_sky_init(comb_tmplmap_morlet, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_BUTTERFLY))

        tmpl = s2_sky_init(comb_tmplmap_butterfly, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_BUBBLE))

        tmpl = s2_sky_init(comb_tmplmap_bubble, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case(trim(TMPL_TYPE_TEXTURE))

        tmpl = s2_sky_init(comb_tmplmap_texture, nside, pix_scheme, &
             lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)

     case default

        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
             's2fil_filter_construct', comment_add='Invalid template type')

     end select

  end if

!  ! Uses cswt template definitions.
!  select case(trim(tmpl_type))
!
!     case(trim(TMPL_TYPE_GAUSSIAN))
!        
!        tmpl = s2_sky_init(cswt_tmpl_s2_gaussian, nside, pix_scheme, &
!          lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)
! 
!     case(trim(TMPL_TYPE_MEXHAT))
!        
!        tmpl = s2_sky_init(cswt_tmpl_s2_mexhat, nside, pix_scheme, &
!          lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)
!
!     case(trim(TMPL_TYPE_MORLET))
!        
!        tmpl = s2_sky_init(cswt_tmpl_s2_morlet, nside, pix_scheme, &
!          lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)
!
!     case(trim(TMPL_TYPE_BUTTERFLY))
!        
!        tmpl = s2_sky_init(cswt_tmpl_s2_butterfly, nside, pix_scheme, &
!          lmax, mmax, fun_type_in=S2_SKY_FUN_TYPE_SPHERE)
!
!     case default
!
!        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
!          's2fil_filter_construct', comment_add='Invalid template type')
!
!  end select


  ! --------------------------------------
  ! Initialise background
  ! --------------------------------------

  if(cmb_present) then
     background_cmb = s2_pl_init(filename_background_cmb, background_read_lmin, &
       lmax, NUM_COMMENT_LINES_BKGND_FILE, background_read_scale)
  end if


  ! --------------------------------------
  ! Apply beam if present
  ! --------------------------------------

  if(beam_present .and. cmb_present) then

     ! Convert beam_fwhm to radians.
     beam_fwhm = s2_vect_arcmin_to_rad(beam_fwhm)

     ! Construct beam.
     beam = s2_pl_init_guassian(beam_fwhm, lmax)

     ! Convolve background spectrum with beam.
     call s2_pl_conv(background_cmb, beam)

  end if


  ! --------------------------------------
  ! Add noise if required
  ! --------------------------------------

  ! Initialise noise if present.
  if(noise_present) then
     background_noise = s2_pl_init(noise_var, lmax)
  end if

  ! Add cmb and noise background if both present.
  if(cmb_present .and. noise_present) then
     background = s2_pl_add(background_cmb, background_noise)
  elseif(cmb_present) then
     background = s2_pl_init(background_cmb)
  elseif(noise_present) then
     background = s2_pl_init(background_noise)
  else
     call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
       's2fil_filter_construct', comment_add='Background not defined')
  end if

!!$! Print spectra to check units consistent.
!!$do el = 0,lmax
!!$   write(*,'(i5,4f12.5)') el, &
!!$        s2_pl_get_spec_l(background_cmb, el)*el*(el+1)/2/PI, &
!!$        s2_pl_get_spec_l(background_cmb, el), &
!!$        s2_pl_get_spec_l(background_noise, el), &
!!$        s2_pl_get_spec_l(background, el)
!!$end do

  ! When considering beam must convolve template at each scale with beam
  ! (done in s2fil_filter_comp_filter routine).


  ! --------------------------------------
  ! Initialise types
  ! --------------------------------------
 
  ! Set filter type.
  select case(trim(filter_type_string))

     case(trim(S2FIL_FILTER_TYPE_STR_MF))
        filter_type = S2FIL_FILTER_TYPE_MF

     case(trim(S2FIL_FILTER_TYPE_STR_SAF))
        filter_type = S2FIL_FILTER_TYPE_SAF

     case default
        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
          's2fil_filter_construct', comment_add='Invalid filter type string')

  end select

  ! Set scale type.
  select case(trim(scale_type_string))

     case(trim(S2FIL_FILTER_SCALE_TYPE_STR_TMPL))
        scale_type = S2FIL_FILTER_SCALE_TYPE_TMPL

     case(trim(S2FIL_FILTER_SCALE_TYPE_STR_FILTER))
        scale_type = S2FIL_FILTER_SCALE_TYPE_FILTER

     case default
        call s2fil_error(S2FIL_ERROR_PROG_INPUT_INVALID, &
          's2fil_filter_construct', comment_add='Invalid scale type string')

  end select


  ! --------------------------------------
  ! Construct and save filter
  ! --------------------------------------

  filter = s2fil_filter_init(tmpl, background, scale, &
    filter_type, scale_type, mmax, norm_pres_dil, beam, filter_heu)

  call s2fil_filter_io_fits_write(filename_filter, filter)


  ! --------------------------------------
  ! Free memory
  ! --------------------------------------

  deallocate(scale)
  call s2_sky_free(tmpl)
  if(cmb_present) call s2_pl_free(background_cmb)
  if(noise_present) call s2_pl_free(background_noise)
  if(beam_present) call s2_pl_free(beam)
  if(param_file_present) deallocate(tmpl_params)
  call s2_pl_free(background)
  call s2fil_filter_free(filter)


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
            write(*,'(a,a)') 'Usage: s2fil_filter_construct ', &
              '[-out filename_filter]'
            write(*,'(a,a)') '                              ', &
              '[-nside nside]'
            write(*,'(a,a)') '                              ', &
              '[-filename_dil filename_dil]'
            write(*,'(a,a)') '                              ', &
              '[-bkgnd_cmb filename_background_cmb]'
            write(*,'(a,a)') '                              ', &
              '[-noise_var noise_var]'
            write(*,'(a,a)') '                              ', &
              '[-tmpl tmpl_type]'
            write(*,'(a,a)') '                              ', &
              '[-tmpl_param tmpl_param_file]'
            write(*,'(a,a)') '                              ', &
              '[-beam_fwhm beam_fwhm (arcmin)]'
            write(*,'(a,a)') '                              ', &
              '[-filter_heu filter_heu]'
            write(*,'(a,a)') '                              ', &
              '[-filter_type filter_type_string]'
            write(*,'(a,a)') '                              ', &
              '[-scale_type scale_type_string]'
            stop

          case ('-out')
            filename_filter = trim(arg)
  
          case ('-nside')
            read(arg,*) nside

          case ('-filename_dil')
            filename_dil = trim(arg)

          case ('-bkgnd_cmb')
            filename_background_cmb = trim(arg)
            cmb_present = .true.

          case ('-noise_var')
            read(arg,*) noise_var
            noise_present = .true.

          case ('-tmpl')
            tmpl_type = trim(arg)

          case ('-tmpl_param')
            tmpl_param_file = trim(arg)
            param_file_present = .true.

          case ('-beam_fwhm')
            read(arg,*) beam_fwhm
            beam_present = .true.

          case ('-filter_heu')
            read(arg,*) filter_heu

          case ('-filter_type')
            filter_type_string = trim(arg)
         
         case ('-scale_type')
            scale_type_string = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program s2fil_filter_construct
