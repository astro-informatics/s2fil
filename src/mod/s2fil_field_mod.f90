!------------------------------------------------------------------------------
! s2fil_field_mod -- S2FIL library field class
!
!! Functionality to compute and store the coefficients of a filtered field from
!! the optimal filter and data sky map.  Interfaces with the CSWT library to
!! actually compute the `wavelet' coefficients.
!!
!! Notes:
!!   - Note dilation values stored in tr structures are all unitary, since 
!!     filters already dilated to appropriate size.  To find actual 
!!     corresponding dilation look in filter%scale array.

!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 April 2005
!
! Revisions:
!   April 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module s2fil_field_mod

  use s2_types_mod
  use s2_sky_mod
  use s2fil_error_mod
  use s2fil_types_mod
  use s2fil_filter_mod
  use cswt_tr_mod
  use cswt_swav_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    s2fil_field_init, &
    s2fil_field_free, &
    s2fil_field_compute, &
    s2fil_field_loc_thres, &
    s2fil_field_io_fits_write, &
    s2fil_field_io_txt_dots_write, &
    s2fil_field_io_txt_dots_read, &
    s2fil_field_get_init, &
    s2fil_field_get_n_gamma, &
    s2fil_field_get_computed_status, &
    s2fil_field_get_n_alpha, &
    s2fil_field_get_n_beta, &
    s2fil_field_get_lmax, &
    s2fil_field_get_mmax, &
    s2fil_field_get_n_scale, &
    s2fil_field_get_nside, &
    s2fil_field_get_scale_val  


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface s2fil_field_init
     module procedure &
       s2fil_field_init_data, & 
       s2fil_field_init_file, &
       s2fil_field_init_copy
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! None.


  !---------------------------------------
  ! Data types
  !---------------------------------------
 
  type, public :: s2fil_field
     private
     logical :: init = .false.
     integer :: n_gamma = 0
     type(s2fil_filter) :: filter
     type(s2_sky) :: sky
     type(cswt_tr), allocatable :: tr(:)
     logical :: computed_status = .false.
  end type s2fil_field


  !----------------------------------------------------------------------------
  
  contains


    !--------------------------------------------------------------------------
    ! s2fil_field_init_data
    !
    !! Initialise an uncomputed s2fil filtered field from data contained in 
    !! files.  Note that space for the filtered field coefficient tr 
    !! structures is allocated here but the tr structures themselves are not
    !! computed yet. 
    !!
    !! Variables:
    !!   - filename_filter: Name of full s2fil_filter fits file containing all
    !!     optimal filter data (including filenames of associated filter 
    !!     sky, template and background files).
    !!   - filename_sky: Name of data sky file (either fits map file or fits
    !!     full s2_sky file as defined by sky_file_type).
    !!   - sky_file_type: Type of data sky file read (options specified in 
    !!     s2_sky_mod as either S2_SKY_FILE_TYPE_MAP or 
    !!     S2_SKY_FILE_TYPE_SKY).
    !!   - n_gamma: N_gamma (number of orientations) considered for the 
    !!     filtered field map.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_init_data(filename_filter, filename_sky, &
      sky_file_type, n_gamma) result(field)

      character(len=*), intent(in) :: filename_filter
      character(len=*), intent(in) :: filename_sky
      integer, intent(in) :: sky_file_type
      integer, intent(in) :: n_gamma
      type(s2fil_field) :: field

      integer :: i_scale, n_scale, lmax, mmax, fail
      real(s2_sp) :: dilation_unit(1,CSWT_SWAV_DILATION_DIM) = 1.0e0
      type(s2_sky) :: temp_sky
      type(cswt_swav) :: swav

      ! Check object not already initialised.
      if(field%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_field_init_data')
        return
      end if

      ! Initialise filter from files.
      field%filter = s2fil_filter_init(filename_filter)

      ! Read sky to analyse.
      field%sky = s2_sky_init(filename_sky, sky_file_type)

      ! Initialise field n_gamma variable.
      field%n_gamma = n_gamma

      ! Set local size variables.
      n_scale = s2fil_filter_get_n_scale(field%filter)
      lmax = s2fil_filter_get_lmax(field%filter)
      mmax = s2fil_filter_get_mmax(field%filter)

      ! Allocate space for tr structures.
      allocate(field%tr(1:n_scale), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, 's2fil_field_init_data')
      end if

      ! Initialse tr structures for each scale.
      do i_scale = 1,n_scale

         ! Initialse wavelet from filter sky.
         temp_sky = s2fil_filter_get_sky(field%filter, i_scale)
         swav = cswt_swav_init(temp_sky)
         call s2_sky_free(temp_sky)

         ! Initialise tr structure with unit dilations.
         ! (Filters already dilated so no need for further dilation).
         field%tr(i_scale) &
           = cswt_tr_init(swav, dilation_unit, lmax, mmax, field%n_gamma)

         ! Free spherical wavelet object used to initialise tr structure.
         call cswt_swav_free(swav)

      end do

      ! Set status flags.
      field%init = .true.
      field%computed_status = .false.

      ! Check filters and sky consistent sizes.
      ! Set sky sizes to same as filter if undefined.
      ! Must wait until object initialised before running check.
      call s2fil_field_ckset_size(field)

    end function s2fil_field_init_data


    !--------------------------------------------------------------------------
    ! s2fil_filter_init_file
    !
    !! Initialise a field structure from a field fits file.  Wrapper for 
    !! routine s2fil_field_io_fits_read.
    !!
    !! Notes:
    !!  - Much of the field data is not stored directly in the s2fil field 
    !!    fits file but rather in associated child files, the filenames of 
    !!    which are stored in the field file header.
    !!
    !! Variables:
    !!  - filename: The fits file to read the field from.
    !!  - field: Filtered field constructed with the data read
    !!    from the file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_init_file(filename, filename_filter) result(field)

      character(len=*), intent(in) :: filename
      type(s2fil_field) :: field
      character(len=*), intent(in), optional :: filename_filter

      ! Check object not already initialised.
      if(field%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_field_init_file')
        return
      end if

      ! Read fiield file.
      call s2fil_field_io_fits_read(filename, field, filename_filter)

      ! All status flags set in s2fil_filter_io_fits_read routine.

    end function s2fil_field_init_file


    !--------------------------------------------------------------------------
    ! s2fil_field_init_copy
    !
    !! Initialse a field structure as a copy of another field.
    !!
    !! Variables:
    !!   - orig: The original field to copy.
    !!   - copy: The initialised field copied from the original.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_init_copy(orig) result(copy)

      type(s2fil_field), intent(in) :: orig
      type(s2fil_field) :: copy

      integer :: n_scale, i_scale, fail

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_field_init_copy')
      end if

      ! Copy object attributes.

      copy%n_gamma = orig%n_gamma
      copy%sky = s2_sky_init(orig%sky)
      copy%filter = s2fil_filter_init(orig%filter)
      copy%computed_status = orig%computed_status

      ! Get local n_scale variable.
      n_scale = s2fil_field_get_n_scale(orig)

      ! Allocate space for tr structures.
      allocate(copy%tr(1:n_scale), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, 's2fil_field_init_data')
      end if
      
      ! Initialse tr structures as copies for each scale.
      do i_scale = 1,n_scale
         copy%tr(i_scale) = cswt_tr_init(orig%tr(i_scale))
      end do

      ! Set copy as initialised.
      copy%init = .true.

    end function s2fil_field_init_copy


    !--------------------------------------------------------------------------
    ! s2fil_field_free
    !
    !! Free all data associated with an initialised field and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - field: The s2fil field object to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_free(field)
      
      type(s2fil_field), intent(inout) :: field

      integer :: n_scale, i_scale

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_free')
      end if

      n_scale = s2fil_field_get_n_scale(field)

      ! Free cswt tr structures for each scale.
      do i_scale = 1,n_scale
         call cswt_tr_free(field%tr(i_scale))
      end do

      ! Free other structures.
      call s2fil_filter_free(field%filter)
      call s2_sky_free(field%sky)

      ! Reset other variables.
      field%n_gamma = 0
      field%computed_status = .false.
      field%init = .false.

    end subroutine s2fil_field_free


    !--------------------------------------------------------------------------
    ! s2fil_field_compute
    !
    !! Compute filter field wavelet coefficient tr structures.
    !!
    !! Notes:
    !!   - Sky alms are computed before calling cswt_tr_analysis, so that 
    !!     they need only be computed once for all scales.
    !!   - Object must already be initialised with space for tr structures 
    !!     allocated for all scales.
    !!
    !! Variables:
    !!   - field: Field to compute tr structures for.  
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_compute(field)

      type(s2fil_field), intent(inout) :: field

      integer :: i_scale
      
      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_compute')
        return
      end if
      
      ! Check filtered field not already computed.
      if(field%computed_status) then
        call s2fil_error(S2FIL_ERROR_FIELD_TR_DEF, 's2fil_field_compute')
        return
      end if

      ! Compute alms for sky so only need to compute once.
      ! (Sky lmax and mmax, if not initially defined, will have been set to
      ! filter size values by s2fil_field_ckset_size.)
      call s2_sky_compute_alm(field%sky)

      ! Perform filtering for each scale.
      do i_scale = 1,s2fil_field_get_n_scale(field)
         call cswt_tr_analysis(field%tr(i_scale), field%sky, method=CSWT_TR_METHOD_FAST_FFT)
      end do

      ! Set field as computed.
      field%computed_status = .true.

    end subroutine s2fil_field_compute


    !--------------------------------------------------------------------------
    ! s2fil_field_ckset_size
    !
    !! Check and set field sizes.  Ensure filter and sky sizes are consistent.
    !! If sky sizes are not defined (as most likely for lmax and mmax if just
    !! read from map file) then set to sizes defined by filter structure.
    !!
    !! Variables:
    !!   - field: Field to check filter and sky sizes and set if necessary.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_ckset_size(field)

      type(s2fil_field), intent(inout) :: field

      integer :: nside_filter, lmax_filter, mmax_filter
      integer :: nside_sky, lmax_sky, mmax_sky
      
      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_ckset_size')
      end if
      
      ! Get filter sizes.
      nside_filter = s2fil_filter_get_nside(field%filter)
      lmax_filter = s2fil_filter_get_lmax(field%filter)
      mmax_filter = s2fil_filter_get_mmax(field%filter)

      ! Get sky sizes.
      nside_sky = s2_sky_get_nside(field%sky)
      lmax_sky = s2_sky_get_lmax(field%sky)
      mmax_sky = s2_sky_get_mmax(field%sky)

      ! If size variables not defined for sky then set same as filter.
      if(nside_sky == 0) then
         nside_sky = nside_filter
         call s2_sky_set_nside(field%sky, nside_filter)
      end if
      if(lmax_sky == 0 .and. mmax_sky == 0) then
         lmax_sky = lmax_filter
         mmax_sky = mmax_filter
         call s2_sky_set_lmax(field%sky, lmax_filter, mmax_filter)
      end if

      ! Check sizes consistent.
      if(nside_filter /= nside_sky &
         .or. lmax_filter /= lmax_sky &
         .or. mmax_filter /= mmax_sky) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, 's2fil_field_ckset_size', &
           comment_add='Filter and data sky sizes are inconsistent')
      end if

    end subroutine s2fil_field_ckset_size


    !--------------------------------------------------------------------------
    ! s2fil_field_loc_thres
    !
    !! Find localised regions in filtered field using thresholding strategy.
    !!
    !! Notes:
    !!   - *Currently only set up to consider only one dilation 
    !!     and orientation.*
    !!   - Data arrays of maximum local values and locations are allocated
    !!     herein, since no prior knowledge of number of objects contained 
    !!     in map.  These arrays *must* be freed by the calling routine.
    !!   - Sigma computed across all scales and orientations.
    !!
    !! Variables:
    !!   - field: Field containing filtered coefficients to detect local 
    !!     maxima in.
    !!   - nsigma: Number of sigma level to perform thresholding at.
    !!   - max_val: List of local maximum values found.
    !!   - max_loc: List of locations of local maximum values found.
    !!   - [filename_thres]: If present save cswt file of thresholded
    !!     coefficients.
    !!   - [filename_connected]: If present save cswt file of connected 
    !!     components.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_loc_thres(field, nsigma, n_regions, &
         max_val, max_loc, max_siz, filename_thres, filename_connected)

      type(s2fil_field), intent(in) :: field
      real(s2_sp), intent(in) :: nsigma
      integer, allocatable, intent(out) :: n_regions(:)
      real(s2_sp), allocatable, intent(out) :: max_val(:,:)
      integer, allocatable, intent(out) :: max_loc(:,:,:)
      integer, allocatable, intent(out) :: max_siz(:,:)
      character(len=*), intent(in), optional :: filename_thres
      character(len=*), intent(in), optional :: filename_connected

      type(cswt_tr) :: tr_merg
!      real(s2_sp) :: dil(1:S2FIL_SCALE_DIM2_SIZE)
!      real(s2_sp) :: morph_dil_angle
!      real(s2_sp), parameter :: MORPH_DIL_EXTEND_SIZE = 1.0e0

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_loc_thres')
      end if

      ! Check tr structure computed.
      if(.not. field%computed_status) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_loc_thres', &
          comment_add='Tr strucutres not computed')
      end if

      ! Merge tr structures into one so can consider all dilations at once
      ! in cswt_tr module.
      tr_merg = cswt_tr_init(field%tr)

      ! Threshold merged tr structure.
      call cswt_tr_wcoeff_thres_nsigma(tr_merg, nsigma, &
           CSWT_TR_THRES_NSIGMA_MODE_ABOVE)

      ! Save thresholded coefficients if filename_thres present.
      if(present(filename_thres)) then
         call cswt_tr_io_fits_write_wcoeff(filename_thres, tr_merg)
      end if

      ! Find values and locations of local maxima.
      call cswt_tr_localmax_ab(tr_merg, 1, n_regions, max_val, max_loc, &
           max_siz, filename_connected)

      ! Find across orietations also (instead of above call).
      ! Not completed.  Just compare distances to all posistions found
      ! using above technique.
!      call s2fil_filter_get_scale(field%filter, 1, dil)
!      morph_dil_angle = MORPH_DIL_EXTEND_SIZE * &
!           4.0e0 * atan(maxval(dil) / sqrt(2.0e0))
!      write(*,*) 'morph_dil_angle: ', morph_dil_angle
!      call cswt_tr_localmax_abg(tr_merg, 1, n_regions, max_val, max_loc, &
!           morph_dil_angle, filename_connected)

      ! Free temporary merged tr structure.
      call cswt_tr_free(tr_merg)

    end subroutine s2fil_field_loc_thres


    !--------------------------------------------------------------------------
    ! File IO routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2fil_field_io_fits_write
    !
    !! Write filtered field structure to a fits file.
    !!
    !! Notes:
    !!   - The field file only contains a header and no data directly (except
    !!     scales).  All other data is stored in associated child files, the
    !!     filenames of which are stored in the filtered field file header.
    !!   - An option is given for saving the filter file since this will most
    !!     likely already have been read from an input file and so, to avoid 
    !!     duplication, another output need not be written.
    !!   - The sky data file is written (and may be duplicated), however, 
    !!     since this will have alms computed and will be saved in the full 
    !!     s2_sky fits format.
    !!
    !! Variables:
    !!  - filename: Name of fits file to write to.
    !!  - field: The filtered field structure containing all data to be 
    !!    written to fits file.
    !!  - [comment]: Comment to append to header of output fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_fits_write(filename, field, write_filter, &
      comment)

      character(len=*), intent(in) :: filename
      type(s2fil_field), intent(in) :: field
      logical, intent(in), optional :: write_filter
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      logical :: simple, extend, file_exists
      integer :: naxis
      integer :: naxes(1)
      integer :: tfields, nrows, varidat
      character(len=32) :: ttype(2), tform(2), tunit(2), extname
      integer :: frow, felem, colnum

      character(len=S2_STRING_LEN) :: comment_child_file, filename_sky, &
           filename_tr, file_key, filename_filter
      integer :: i_scale, n_scale, fail
      logical :: write_filter_status = .false.
      real(s2_sp), allocatable :: scales(:,:)      

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_io_fits_write')
      end if 

      ! Check tr structure computed.
      if(.not. field%computed_status) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_io_fits_write', &
          comment_add='Tr strucutres not computed')
      end if

      ! Set write filter status.
      if(present(write_filter)) write_filter_status = write_filter
      
      ! Define FITS parameters.

      bitpix=-32 ! Real single precision.
      status=0   ! Initialse error status to zero.

      ! Check if file already exists.
      call s2fil_field_io_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_EXISTS, &
              's2fil_field_io_fits_write')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.

      call ftpcom(unit, &
        '  S2fil field file created by s2fil-0.1',status)
      call ftpcom(unit, &
        '  Primary extension empty',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if

      call ftpkyj(unit,'NSIDE', s2fil_field_get_nside(field), &
        'max spherical harmonic l considered',status)
      call ftpkyj(unit,'LMAX', s2fil_field_get_lmax(field), &
        'max spherical harmonic l considered',status)
      call ftpkyj(unit,'MMAX', s2fil_field_get_mmax(field), &
        'max spherical harmonic m considered',status)
      call ftpkyj(unit,'NSCALE', s2fil_field_get_n_scale(field), &
        'number of scales considered',status)
      call ftpkyj(unit,'NALPHA', s2fil_field_get_n_alpha(field), &
        'max spherical harmonic m considered',status)
      call ftpkyj(unit,'NBETA', s2fil_field_get_n_beta(field), &
        'max spherical harmonic m considered',status)
      call ftpkyj(unit,'NGAMMA', s2fil_field_get_n_gamma(field), &
        'max spherical harmonic m considered',status)
      call ftpkyl(unit,'NORMDIL', &
        s2fil_filter_get_norm_pres_dil(field%filter), &
        'norm preserving dilation status',status)

      call ftpkyj(unit,'FTYPE', &
        s2fil_filter_get_filter_type(field%filter), &
        'filter type flag',status)
      select case(s2fil_filter_get_filter_type(field%filter)) 
         case(S2FIL_FILTER_TYPE_MF)
            call ftpkys(unit,'FTYPESTR', S2FIL_FILTER_TYPE_STR_MF, &
              'string describing filter type',status)            
         case(S2FIL_FILTER_TYPE_SAF)
            call ftpkys(unit,'FTYPESTR', S2FIL_FILTER_TYPE_STR_SAF, &
              'string describing filter type',status)
      end select

      call ftpkyj(unit,'STYPE', &
        s2fil_filter_get_scale_type(field%filter), &
        'scale type flag',status)
      select case(s2fil_filter_get_scale_type(field%filter)) 
         case(S2FIL_FILTER_SCALE_TYPE_TMPL)
            call ftpkys(unit,'STYPESTR', &
              S2FIL_FILTER_SCALE_TYPE_STR_TMPL, &
              'string describing scale type',status)            
         case(S2FIL_FILTER_SCALE_TYPE_FILTER)
            call ftpkys(unit,'STYPESTR', &
               S2FIL_FILTER_SCALE_TYPE_STR_FILTER, &
              'string describing scale type',status)
      end select

      ! Get n_scale for multiple local use below.
      n_scale = s2fil_field_get_n_scale(field)

      ! Write child files.

      write(comment_child_file, '(a,a)') &
        '  Part of parent file structure: ', &
        trim(filename)

      ! Write data sky.
      write(filename_sky,'(a,a)') filename(1:len(trim(filename))-4), &
        '_data.sky'
      call s2_sky_io_fits_write(filename_sky, field%sky, &
        trim(comment_child_file))
      call ftpkys(unit,'SKY', trim(filename_sky), &
        'name of data sky file',status)  

      ! Write tr structures.
      do i_scale = 1,n_scale

         write(filename_tr,'(a,a,i2.2,a)') &
           filename(1:len(trim(filename))-4), &
           '_tr', i_scale, '.cswt'

         call cswt_tr_io_fits_write_wcoeff(filename_tr, &
           field%tr(i_scale), trim(comment_child_file))

         write(file_key, '(a,i2.2)') 'TR', i_scale
         call ftpkys(unit, trim(file_key), trim(filename_tr), &
              'name of filtered field coefficient file', status)

      end do

      ! Write filter (optional).
      if(write_filter_status) then

         write(filename_filter,'(a,a)') filename(1:len(trim(filename))-4), &
           '_filter.fil'

         call s2fil_filter_io_fits_write(filename_filter, field%filter, &
           trim(comment_child_file))

         call ftpkys(unit,'FILTER', trim(filename_filter), &
           'name of filter file',status)  

      else

         call ftpkys(unit,'FILTER', 'Not saved', &
           'name of filter file',status)  

      end if

      ! Save scales in new extension.

      ! Insert binary table extension for scales.
      extname='SCALE'
      ttype(1)='SC1'
      ttype(2)='SC2'
      tform(1)='1E'
      tform(2)='1E'
      tunit(1)='direct (rad)'
      tunit(2)='direct (rad)'
      tfields=2
      nrows=n_scale
      varidat=0
      call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

      ! Get scales that are stored in filter.
      allocate(scales(1:n_scale, S2FIL_SCALE_DIM2_SIZE), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
           's2fil_field_io_fits_write')
      end if
      call s2fil_filter_get_scale(field%filter, scales)

      ! Write dilations to binary table.
      frow=1
      felem=1
      colnum=1
      call ftpcle(unit,colnum,frow,felem,nrows,scales(:,1),status)
      colnum=2
      call ftpcle(unit,colnum,frow,felem,nrows,scales(:,2),status)

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2fil_field_io_fits_error_check(status, .true.)

      ! Free temporary memory.
      deallocate(scales)

    end subroutine s2fil_field_io_fits_write


    !--------------------------------------------------------------------------
    ! s2fil_field_io_fits_read
    ! 
    !! Read a s2fil_field fits file and all associated child files (i.e. data
    !! sky, tr structure coefficients and filter files.)
    !!
    !! Notes:
    !!   - Child files (as described in the filter file header) must be 
    !!     located in the same directory.
    !!
    !! Variables:
    !!   - filename: Name of primary s2fil_filter fits file.
    !!   - filter: Returned filter structure initialised with the data 
    !!     contained in the filter fits file and the associated child files.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------


    subroutine s2fil_field_io_fits_read(filename, field, filename_filter_in)

      character(len=*), intent(in) :: filename
      type(s2fil_field), intent(out) :: field
      character(len=*), intent(in), optional :: filename_filter_in

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: ihdu, hdutype, naxis
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: colnum, frow, felem, nelem
      real(s2_sp) :: nullval

      integer :: nside, lmax, mmax, n_scale, n_alpha, n_beta, n_gamma
      integer :: filter_type, scale_type, i_scale, fail
      logical :: norm_pres_dil
      character(len=S2_STRING_LEN) :: filename_sky, filename_tr, &
        filename_filter, file_key
      real(s2_sp), allocatable :: scales_field(:,:), scales_filter(:,:)

      real(s2_sp), parameter :: TOL = 1.0e-5

      ! Check object not already initialised.
      if(field%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_field_io_fits_read')
        return
      end if

      ! Initialse error status to zero.
      status=0

      ! Check if file already exists.
      call s2fil_field_io_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)


      !---------------------------------------
      ! Read primary header variables
      !---------------------------------------

      call ftgkyj(unit, 'NSIDE', nside, comment, status)
      call ftgkyj(unit, 'LMAX', lmax, comment, status)
      call ftgkyj(unit, 'MMAX', mmax, comment, status)

      call ftgkyj(unit, 'NSCALE', n_scale, comment, status)
      call ftgkyj(unit, 'NALPHA', n_alpha, comment, status)
      call ftgkyj(unit, 'NBETA', n_beta, comment, status)
      call ftgkyj(unit, 'NGAMMA', n_gamma, comment, status)
      call ftgkyl(unit, 'NORMDIL', norm_pres_dil, comment, status)

      call ftgkyj(unit, 'FTYPE', filter_type, comment, status)
      call ftgkyj(unit, 'STYPE', scale_type, comment, status)

      ! Check filter_type valid.
      if(filter_type /= S2FIL_FILTER_TYPE_MF &
        .and. filter_type /= S2FIL_FILTER_TYPE_SAF) then
         call s2fil_error(S2FIL_ERROR_FILTER_TYPE_INVALID, &
           's2fil_field_io_fits_read')
      end if

      ! Check scale_type valid.
      if(scale_type /= S2FIL_FILTER_SCALE_TYPE_TMPL &
        .and. scale_type /= S2FIL_FILTER_SCALE_TYPE_FILTER) then
         call s2fil_error(S2FIL_ERROR_SCALE_TYPE_INVALID, &
           's2fil_field_io_fits_read')
      end if

      ! Check correct number of HDUs in input file.
      hdunum = 2   ! Primary header plus scales.
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
       call s2fil_error(S2FIL_ERROR_SCALE_TYPE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='Invalid number of extensions')
      end if

      ! Set n_gamma
      field%n_gamma = n_gamma

      !---------------------------------------
      ! Read data sky and tr structures
      !---------------------------------------

      ! Get data sky filename from header then read data from 
      ! sky file.
      call ftgkys(unit, 'SKY', filename_sky, comment, status)
      field%sky = s2_sky_init(filename_sky, S2_SKY_FILE_TYPE_SKY)

      ! Check sizes consistent.
      if(s2_sky_get_nside(field%sky) /= nside &
           .or. s2_sky_get_lmax(field%sky) /= lmax &
           .or. s2_sky_get_mmax(field%sky) /= mmax) then
        call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
          's2fil_field_io_fits_read', &
          comment_add='Data sky size inconsistent')
      end if

      ! Allocate space for tr structure for each scale.
      allocate(field%tr(1:n_scale), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
           's2fil_field_io_fits_read')
      end if

      ! Read tr structure fits files containing coefficients 
      ! for each scale.
      do i_scale = 1,n_scale

         ! Get tr structure filename from field fits file header.
         write(file_key,'(a,i2.2)') 'TR', i_scale
         call ftgkys(unit, trim(file_key), filename_tr, comment, status)

         ! Read tr structure coefficients from file.
         field%tr(i_scale) = cswt_tr_init(filename_tr)

         ! Check sizes consistent.
         if(cswt_tr_get_lmax(field%tr(i_scale)) /= lmax &
              .or. cswt_tr_get_mmax(field%tr(i_scale)) /= mmax &
              .or. cswt_tr_get_n_alpha(field%tr(i_scale)) /= n_alpha &
              .or. cswt_tr_get_n_beta(field%tr(i_scale)) /= n_beta &
              .or. cswt_tr_get_n_gamma(field%tr(i_scale)) /= n_gamma) then
            call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
              's2fil_field_io_fits_read', &
              comment_add='Tr structure size inconsistent')
        end if

      end do
        
 
      !---------------------------------------
      ! Read filter
      !---------------------------------------

      ! Get appropriate filter filename, either from passed argument 
      ! or file header.

      if(present(filename_filter_in)) then

         ! Read filename from s2fil field file and check filter file not 
         ! specified there.
         call ftgkys(unit, 'FILTER', filename_filter, comment, status)
         if(trim(filename_filter) /= 'Not saved') then
            call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
              's2fil_field_io_fits_read', &
              comment_add='Filter specified in field file but read from elsewhere')
         end if

         ! Set filter filename as passed argument.
         filename_filter = filename_filter_in

      else

         ! Read filter filename from field file header.
         call ftgkys(unit, 'FILTER', filename_filter, comment, status)

         ! Check valid filename, i.e. not 'Not saved.'
         if(trim(filename_filter) == 'Not saved') then
            call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
              's2fil_field_io_fits_read', &
              comment_add='Filter not specified in field file')
         end if

      end if

      ! Read filter from file.
      field%filter = s2fil_filter_init(filename_filter)

      ! Check filter sizes and variables consistent.
      if(s2fil_filter_get_lmax(field%filter) /= lmax &
           .or. s2fil_filter_get_mmax(field%filter) /= mmax &
           .or. s2fil_filter_get_nside(field%filter) /= nside &
           .or. s2fil_filter_get_n_scale(field%filter) /= n_scale &
           .or. s2fil_filter_get_filter_type(field%filter) /= filter_type &
           .or. s2fil_filter_get_scale_type(field%filter) /= scale_type &
           .or. s2fil_filter_get_norm_pres_dil(field%filter) .neqv. norm_pres_dil) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
              's2fil_field_io_fits_read', &
              comment_add='Filter inconsistent with field file header variables')
      end if


      !---------------------------------------
      ! Read scales from next extension
      !---------------------------------------

      ! Allocate space for scales, both those stored in the field header and 
      ! those stored in the filter.
      allocate(scales_field(n_scale, S2FIL_SCALE_DIM2_SIZE), stat=fail)
      allocate(scales_filter(n_scale, S2FIL_SCALE_DIM2_SIZE), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
            's2fil_field_io_fits_read')
      end if

      ! Move to second extension (binary table containing scale values).
      ihdu = 2 
      call ftmahd(unit, ihdu, hdutype, status)

      ! Check correct hdutype.
      if(hdutype /= 2) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='Scales not stored in binary table')
      end if

      ! Read header NAXIS2 and check same as n_a.
      call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
      if(naxis /= n_scale) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='Inconsistent number of scales')
      end if

      ! Read dilation values from binary table.
      frow=1
      felem=1
      nelem=n_scale
      nullval = -999      ! Arbitrary since will stop and return error 
                          ! if null values detected.
      colnum=1      
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           scales_field(:,1),anynull,status)
      colnum = 2
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           scales_field(:,2),anynull,status)
      if(anynull) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='Null scale values contained in file')
      end if

      ! Don't need to be saved anywhere since should be same as values stored
      ! in filter object.  Just check they are indeed the same.

      ! Get scales from filter.
      call s2fil_filter_get_scale(field%filter, scales_filter)

      ! Check scales consistent.
      if( sum(abs(scales_field - scales_filter)) > TOL ) then
         call s2fil_error(S2FIL_ERROR_FIELD_FILE_INVALID, &
           's2fil_field_io_fits_read', &
           comment_add='Field and filter scales inconsistent')
      end if


      ! --------------------------------------
      ! Tidy up
      ! --------------------------------------

      ! Set status to initialised.
      field%computed_status = .true.
      field%init = .true.

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2fil_field_io_fits_error_check(status, .true.)

      ! Free temporary memory used.
      deallocate(scales_field, scales_filter)

    end subroutine s2fil_field_io_fits_read


    !--------------------------------------------------------------------------
    ! s2fil_field_io_fits_error_check
    !
    !! Check if a fits error has occured and print error message.  Halt
    !! program execution if halt flag is set.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - status: Fits integer status code.
    !!   - halt: Logical to indicate whether to halt program execution if an 
    !!     error is detected.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_fits_error_check(status, halt)

      integer, intent(inout) :: status
      logical, intent(in) :: halt

      character(len=30) :: errtext
      character(len=80) :: errmessage

      !  Check if status is OK (no error); if so, simply return.
      if (status .le. 0) return

      ! The FTGERR subroutine returns a descriptive 30-character text 
      ! string that corresponds to the integer error status number.  
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      ! The FTGMSG subroutine retrieves the oldest message from
      ! the stack and shifts any remaining messages on the stack down one
      ! position.  FTGMSG is called repeatedly until a blank message is
      ! returned, which indicates that the stack is empty.  Each error message
      ! may be up to 80 characters in length. 
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write(*,*) trim(errmessage)
          call ftgmsg(errmessage)
      end do

      if(halt) stop

    end subroutine s2fil_field_io_fits_error_check


    !--------------------------------------------------------------------------
    ! s2fil_field_io_fits_exists
    !
    !! Check if a fits file exists.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to check existence of.
    !!   - status: Fits integer status code.
    !!   - exists: Logical indicating whether the fits file already exists.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_fits_exists(filename, status, exists)

      character(len=*), intent(in) :: filename
      integer, intent(inout) :: status
      logical, intent(out) :: exists

      integer :: unit, blocksize
      logical :: halt

      ! Simply return if status is already greater than zero.
      if (status .gt. 0) return

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      call ftopen(unit, filename, 1, blocksize, status)

      ! Check status of opening file.
      if(status == 0) then

        ! File was opened.  Close it and set exists flag accordingly.
        call ftclos(unit, status)
        exists = .true.

      else if (status == 104) then
        
        ! File does not exist.  Reset status and set exists flag accordingly.
         status = 0
         exists = .false.

      else

        ! Some other error occured while opening file.
        halt = .false.
        call s2fil_field_io_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine s2fil_field_io_fits_exists


    !--------------------------------------------------------------------------
    ! s2fil_field_io_fits_del
    !
    !! Delete a fits file.
    !!
    !! Notes:
    !!   - Copied form cswt_tr_mod module.
    !!
    !! Variables:
    !!   - filename: Name of fits file to detele.
    !!   - status: Fits integer status code.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_fits_del(filename, status)

      character(len=*), intent(in) :: filename
      integer, intent(inout) ::  status

      integer :: unit, blocksize

      ! Simply return if status is greater than zero.
      if (status .gt. 0)return

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Try to open the file, to see if it exists.
      call ftopen(unit,filename,1,blocksize,status)

      if(status .eq. 0) then
         ! File was opened;  so now delete it.
         call ftdelt(unit,status)
      else if(status .eq. 103) then
         ! File doesn't exist, so just reset status to zero and clear errors.
          status=0
          call ftcmsg
      else
         ! There was some other error opening the file; delete the file anyway.
         status=0
         call ftcmsg
         call ftdelt(unit,status)
      end if

      ! Free the unit number for later reuse.
      call ftfiou(unit, status)

    end subroutine s2fil_field_io_fits_del


    !--------------------------------------------------------------------------
    ! s2fil_field_io_txt_dots_write
    !
    !! Write source positions to an output dot text file.
    !!
    !! Variables:
    !!   - filename: Name of text file to write source positions to.
    !!   - alpha: Alpha (phi) spherical coordinate array of source positions 
    !!     to write.
    !!   - beta: Beta (theta) spherical coordinate array of source positions 
    !!     to write.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_txt_dots_write(filename, alpha, beta)

      character(S2_STRING_LEN), intent(in) :: filename
      real(s2_sp), intent(in) :: alpha(:)
      real(s2_sp), intent(in) :: beta(:)

      character(len=1), parameter :: COMMENT_CHAR = '#'
      integer :: fileid = 13, i_source

      ! Check size of alpha and beta arrays are the same.
      if(size(alpha) /= size(beta)) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, 's2fil_field_io_txt_dots_write', &
           comment_add='Alpha and beta arrays of inconsistent size')
      end if

      open(unit=fileid, file=filename, status='replace', action='write')

      write(fileid,'(a,a)') COMMENT_CHAR, ' Localised source positions'
      write(fileid,'(a,a)') COMMENT_CHAR, ' Written by s2fil-0.1'

      write(fileid,'(a,i5)') 'n_sources: ', size(alpha)

      do i_source = 1,size(alpha)
         write(fileid,'(f7.5,f9.5)') alpha(i_source), beta(i_source)      
      end do

      close(fileid)

    end subroutine s2fil_field_io_txt_dots_write


    !--------------------------------------------------------------------------
    ! s2fil_field_io_txt_dots_read
    !
    !! Read source positions from an input dot text file.
    !!
    !! Notes:
    !!   - Alpha and beta arrays are allocated herein and must be freed by 
    !!     the calling routine.
    !!
    !! Variables:
    !!   - filename: Name of text file to read source positions from.
    !!   - alpha: Alpha (phi) spherical coordinate array of source positions 
    !!     read.
    !!   - beta: Beta (theta) spherical coordinate array of source positions 
    !!     read.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_io_txt_dots_read(filename, alpha, beta)

      character(S2_STRING_LEN), intent(in) :: filename
      real(s2_sp), intent(out), allocatable :: alpha(:), beta(:)

      character(len=1), parameter :: COMMENT_CHAR = '#'
      integer :: fileid = 14, n_source = 0, i_source, fail
      character(len=S2_STRING_LEN) :: line, line2

      open(fileid, file=filename, form='formatted', status='old')

      ! Ignore leading comment lines.
      line = COMMENT_CHAR
      do while(line(1:1) == COMMENT_CHAR)
         read(fileid,'(a)') line
      end do

      ! Read number of source positions from last line read (that is actually
      ! not a comment line).
      read(line, *) line2, n_source

      ! Allocate space for alpha and beta arrays.
      allocate(alpha(1:n_source), stat=fail)
      allocate(beta(1:n_source), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, 's2fil_field_io_txt_dots_read')
      end if

      ! Read source positions.
      do i_source = 1,n_source
         read(fileid,*) alpha(i_source), beta(i_source)
      end do
         
      close(fileid)

    end subroutine s2fil_field_io_txt_dots_read


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2fil_field_get_init
    !
    !! Get init variable from the passed field.
    !!
    !! Variables:
    !!   - field: Field object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_init(field) result(init)
      
      type(s2fil_field), intent(in) :: field
      logical :: init

      init = field%init

    end function s2fil_field_get_init

 
    !--------------------------------------------------------------------------
    ! s2fil_field_get_n_gamma
    !
    !! Get n_gamma variable from the passed field.
    !!
    !! Variables:
    !!   - field: Field object to get the variable of.
    !!   - n_gamma: Object n_gamma variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_n_gamma(field) result(n_gamma)
      
      type(s2fil_field), intent(in) :: field
      integer :: n_gamma

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_gamma')
      end if

      n_gamma = field%n_gamma

    end function s2fil_field_get_n_gamma


    !--------------------------------------------------------------------------
    ! s2fil_field_get_computed_status
    !
    !! Get computed_status variable from the passed field.
    !!
    !! Variables:
    !!   - field: Field object to get the variable of.
    !!   - computed_status: Object computed_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_computed_status(field) result(computed_status)
      
      type(s2fil_field), intent(in) :: field
      logical :: computed_status

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_computed_status')
      end if

      computed_status = field%computed_status

    end function s2fil_field_get_computed_status


    !--------------------------------------------------------------------------
    ! Get hidden variable routines
    ! (Variables stored in child data structures)
    !--------------------------------------------------------------------------
   
    !--------------------------------------------------------------------------
    ! s2fil_field_get_n_alpha
    !
    !! Get n_alpha variable from the passed field.  Note n_alpha stored in 
    !! tr structure.
    !!
    !! Variables:
    !!   - field: Field object to get the n_alpha of.
    !!   - n_alpha: N_alpha variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_n_alpha(field) result(n_alpha)
      
      type(s2fil_field), intent(in) :: field
      integer :: n_alpha

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_alpha')
      end if

      ! Check tr structure computed.
      if(.not. field%computed_status) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_alpha', &
          comment_add='Tr structure not initialised')
     end if

      n_alpha = cswt_tr_get_n_alpha(field%tr(1))

    end function s2fil_field_get_n_alpha


    !--------------------------------------------------------------------------
    ! s2fil_field_get_n_beta
    !
    !! Get n_beta variable from the passed field.  Note n_beta stored in 
    !! tr structure.
    !!
    !! Variables:
    !!   - field: Field object to get the n_beta of.
    !!   - n_beta: N_beta variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_n_beta(field) result(n_beta)
      
      type(s2fil_field), intent(in) :: field
      integer :: n_beta

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_beta')
      end if

      ! Check tr structure computed.
      if(.not. field%computed_status) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_beta', &
          comment_add='Tr structure not initialised')
     end if

      n_beta = cswt_tr_get_n_beta(field%tr(1))

    end function s2fil_field_get_n_beta


    !--------------------------------------------------------------------------
    ! s2fil_field_get_lmax
    !
    !! Get lmax variable from the passed field.  Note lmax stored in 
    !! filter object.
    !!
    !! Variables:
    !!   - field: Field object to get the lmax of.
    !!   - lmax: Lmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_lmax(field) result(lmax)
      
      type(s2fil_field), intent(in) :: field
      integer :: lmax

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_lmax')
      end if

      lmax = s2fil_filter_get_lmax(field%filter)

    end function s2fil_field_get_lmax


    !--------------------------------------------------------------------------
    ! s2fil_field_get_mmax
    !
    !! Get mmax variable from the passed field.  Note mmax stored in 
    !! filter object.
    !!
    !! Variables:
    !!   - field: Field object to get the mmax of.
    !!   - mmax: Mmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_mmax(field) result(mmax)
      
      type(s2fil_field), intent(in) :: field
      integer :: mmax

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_mmax')
      end if

      mmax = s2fil_filter_get_mmax(field%filter)

    end function s2fil_field_get_mmax


    !--------------------------------------------------------------------------
    ! s2fil_field_get_n_scale
    !
    !! Get n_scale variable from the passed field.  Note nside stored in 
    !! filter object.
    !!
    !! Variables:
    !!   - field: Field object to get the n_scale of.
    !!   - n_scale: N_scale variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_n_scale(field) result(n_scale)
      
      type(s2fil_field), intent(in) :: field
      integer :: n_scale

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_n_scale')
      end if

      n_scale = s2fil_filter_get_n_scale(field%filter)

    end function s2fil_field_get_n_scale


    !--------------------------------------------------------------------------
    ! s2fil_field_get_nside
    !
    !! Get nside variable from the passed field.  Note nside stored in 
    !! filter object (where it is in turn stored in the tmpl sky)
    !!
    !! Variables:
    !!   - field: Field object to get the nside of.
    !!   - nside: Nside variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_field_get_nside(field) result(nside)
      
      type(s2fil_field), intent(in) :: field
      integer :: nside

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_nside')
      end if

      nside = s2fil_filter_get_nside(field%filter)

    end function s2fil_field_get_nside


    !--------------------------------------------------------------------------
    ! s2fil_field_get_scale_val
    !
    !! Get value from filter scale array corresponding to the index iscale. 
    !! Note value copied is 2D array (since 2D dilation/scale).
    !!
    !! Notes:
    !!   - Space for one_scale must be allocated and deallocated by calling
    !!     routine.
    !!
    !! Variables:
    !!   - field: Field containing filter objects to get the variable of.
    !!   - iscale: Index of scale to get.
    !!   - one_scale: The single scale (2D for 2D dilation) returned 
    !!     corresponding to the specified scale.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2006
    !
    ! Revisions:
    !   April 2006 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_field_get_scale_val(field, iscale, one_scale)
      
      type(s2fil_field), intent(in) :: field
      integer, intent(in) :: iscale
      real(s2_sp), intent(out) :: one_scale(:)

      ! Check object initialised.
      if(.not. field%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_field_get_scale_val')
      end if

      ! Check sizes consistent.
      if(size(one_scale) /= S2FIL_SCALE_DIM2_SIZE) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
           's2fil_field_get_scale_val')
      end if

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > s2fil_filter_get_n_scale(field%filter)) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_field_get_scale_val', comment_add='iscale outside range')
      end if

      call s2fil_filter_get_scale(field%filter, iscale, one_scale)

    end subroutine s2fil_field_get_scale_val


end module s2fil_field_mod
