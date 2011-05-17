!------------------------------------------------------------------------------
! s2fil_filter_mod -- S2FIL library filter class
!
!! Functionality to compute (and store) optimal filters from a background 
!! noise process and template function defined on the sky.  Both spherical 
!! directional matched filters and scale adaptive filters may be constructed.
!! Filters at various scales may be constrcuted either by dilating the 
!! original template (`correct' approach) or by dilating the first computed
!! optimal filter.
!!
!! Notes:
!!   - Regarding beams: See comment in s2fil_filter_init_data comment header.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 January 2005
!
! Revisions:
!   January 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module s2fil_filter_mod

  use s2_types_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2fil_types_mod
  use s2fil_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

   public :: &
    s2fil_filter_init, &
    s2fil_filter_free, &
    s2fil_filter_io_fits_write, &
    s2fil_filter_write_filter_map, &
    s2fil_filter_get_init, &
    s2fil_filter_get_lmax, &
    s2fil_filter_get_mmax, &
    s2fil_filter_get_nside, &
    s2fil_filter_get_n_scale, &
    s2fil_filter_get_norm_pres_dil, &
    s2fil_filter_get_filter_type, &
    s2fil_filter_get_scale_type, &
    s2fil_filter_get_background, &
    s2fil_filter_get_beam, &
    s2fil_filter_get_tmpl, &
    s2fil_filter_get_scale, &
    s2fil_filter_get_sky, &
    s2fil_filter_get_sky_status, &
    s2fil_filter_get_filter_cl, &
    s2fil_filter_get_tmpl_cl


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface s2fil_filter_init
     module procedure &
       s2fil_filter_init_data, &
       s2fil_filter_init_copy, &
       s2fil_filter_init_file
  end interface

  interface s2fil_filter_get_scale
     module procedure &
       s2fil_filter_get_scale_array, &
       s2fil_filter_get_scale_val
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Filter type: Matched filter.
  integer, public, parameter :: S2FIL_FILTER_TYPE_MF = 1

  !! Filter type: Scale adaptive filter.
  integer, public, parameter :: S2FIL_FILTER_TYPE_SAF = 2

  !! Matched filter string description.
  character(len=*), public, parameter :: S2FIL_FILTER_TYPE_STR_MF = 'mf'

  !! Scale adaptive filter string description.
  character(len=*), public, parameter :: S2FIL_FILTER_TYPE_STR_SAF = 'saf'

  !! Scale type: Scale template before 
  integer, public, parameter :: S2FIL_FILTER_SCALE_TYPE_TMPL = 1

  !! Scale type: Scale filter
  integer, public, parameter :: S2FIL_FILTER_SCALE_TYPE_FILTER = 2

  !! Scale type template string description.
  character(len=*), public, parameter :: &
    S2FIL_FILTER_SCALE_TYPE_STR_TMPL = 'tmpl'

  !! Scale type filter string description.
  character(len=*), public, parameter :: &
    S2FIL_FILTER_SCALE_TYPE_STR_FILTER = 'filter'


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: s2fil_filter
     private
     logical :: init = .false.
     integer :: lmax = 0
     integer :: mmax = 0
     integer :: n_scale = 0
     logical :: norm_pres_dil = S2FIL_NORM_PRES_DIL_DEFAULT
     integer :: filter_type = S2FIL_FILTER_TYPE_MF
     integer :: scale_type = S2FIL_FILTER_SCALE_TYPE_TMPL
     type(s2_pl) :: background
     type(s2_pl) :: beam
     type(s2_sky) :: tmpl
     real(s2_sp), allocatable :: scale(:,:)
     type(s2_sky), allocatable :: sky(:)
     logical :: sky_status = .false.
     logical :: beam_status = .false.
  end type s2fil_filter

 
  !----------------------------------------------------------------------------
  
  contains

  
    !--------------------------------------------------------------------------
    ! s2fil_filter_init_data
    !
    !! Initialise a filter structure directly from data.
    !!
    !! Notes:
    !!   - Tmpl is just a function defined on the sphere and does not
    !!     have an inherent scale(/dilation) or rotation. 
    !!     The _template function is assumed to be of scale one and_
    !!     _zero rotation_.  Filters for all other scales are constructed 
    !!     from the template and its assumed unity scale.
    !!   - Filter lmax is set as lmax of background.
    !!   - If mmax is present and less than lmax, then set, else mmax is
    !!     set to lmax.
    !!   - Regarding beams: The filter template attribute is saved for a 
    !!     dilation of unity.  When used the template is scaled to the 
    !!     appropriate size.  When a beam is present the beam must be applied 
    !!     each time the template is used once the scaling has been performed.
    !!     The beam must already be applied to the background before it is 
    !!     used to initialise the filter object.  That is, the filter 
    !!     background attribute saved has already had a beam applied.  This 
    !!     is because the background may consist of a cmb (for which the beam
    !!     must be applied) plus noise (for which the beam is *not* applied 
    !!     to).
    !!
    !! Variables:
    !!  - tmpl: Template compact object defined on the sky.
    !!  - background: Background power spectrum template object is embedded in.
    !!  - scales: Scales of tmpl (dilations) to construct filters for.
    !!  - [mmax]: Optional mmax which is only used if it's below the 
    !!    background lmax.  Otherwise background lmax is used.
    !!  - [norm_pres_dil]: Optional to specify whether norm preserving 
    !!    dilations are to be performed or not.  Over-rides default 
    !!    filter%norm_pres_dil if included.
    !!  - filter: Optimal filter constructed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_init_data(tmpl, background, scale, filter_type, &
      scale_type, mmax, norm_pres_dil, beam, heuristic) result(filter)

      type(s2_sky), intent(in) :: tmpl
      type(s2_pl), intent(in) :: background
      real(s2_sp), intent(in) :: scale(:,:)
      integer, intent(in) :: filter_type
      integer, intent(in) :: scale_type
      integer, intent(in), optional :: mmax
      logical, intent(in), optional :: norm_pres_dil
      type(s2_pl), intent(in), optional :: beam
      logical, intent(in), optional :: heuristic
      type(s2fil_filter) :: filter

      real(s2_sp), parameter :: SCALE_TOL = 1e-5
      integer :: fail

      ! Check object not already initialised.
      if(filter%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_init_data')
      end if

      ! Check scale vector is correct size for list of 2D dilations.
      if(size(scale,2) /= S2FIL_SCALE_DIM2_SIZE) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
           's2fil_filter_init_data')
      end if

      ! Set number of scales.
      filter%n_scale = size(scale,1)
      if(filter%n_scale <= 0) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
           's2fil_filter_init_data')
      end if 

      ! Check filter_type valid.
      if(filter_type /= S2FIL_FILTER_TYPE_MF &
        .and. filter_type /= S2FIL_FILTER_TYPE_SAF) then
         call s2fil_error(S2FIL_ERROR_FILTER_TYPE_INVALID, &
           's2fil_filter_init_data')
      end if

      ! Check scale_type valid.
      if(scale_type /= S2FIL_FILTER_SCALE_TYPE_TMPL &
        .and. scale_type /= S2FIL_FILTER_SCALE_TYPE_FILTER) then
         call s2fil_error(S2FIL_ERROR_SCALE_TYPE_INVALID, &
           's2fil_filter_init_data')
      end if

      ! Set filter beam status.
      if(present(beam)) filter%beam_status = .true.

      ! Check background and beam have same lmax.
      if(filter%beam_status) then
         if(s2_pl_get_lmax(background) /= s2_pl_get_lmax(beam)) then
            call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
              's2fil_filter_init_data', &
              comment_add='Beam and background have inconsistent size')
         end if
      end if

      ! Initialise filter variables.
      filter%filter_type = filter_type
      filter%scale_type = scale_type
      if(present(norm_pres_dil)) filter%norm_pres_dil = norm_pres_dil
      filter%lmax = s2_pl_get_lmax(background)
      if(present(mmax)) then
         if(mmax < filter%lmax) then
            filter%mmax = mmax
         else
            filter%mmax = filter%lmax
         end if
      else
         filter%mmax = filter%lmax
      end if
      filter%background = s2_pl_init(background)
      if(filter%beam_status) filter%beam = s2_pl_init(beam)
      allocate(filter%scale(filter%n_scale,S2FIL_SCALE_DIM2_SIZE), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, 's2fil_filter_init_data')
      end if
      filter%scale = scale
      filter%tmpl = s2_sky_init(tmpl)      

      ! Set lmax and mmax within tmpl sky also so when go to compute
      ! alms correct sizes set.
      call s2_sky_set_lmax(filter%tmpl, filter%lmax, filter%mmax)

      ! Set filter init status.
      ! Note this must be set before call s2fil_filter_comp_filter.
      filter%sky_status = .false.
      filter%init = .true.

      ! Compute filters for each scales.
      call s2fil_filter_comp_filter(filter, heuristic)

      ! Now set filter status.
      filter%sky_status = .true.

    end function s2fil_filter_init_data


    !--------------------------------------------------------------------------
    ! s2fil_filter_init_copy
    !
    !! Initialise a filter structure as a copy of another filter.
    !!
    !! Variables:
    !!  - orig: The original filter to copy.
    !!  - copy: The initialised filter copied from the original.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_init_copy(orig) result(copy)

      type(s2fil_filter), intent(in) :: orig
      type(s2fil_filter) :: copy

      integer :: iscale, fail

      ! Check original object initialised.
      if(.not. orig%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_init_copy')
      end if

      ! Copy object attributes.

      copy%tmpl = s2_sky_init(orig%tmpl)
      copy%background = s2_pl_init(orig%background)
      copy%filter_type = orig%filter_type
      copy%scale_type = orig%scale_type
      copy%n_scale = orig%n_scale
      copy%lmax = orig%lmax
      copy%mmax = orig%mmax
      copy%norm_pres_dil = orig%norm_pres_dil
      copy%beam_status = orig%beam_status
      if(copy%beam_status) copy%beam = s2_pl_init(orig%beam)

      ! Allocate space for and copy scales.
      allocate(copy%scale(copy%n_scale,S2FIL_SCALE_DIM2_SIZE), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, 's2fil_filter_init_copy')
      end if
      copy%scale = orig%scale

      copy%sky_status = orig%sky_status

      if(copy%sky_status) then

         ! Allocate space for and copy filters for each scale.

         allocate(copy%sky(copy%n_scale), stat=fail)
         if(fail /= 0) then
            call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
                 's2fil_filter_init_copy', &
                 comment_add='Failed whilst allocating space for filters')
         end if

         do iscale = 1,copy%n_scale
            copy%sky(iscale) = s2_sky_init(orig%sky(iscale))
         end do

      end if

      copy%init = .true.

    end function s2fil_filter_init_copy


    !--------------------------------------------------------------------------
    ! s2fil_filter_init_file
    !
    !! Initialise a filter structure from a filter file.
    !!
    !! Variables:
    !!  - filename: The fits file to read the filter from.
    !!  - filter: Optimal filter constructed with the data read
    !!    from the file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_init_file(filename) result(filter)

      character(len=*), intent(in) :: filename
      type(s2fil_filter) :: filter

      ! Check object not already initialised.
      if(filter%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_init_data')
      end if

      ! Read filter file.
      call s2fil_filter_io_fits_read(filename, filter)

      ! All status flags set in s2fil_filter_io_fits_read routine.

    end function s2fil_filter_init_file


    !--------------------------------------------------------------------------
    ! s2fil_filter_comp_filter
    !
    !! Compute the filter structure filter variable. Space for the filters 
    !! is allocated here.
    !!
    !! Notes:
    !!   - All other filter variables must already be allocated and defined.
    !!   - Filter status should not be set before calling functions, but 
    !!     sould be set be calling routine following call to this function.
    !!   - Filter nside (for each scale) set from tmpl nside.
    !!   - Start from *l=2*, thus monopole and dipole not considered.  
    !!     This is usually the case of interest.
    !!
    !! Variables:
    !!  - filter: Optimal filter object that filter variable is 
    !!    computed for.
    !!  - [heu]: Logical specifying whether to adopt heuristic lmax level.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_comp_filter(filter, heu) 

      type(s2fil_filter), intent(inout) :: filter
      logical, intent(in), optional :: heu

      real(s2_sp), parameter :: ZERO_TOL = 1.0e-10

      complex(s2_spc), allocatable :: filter_alm(:,:)
      complex(s2_spc), allocatable :: tmpl_alm(:,:)
      real(s2_sp), allocatable :: bk_pl(:)
      type(s2_sky) :: temp_tmpl, temp_filter_sky
      integer :: l, m, iscale, fail, nside    

      real(s2_sp) :: dil1, dil2
      complex(s2_spc) :: b, F
      real(s2_sp) :: a, c, delta

      integer :: tmpl_lmax
!      real(s2_sp), parameter :: CUTOFF_THRES = 0.9999e0
      real(s2_sp) :: CUTOFF_THRES = 0.9999e0
      logical :: heuristic = .true.
      real(s2_sp), allocatable :: cl_tmpl_vals(:)
      type(s2_pl) :: cl_tmpl
      real(s2_sp) :: cutoff

      ! Set heuristic technique status.
      if(present(heu)) heuristic = heu

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_comp_filter')
      end if 

      ! Check sky_status not set.
      if(filter%sky_status) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_comp_filter', &
          comment_add='Attempting to recompute filter')
      end if 

      ! Allocate space for filter for each scale.
      allocate(filter%sky(filter%n_scale), stat=fail)
      if(fail /= 0) then               
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
           's2fil_filter_comp_filter', &
           comment_add='Failed whilst allocating space for filters')
      end if

      ! Allocate space for temporary arrays used when calculating
      ! optimal filter (tmpl_alm, bk_pl and filter_alm).
      allocate(tmpl_alm(0:filter%lmax, 0:filter%mmax), stat=fail)
      allocate(bk_pl(0:filter%lmax), stat=fail)
      allocate(filter_alm(0:filter%lmax, 0:filter%mmax), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
           's2fil_filter_comp_filter')
      end if

      ! Get background spectrum (same for each scale so can do 
      ! outside loop).
      call s2_pl_get_spec(filter%background, bk_pl)
      
      ! Compute optimal filter for each scale.
      do iscale = 1,filter%n_scale

         ! Get temporary template to dilate without altering original 
         ! template.
         temp_tmpl = s2_sky_init(filter%tmpl)

         ! Dilate by scale corresponding to iscale.
         call s2_sky_dilate(temp_tmpl, filter%scale(iscale,1), &
           filter%scale(iscale,2), filter%norm_pres_dil)

         ! Compute alms of dilated template.
         call s2_sky_compute_alm(temp_tmpl, filter%lmax, filter%mmax)

         ! Convolve template at current scale with beam if present.
         if(filter%beam_status) then
            call s2_sky_conv(temp_tmpl, filter%beam)
         end if
      
         ! Get dilated template alms.
         call s2_sky_get_alm(temp_tmpl, tmpl_alm)

         ! Reset filter alms.
         filter_alm = cmplx(0.0e0,0.0e0)


         ! --------------------------------------
         ! Find band limit of template.
         ! I.e. min l that for which template has 
         ! no power above.
         ! --------------------------------------

         ! Compute template cl and extract spectrum values.
         cl_tmpl = s2_sky_get_cl(temp_tmpl)
         allocate(cl_tmpl_vals(0:s2_pl_get_lmax(cl_tmpl)), stat=fail)
         if(fail /= 0) then
            call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
              's2fil_filter_comp_filter')
         end if
         call s2_pl_get_spec(cl_tmpl, cl_tmpl_vals)
         
         ! Replace tmpl cl vals with cummulative sum.
         do l = 1,s2_pl_get_lmax(cl_tmpl)  ! Note starts from 1.
            cl_tmpl_vals(l) = (2*l+1) *cl_tmpl_vals(l) + cl_tmpl_vals(l-1) 
         end do

         ! Set cutoff.
         cutoff = CUTOFF_THRES * cl_tmpl_vals(s2_pl_get_lmax(cl_tmpl))
        
         ! Step through ls until find first value in cummulative array that
         ! is above the cutoff.
         do l = 0,s2_pl_get_lmax(cl_tmpl)  
            tmpl_lmax = l
            if(cl_tmpl_vals(l) >= cutoff) exit   ! Stop loop, found tmpl_lmax.
         end do

         deallocate(cl_tmpl_vals)
         call s2_pl_free(cl_tmpl)

         ! Set lmax to consider up to depending on whether heuristic approach
         ! selected or not.
         if(heuristic) then
            tmpl_lmax = min(floor(tmpl_lmax * 1.5e0), filter%lmax)
            write(*,'(a,a,i5)') 's2fil_filter_compute_filter> ', &
                 'Considering up to heuristically determined lmax =', &
                 tmpl_lmax      
         else
            tmpl_lmax = filter%lmax
         end if


         ! --------------------------------------
         ! Compute filter variables required.
         ! --------------------------------------

         ! Compute 'a' variable (required for both MF and SAF).
         a = 0.0e0
         do l = 2,tmpl_lmax  !filter%lmax

            ! Check won't divide by zero.
            if(abs(bk_pl(l)) < ZERO_TOL) then 
               call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                  's2fil_filter_comp_filter', &
                  comment_add='Background pl spectrum close to zero')
               stop
            end if

            a = a + abs(tmpl_alm(l,0))**2.0e0 / bk_pl(l)           
            do m = 1,min(l,filter%mmax)
               a = a + 2*abs(tmpl_alm(l,m))**2.0e0 / bk_pl(l)           
            end do
         end do

         ! Compute b, c and delta if computing SAF (not required for MF).
         if(filter%filter_type == S2FIL_FILTER_TYPE_SAF) then

            ! Compute SAF for dilation with cocycle.             
            ! (p=2)
            if(filter%norm_pres_dil) then
            
               ! Compute 'b' SAF filter variable.
               b = cmplx(0.0e0,0.0e0)
               do l = 2,tmpl_lmax  !filter%lmax

                  ! Check won't divide by zero.
                  if(abs(bk_pl(l)) < ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add='Background pl spectrum close to zero')
                     stop
                  end if

                  m = 0
                  F = l * conjg(tmpl_alm(l,m)) &
                       - sqrt(real(l**2 - m**2,s2_sp)) * conjg(tmpl_alm(l-1,m))
                  b = b + F * tmpl_alm(l,m) / bk_pl(l)
                  
                  do m = 1,min(l,filter%mmax)
                     F = l * conjg(tmpl_alm(l,m)) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * conjg(tmpl_alm(l-1,m))
                     b = b + 2 * F * tmpl_alm(l,m) / bk_pl(l)
                  end do
               end do
               
               ! Compute 'c' SAF filter variable.
               c = cmplx(0.0e0,0.0e0)
               do l = 2,tmpl_lmax  !filter%lmax
                  
                  ! Check won't divide by zero.
                  if(abs(bk_pl(l)) < ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add='Background pl spectrum close to zero')
                     stop
                  end if
                  
                  m = 0
                  F = l * tmpl_alm(l,m) &
                       - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                  c = c + abs(F)**2.0e0 /  bk_pl(l)
                  
                  do m = 1,min(l,filter%mmax)
                     F = l * tmpl_alm(l,m) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                     c = c + 2 * abs(F)**2.0e0 /  bk_pl(l)
                  end do
               end do
               
               ! Compute 'delta' SAF filter variable.
               delta = a * c - abs(b)**2.0e0
               
            ! Compute SAF for dilation without cocyle.
            ! (p=infinity)            
            else

               ! Compute 'b' SAF filter variable.
               b = cmplx(0.0e0,0.0e0)
               do l = 2,tmpl_lmax  !filter%lmax

                  ! Check won't divide by zero.
                  if(abs(bk_pl(l)) < ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add='Background pl spectrum close to zero')
                     stop
                  end if

                  m = 0
                  F = (l-1) * conjg(tmpl_alm(l,m)) &
                       - sqrt(real(l**2 - m**2,s2_sp)) * conjg(tmpl_alm(l-1,m))
                  b = b + F * tmpl_alm(l,m) / bk_pl(l)
                  
                  do m = 1,min(l,filter%mmax)
                     F = (l-1) * conjg(tmpl_alm(l,m)) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * conjg(tmpl_alm(l-1,m))
                     b = b + 2 * F * tmpl_alm(l,m) / bk_pl(l)
                  end do
               end do
               
               ! Compute 'c' SAF filter variable.
               c = cmplx(0.0e0,0.0e0)
               do l = 2,tmpl_lmax  !filter%lmax
                  
                  ! Check won't divide by zero.
                  if(abs(bk_pl(l)) < ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add='Background pl spectrum close to zero')
                     stop
                  end if
                  
                  m = 0
                  F = (l-1) * tmpl_alm(l,m) &
                       - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                  c = c + abs(F)**2.0e0 /  bk_pl(l)
                  
                  do m = 1,min(l,filter%mmax)
                     F = (l-1) * tmpl_alm(l,m) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                     c = c + 2 * abs(F)**2.0e0 /  bk_pl(l)
                  end do
               end do
               
               ! Compute 'delta' SAF filter variable.
               delta = a * c - abs(b)**2.0e0

            end if

         end if


         ! --------------------------------------
         ! Compute filter alms.
         ! --------------------------------------

         if(filter%filter_type == S2FIL_FILTER_TYPE_MF) then

            filter_alm = cmplx(0.0e0, 0.0e0)

            do l = 2,tmpl_lmax  !filter%lmax

               ! Check won't divide by zero.
               if(abs(bk_pl(l)) < ZERO_TOL .or. abs(a) < ZERO_TOL) then 
                  call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                       's2fil_filter_comp_filter', &
                       comment_add=&
                       'Background pl spectrum or ''a'' close to zero')
               end if

               do m = 0,min(l,filter%mmax)
                  filter_alm(l,m) = tmpl_alm(l,m) / (a * bk_pl(l))
               end do

            end do
          
            ! Save computed filter for this scale.
            filter%sky(iscale) = s2_sky_init(filter_alm, filter%lmax, &
              filter%mmax, s2fil_filter_get_nside(filter))
            
         else if(filter%filter_type == S2FIL_FILTER_TYPE_SAF) then

            filter_alm = cmplx(0.0e0, 0.0e0)

            if(filter%norm_pres_dil) then
               
               ! Compute SAF for dilation with cocycle.             
               ! (p=2)

               do l = 2,tmpl_lmax  !filter%lmax

                  ! Check won't divide by zero.
                  if(abs(bk_pl(l))<ZERO_TOL .or. abs(delta)<ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add=&
                          'Background pl spectrum or ''delta'' close to zero')
                  end if

                  do m = 0,min(l,filter%mmax)                     
                      F = l * tmpl_alm(l,m) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                      filter_alm(l,m) = (c * tmpl_alm(l,m) - b * F) &
                           / (delta * bk_pl(l))
                  end do

               end do

               ! Save computed filter for this scale.
               filter%sky(iscale) = s2_sky_init(filter_alm, filter%lmax, &
                    filter%mmax, s2fil_filter_get_nside(filter))

            else

               ! Compute SAF for dilation without cocyle.
               ! (p=infinity)

               do l = 2,tmpl_lmax  !filter%lmax
                     
                  ! Check won't divide by zero.
                  if(abs(bk_pl(l))<ZERO_TOL .or. abs(delta)<ZERO_TOL) then 
                     call s2fil_error(S2FIL_ERROR_DIV_BY_ZERO, &
                          's2fil_filter_comp_filter', &
                          comment_add=&
                          'Background pl spectrum or ''delta'' close to zero')
                  end if

                  do m = 0,min(l,filter%mmax)                     
                      F = (l-1) * tmpl_alm(l,m) &
                          - sqrt(real(l**2 - m**2,s2_sp)) * tmpl_alm(l-1,m)
                      filter_alm(l,m) = (c * tmpl_alm(l,m) - b * F) &
                           / (delta * bk_pl(l))
                  end do

               end do

               ! Save computed filter for this scale.
               filter%sky(iscale) = s2_sky_init(filter_alm, filter%lmax, &
                    filter%mmax, s2fil_filter_get_nside(filter))

            end if

         end if

         ! Free temporary template used for this scale.
         call s2_sky_free(temp_tmpl)

         ! If scale filter rather than templates, then stop after first 
         ! filter constructed (i.e. exit do loop).
         if(filter%scale_type == S2FIL_FILTER_SCALE_TYPE_FILTER) exit

      end do


      ! --------------------------------------
      ! Scale filters if required.
      ! --------------------------------------

      ! If scaling filters then loop over dilations will have stopped
      ! after first scale.  Compute other scales here.
      ! If scaling templates, then previous loop will have covered 
      ! all dilation and we're finished.
      if(filter%scale_type == S2FIL_FILTER_SCALE_TYPE_FILTER &
         .and. filter%n_scale > 1) then

         ! Get nside from tmpl.
         nside = s2fil_filter_get_nside(filter)

         ! Compute filter map from alms for first filter.
         ! (Need map to perform dilations.)
         call s2_sky_compute_map(filter%sky(1), nside)
         
         do iscale = 2,filter%n_scale

            ! Create temporary copy of first filter sky that has 
            ! already been computed.          
            temp_filter_sky = s2_sky_init(filter%sky(1))            

            ! Create current filter by dilating first filter in real space.
            ! First filter constructed at first scale, so only need
            ! to dilate from this to current scale.
            dil1 = filter%scale(iscale,1) / filter%scale(1,1)
            dil2 = filter%scale(iscale,2) / filter%scale(1,2)
            call s2_sky_dilate(temp_filter_sky, dil1, dil2, &
              filter%norm_pres_dil)

            ! Compute alms of current dilated filter.
            call s2_sky_compute_alm(temp_filter_sky, filter%lmax, filter%mmax)

            ! Save computed filter for this scale.
            filter%sky(iscale) = s2_sky_init(temp_filter_sky)

            ! Free temp filter sky.
            call s2_sky_free(temp_filter_sky)

         end do

      end if

      ! Free memory used.
      deallocate(filter_alm)
      deallocate(tmpl_alm)
      deallocate(bk_pl)

    end subroutine s2fil_filter_comp_filter


    !--------------------------------------------------------------------------
    ! s2fil_filter_free
    !
    !! Free all data associated with an initialised filter and reset all other
    !! attributes.
    !!
    !! Variables:
    !!  - filter: The filter to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_free(filter)

      type(s2fil_filter), intent(inout) :: filter

      integer :: iscale

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_free')
      end if 

      ! Free space.
      call s2_sky_free(filter%tmpl)
      call s2_pl_free(filter%background)
      if(filter%beam_status) call s2_pl_free(filter%beam)
      deallocate(filter%scale)
      
      ! Free each filter then deallocate filter%sky (if sky_status set).
      if(filter%sky_status) then
         do iscale = 1,filter%n_scale
            call s2_sky_free(filter%sky(iscale))
         end do
         deallocate(filter%sky)
      end if

      ! Reset attributes.
      filter%init = .false.
      filter%sky_status = .false.
      filter%filter_type = S2FIL_FILTER_TYPE_MF
      filter%scale_type = S2FIL_FILTER_SCALE_TYPE_TMPL
      filter%n_scale = 0
      filter%lmax = 0
      filter%mmax = 0
      filter%norm_pres_dil = S2FIL_NORM_PRES_DIL_DEFAULT

    end subroutine s2fil_filter_free


    !--------------------------------------------------------------------------
    ! File IO routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2fil_filter_io_fits_write
    !
    !! Write filter structure to a fits file.
    !!
    !! Variables:
    !!  - filename: Name of fits file to write to.
    !!  - filter: The filter structure containing all data to be written to fits 
    !!    file.
    !!  - [comment]: Comment to append to header of output fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_io_fits_write(filename, filter, comment)

      character(len=*), intent(in) :: filename
      type(s2fil_filter), intent(in) :: filter
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      logical :: simple, extend, file_exists
      integer :: naxis
      integer :: naxes(1)
      integer :: tfields, nrows, varidat
      character(len=32) :: ttype(2), tform(2), tunit(2), extname
      integer :: frow, felem, colnum

      character(len=S2_STRING_LEN) :: filename_tmpl, filename_sky, &
        filename_background, filename_beam, file_key, comment_child_file
      integer :: i_scale

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_io_fits_write')
      end if 

      ! Check filters constructed.
      if(.not. filter%sky_status) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_io_fits_write', &
        comment_add='Attempt to save filter before computed')
      end if

      ! Define FITS parameters.

      bitpix=-32 ! Real single precision.
      status=0   ! Initialse error status to zero.

      ! Check if file already exists.
      call s2fil_filter_io_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_EXISTS, &
              's2fil_filter_io_fits_write')
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
        '  S2fil filter file created by s2fil-0.1',status)
      call ftpcom(unit, &
        '  Primary extension empty',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if

      call ftpkyj(unit,'LMAX', filter%lmax, &
        'max spherical harmonic l considered',status)
      call ftpkyj(unit,'MMAX', filter%mmax, &
        'max spherical harmonic m considered',status)
      call ftpkyj(unit,'NSCALE', filter%n_scale, &
        'number of scales considered',status)
      call ftpkyl(unit,'NORMDIL',filter%norm_pres_dil, &
        'norm preserving dilation status',status)
      call ftpkyl(unit,'BEAMSTS',filter%beam_status, &
        'beam status',status)

      call ftpkyj(unit,'FTYPE', filter%filter_type, &
        'filter type flag',status)
      select case(filter%filter_type) 
         case(S2FIL_FILTER_TYPE_MF)
            call ftpkys(unit,'FTYPESTR', S2FIL_FILTER_TYPE_STR_MF, &
              'string describing filter type',status)            
         case(S2FIL_FILTER_TYPE_SAF)
            call ftpkys(unit,'FTYPESTR', S2FIL_FILTER_TYPE_STR_SAF, &
              'string describing filter type',status)
      end select

      call ftpkyj(unit,'STYPE', filter%scale_type, &
        'scale type flag',status)
      select case(filter%scale_type) 
         case(S2FIL_FILTER_SCALE_TYPE_TMPL)
            call ftpkys(unit,'STYPESTR', &
              S2FIL_FILTER_SCALE_TYPE_STR_TMPL, &
              'string describing scale type',status)            
         case(S2FIL_FILTER_SCALE_TYPE_FILTER)
            call ftpkys(unit,'STYPESTR', &
               S2FIL_FILTER_SCALE_TYPE_STR_FILTER, &
              'string describing scale type',status)
      end select

      ! Save all skies (inc tmpl) as separate files and write names
      ! to primary header.

      write(comment_child_file, '(a,a)') &
        '  Part of parent file structure: ', &
        trim(filename)

      ! Save tmpl sky file.
      write(filename_tmpl,'(a,a)') filename(1:len(trim(filename))-4), &
        '_tmpl.sky'
      call s2_sky_io_fits_write(filename_tmpl, filter%tmpl, &
        trim(comment_child_file))
      call ftpkys(unit,'TMPL', trim(filename_tmpl), &
        'name of tmpl file',status)   

      ! Save skies for each scale.
      do i_scale = 1,filter%n_scale

         write(filename_sky,'(a,a,i2.2,a)') &
           filename(1:len(trim(filename))-4), &
           '_sky', i_scale, '.sky'

         call s2_sky_io_fits_write(filename_sky, filter%sky(i_scale), &
           trim(comment_child_file))

         write(file_key, '(a,i2.2)') 'SKY', i_scale
         call ftpkys(unit, trim(file_key), trim(filename_sky), &
              'name of filter sky file',status)

      end do

      ! Save background pl file.
      write(filename_background,'(a,a)') &
        filename(1:len(trim(filename))-4), '_bkgnd.cl'
      call s2_pl_io_fits_write(filename_background, filter%background, &
        trim(comment_child_file))
      call ftpkys(unit,'BKGND', trim(filename_background), &
        'name of background file',status)  

      ! Save beam pl file.
      if(filter%beam_status) then
         write(filename_beam,'(a,a)') &
           filename(1:len(trim(filename))-4), '_beam.cl'
         call s2_pl_io_fits_write(filename_beam, filter%beam, &
           trim(comment_child_file))
         call ftpkys(unit,'BEAM', trim(filename_beam), &
           'name of beam file',status)  
      end if

      ! Insert binary table extension for scales.
      extname='SCALE'
      ttype(1)='SC1'
      ttype(2)='SC2'
      tform(1)='1E'
      tform(2)='1E'
      tunit(1)='direct (rad)'
      tunit(2)='direct (rad)'
      tfields=2
      nrows=filter%n_scale
      varidat=0
      call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

      ! Write dilations to binary table.
      frow=1
      felem=1
      colnum=1
      call ftpcle(unit,colnum,frow,felem,nrows,filter%scale(:,1),status)
      colnum=2
      call ftpcle(unit,colnum,frow,felem,nrows,filter%scale(:,2),status)

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2fil_filter_io_fits_error_check(status, .true.)

    end subroutine s2fil_filter_io_fits_write


    !--------------------------------------------------------------------------
    ! s2fil_filter_io_fits_read
    ! 
    !! Read a s2fil_filter file and all associated child files (i.e. template,
    !! background and skies).
    !!
    !! Notes:
    !!   - Child template, background and sky files (as described in the 
    !!     filter file header) must be located in the same directory.
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

    subroutine s2fil_filter_io_fits_read(filename, filter)

      character(len=*), intent(in) :: filename
      type(s2fil_filter), intent(out) :: filter

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: ihdu, hdutype, naxis
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: colnum, frow, felem, nelem
      real(s2_sp) :: nullval

      integer :: fail, i_scale
      character(len=S2_STRING_LEN) :: filename_tmpl, filename_bkgnd, &
        filename_sky, filename_beam, file_key

      ! Check object not already initialised.
      if(filter%init) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_init_data')
      end if

      ! Initialse error status to zero.
      status=0

      ! Check if file already exists.
      call s2fil_filter_io_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_filter_io_fits_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)


      ! --------------------------------------
      ! Read primary header variables
      ! --------------------------------------

      call ftgkyj(unit, 'LMAX', filter%lmax, comment, status)
      call ftgkyj(unit, 'MMAX', filter%mmax, comment, status)
      call ftgkyj(unit, 'NSCALE', filter%n_scale, comment, status)
      call ftgkyl(unit, 'NORMDIL', filter%norm_pres_dil, comment, status)
      call ftgkyl(unit, 'BEAMSTS', filter%beam_status, comment, status)
      call ftgkyj(unit, 'FTYPE', filter%filter_type, comment, status)
      call ftgkyj(unit, 'STYPE', filter%scale_type, comment, status)
 
      ! Check filter_type valid.
      if(filter%filter_type /= S2FIL_FILTER_TYPE_MF &
        .and. filter%filter_type /= S2FIL_FILTER_TYPE_SAF) then
         call s2fil_error(S2FIL_ERROR_FILTER_TYPE_INVALID, &
           's2fil_filter_io_fits_read')
      end if

      ! Check scale_type valid.
      if(filter%scale_type /= S2FIL_FILTER_SCALE_TYPE_TMPL &
        .and. filter%scale_type /= S2FIL_FILTER_SCALE_TYPE_FILTER) then
         call s2fil_error(S2FIL_ERROR_SCALE_TYPE_INVALID, &
           's2fil_filter_io_fits_read')
      end if

      ! Check correct number of HDUs in input file.
      hdunum = 2   ! Primary header plus scales.
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
       call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_filter_io_fits_read', &
           comment_add='Invalid number of extensions')
      end if


      ! --------------------------------------
      ! Read sky and pl files
      ! --------------------------------------

      ! Get template filename from header then read template from 
      ! sky file.
      call ftgkys(unit, 'TMPL', filename_tmpl, comment, status)
      filter%tmpl = s2_sky_init(filename_tmpl, S2_SKY_FILE_TYPE_SKY)

      ! Get background filename form header then read background spectrum
      ! from cl file.
      call ftgkys(unit, 'BKGND', filename_bkgnd, comment, status)
      filter%background = s2_pl_init(filename_bkgnd)

      ! Get beam filename form header then read beam spectrum
      ! from cl file.
      if(filter%beam_status) then
         call ftgkys(unit, 'BEAM', filename_beam, comment, status)
         filter%beam = s2_pl_init(filename_beam)
      end if

      ! Allocate space for filter for each scale.
      allocate(filter%sky(filter%n_scale), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
           's2fil_filter_io_fits_read', &
           comment_add='Failed whilst allocating space for filters')
      end if

      ! Read filter sky file for each scale.
      do i_scale = 1, filter%n_scale

         ! Get filter sky filename from filter fits file header.
         write(file_key,'(a,i2.2)') 'SKY', i_scale
         call ftgkys(unit, trim(file_key), filename_sky, comment, status)

         ! Read filter sky from s2_sky file.
         filter%sky(i_scale) = s2_sky_init(filename_sky, &
           S2_SKY_FILE_TYPE_SKY)

      end do


      ! --------------------------------------
      ! Read scales from next extension
      ! --------------------------------------

      ! Allocate space for scales.
      allocate(filter%scale(filter%n_scale, S2FIL_SCALE_DIM2_SIZE), stat=fail)
      if(fail /= 0) then
         call s2fil_error(S2FIL_ERROR_MEM_ALLOC_FAIL, &
            's2fil_filter_io_fits_read')
      end if

      ! Move to second extension (binary table containing scale values).
      ihdu = 2 
      call ftmahd(unit, ihdu, hdutype, status)

      ! Check correct hdutype.
      if(hdutype /= 2) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_filter_io_fits_read', &
           comment_add='Scales not stored in binary table')
      end if

      ! Read header NAXIS2 and check same as n_a.
      call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
      if(naxis/=filter%n_scale) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_filter_io_fits_read', &
           comment_add='Inconsistent number of scales')
      end if

      ! Read dilation values from binary table.
      frow=1
      felem=1
      nelem=filter%n_scale
      nullval = -999      ! Arbitrary since will stop and return error 
                          ! if null values detected.
      colnum=1      
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           filter%scale(:,1),anynull,status)
      colnum = 2
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           filter%scale(:,2),anynull,status)
      if(anynull) then
         call s2fil_error(S2FIL_ERROR_FIL_FILE_INVALID, &
           's2fil_filter_io_fits_read', &
           comment_add='Null scale values contained in file')
      end if


      ! --------------------------------------
      ! Tidy up
      ! --------------------------------------

      ! Set status to initialised.
      filter%sky_status = .true.
      filter%init = .true.

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call s2fil_filter_io_fits_error_check(status, .true.)

    end subroutine s2fil_filter_io_fits_read


    !--------------------------------------------------------------------------
    ! s2fil_filter_io_fits_error_check
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

    subroutine s2fil_filter_io_fits_error_check(status, halt)

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

    end subroutine s2fil_filter_io_fits_error_check


    !--------------------------------------------------------------------------
    ! s2fil_filter_io_fits_exists
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

    subroutine s2fil_filter_io_fits_exists(filename, status, exists)

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
        call s2fil_filter_io_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine s2fil_filter_io_fits_exists


    !--------------------------------------------------------------------------
    ! s2fil_filter_io_fits_del
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

    subroutine s2fil_filter_io_fits_del(filename, status)

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

    end subroutine s2fil_filter_io_fits_del


    !--------------------------------------------------------------------------
    ! s2fil_filter_write_filter_map
    !
    !! Write a fits sky file for the filter at a specified scale.
    !!
    !! Notes:
    !!   - If the map of the filter chosen is not already computed it will 
    !!     be, hence the returned filter has intent inout since filter maps 
    !!     may have been computed.
    !!
    !! Variables:
    !!  - filter: The filter structure containing the filter map field to 
    !!    write to the file.
    !!  - iscale: Scale of the filter (i.e. particular filter) to output.
    !!  - filename: Name of output fits file.
    !!  - [comment]: Optional comment appended to fits file header if present.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   February 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_write_filter_map(filter, iscale, filename, comment)

      type(s2fil_filter), intent(inout) :: filter
      integer, intent(in) :: iscale
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      logical :: map_constructed

     ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_write_filter_map')
      end if 

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > filter%n_scale) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_filter_write_filter_map', &
            comment_add='iscale outside range')
      end if

      ! Check filters constructed.
      if(.not. filter%sky_status) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_write_filter_map', &
        comment_add='Attempt to save filter before computed')
      end if

      ! If map not already computed then compute.
      map_constructed = s2_sky_get_map_status(filter%sky(iscale))
      if(.not. map_constructed) then
         call s2_sky_compute_map(filter%sky(iscale))
      end if

      ! Write map file.
      call s2_sky_write_map_file(filter%sky(iscale), filename, comment)

    end subroutine s2fil_filter_write_filter_map


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! s2fil_filter_get_init
    !
    !! Get init variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_init(filter) result(init)
      
      type(s2fil_filter), intent(in) :: filter
      logical :: init

      init = filter%init

    end function s2fil_filter_get_init


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_lmax
    !
    !! Get lmax variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - lmax: Object lmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_lmax(filter) result(lmax)

      type(s2fil_filter), intent(in) :: filter
      integer :: lmax

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_lmax')
      end if 

      lmax = filter%lmax

    end function s2fil_filter_get_lmax


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_mmax
    !
    !! Get mmax variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - mmax: Object mmax variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_mmax(filter) result(mmax)

      type(s2fil_filter), intent(in) :: filter
      integer :: mmax

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_mmax')
      end if 

      mmax = filter%mmax

    end function s2fil_filter_get_mmax

    
    !--------------------------------------------------------------------------
    ! s2fil_filter_get_nside
    !
    !! Get n_scale variable from the passed filter.  Note nside stored in
    !! tmpl sky.
    !!
    !! Variables:
    !!   - filter: Filter object to get nside of.
    !!   - nside: Nside of the filter (stored in tmpl sky).
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_nside(filter) result(nside)

      type(s2fil_filter), intent(in) :: filter
      integer :: nside

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_nside')
      end if 

      nside = s2_sky_get_nside(filter%tmpl)

    end function s2fil_filter_get_nside


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_n_scale
    !
    !! Get n_scale variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - n_scale: Object n_scale variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_n_scale(filter) result(n_scale)

      type(s2fil_filter), intent(in) :: filter
      integer :: n_scale

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_n_scale')
      end if 

      n_scale = filter%n_scale

    end function s2fil_filter_get_n_scale


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_norm_pres_dil
    !
    !! Get norm_pres_dil variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - norm_pres_dil: Object norm_pres_dil variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_norm_pres_dil(filter) result(norm_pres_dil)

      type(s2fil_filter), intent(in) :: filter
      logical :: norm_pres_dil

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, &
          's2fil_filter_get_norm_pres_dil')
      end if 

      norm_pres_dil = filter%norm_pres_dil

    end function s2fil_filter_get_norm_pres_dil


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_filter_type
    !
    !! Get filter_type variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - filter_type: Object filter_type variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_filter_type(filter) result(filter_type)

      type(s2fil_filter), intent(in) :: filter
      integer :: filter_type

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_filter_type')
      end if 

      filter_type = filter%filter_type

    end function s2fil_filter_get_filter_type


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_scale_type
    !
    !! Get scale_type variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - scale_type: Object scale_type variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 March 2005
    !
    ! Revisions:
    !   March 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_scale_type(filter) result(scale_type)

      type(s2fil_filter), intent(in) :: filter
      integer :: scale_type

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_scale_type')
      end if 

      scale_type = filter%scale_type

    end function s2fil_filter_get_scale_type


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_background
    !
    !! Get background variable from the passed filter.
    !!
    !! Notes:
    !!   - Initialises a new background (s2_pl) object as a copy of the 
    !!     background stored herein.  The returned (copy) background must be
    !!     freed by the calling routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - background: Copy of object background variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_background(filter) result(background)

      type(s2fil_filter), intent(in) :: filter
      type(s2_pl) :: background

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_background')
      end if 

      ! Make a copy for the returned background.
      background = s2_pl_init(filter%background)

    end function s2fil_filter_get_background


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_beam
    !
    !! Get beam variable from the passed filter.
    !!
    !! Notes:
    !!   - Initialises a new beam (s2_pl) object as a copy of the 
    !!     beam stored herein.  The returned (copy) beam must be
    !!     freed by the calling routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - beam: Copy of object beam variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_beam(filter) result(beam)

      type(s2fil_filter), intent(in) :: filter
      type(s2_pl) :: beam

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_beam')
      end if 

      ! Check beam present.
      if(.not. filter%beam_status) then
        call s2fil_error(S2FIL_ERROR_BEAM_NOT_DEF, 's2fil_filter_get_beam')
     end if

      ! Make a copy for the returned beam.
      beam = s2_pl_init(filter%beam)

    end function s2fil_filter_get_beam


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_tmpl
    !
    !! Get tmpl variable from the passed filter.
    !!
    !! Notes:
    !!   - Initialises a new tmpl (s2_sky) object as a copy of the tmpl stored
    !!     herein.  The returned (copy) tmpl must be freed by the calling 
    !!     routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - tmpl: Copy of object tmpl variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_tmpl(filter) result(tmpl)

      type(s2fil_filter), intent(in) :: filter
      type(s2_sky) :: tmpl

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_tmpl')
      end if 

      ! Make a copy for the returned tmpl.
      tmpl = s2_sky_init(filter%tmpl)

    end function s2fil_filter_get_tmpl


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_scale_array
    !
    !! Get copy of scale array from the passed filter.
    !!
    !! Notes:
    !!   - Space for scale array must be allocated and deallocated by calling
    !!     routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - scale: Object scale array returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_get_scale_array(filter, scale)

      type(s2fil_filter), intent(in) :: filter
      real(s2_sp), intent(out) :: scale(:,:)

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_scale_array')
      end if 

      ! Check sizes consistent.
      if(size(scale,1) /= filter%n_scale .or. &
         size(scale,2) /= S2FIL_SCALE_DIM2_SIZE) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
           's2fil_filter_get_scale_array')
      end if

      scale = filter%scale

    end subroutine s2fil_filter_get_scale_array


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_scale_val
    !
    !! Get value from filter scale array corresponding to the index iscale. 
    !! Note value copied is 2D array (since 2D dilation/scale).
    !!
    !! Notes:
    !!   - Space for one_scale must be allocated and deallocated by calling
    !!     routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - iscale: Index of scale to get.
    !!   - one_scale: The single scale (2D for 2D dilation) returned 
    !!     corresponding to the specified scale.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_filter_get_scale_val(filter, iscale, one_scale)

      type(s2fil_filter), intent(in) :: filter
      integer, intent(in) :: iscale
      real(s2_sp), intent(out) :: one_scale(:)

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_scale_val')
      end if 

      ! Check sizes consistent.
      if(size(one_scale) /= S2FIL_SCALE_DIM2_SIZE) then
         call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
           's2fil_filter_get_scale_val')
      end if

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > filter%n_scale) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_filter_get_scale_val', comment_add='iscale outside range')
      end if

      one_scale = filter%scale(iscale,:)
       
    end subroutine s2fil_filter_get_scale_val
       
    
    !--------------------------------------------------------------------------
    ! s2fil_filter_get_sky
    !
    !! Get filter variable(/sky) corresponding to the iscale index from the
    !! passed filter.
    !!
    !! Notes:
    !!   - Initialises a new filter (s2_sky) object as a copy of the 
    !!     filter stored herein for the specified scale.  The returned 
    !!     (copy) filter must be freed by the calling routine.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - iscale: Index of scale to get.
    !!   - filter: Copy of object filter variable returned corresponding 
    !!     to the specified scale.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_sky(filter, iscale) result(sky)

      type(s2fil_filter), intent(in) :: filter
      integer, intent(in) :: iscale
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_filter')
      end if 

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > filter%n_scale) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_filter_get_filter', comment_add='iscale outside range')
      end if

      ! Check filters constructed.
      if(.not. filter%sky_status) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_get_filter', &
        comment_add='Attempt to get filter before computed')
      end if

      ! Make a copy for the returned filter.
      sky = s2_sky_init(filter%sky(iscale))

    end function s2fil_filter_get_sky


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_sky_status
    !
    !! Get sky_status variable from the passed filter.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - sky_status: Object sky_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_sky_status(filter) result(sky_status)

      type(s2fil_filter), intent(in) :: filter
      logical :: sky_status

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_sky_status')
      end if 

      sky_status = filter%sky_status

    end function s2fil_filter_get_sky_status


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_filter_cl
    !
    !! Compute and return filter cl spectrum for the iscale index specified.
    !!
    !! Notes:
    !!   - Initialises a new filter_cl (s2_pl) object.  The returned 
    !!     filter_cl must be freed by the calling routine.
    !!   - Alms of filter must already be computed (which they will be since 
    !!     the optimal filters are constructed in harmonic space).
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - iscale: Index of scale to get.
    !!   - filter_cl: Object filter cl spectrum corresponding 
    !!     to the specified scale.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_filter_cl(filter, iscale) result(filter_cl)

      type(s2fil_filter), intent(in) :: filter
      integer, intent(in) :: iscale
      type(s2_pl) :: filter_cl

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_filter_cl')
      end if 

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > filter%n_scale) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_filter_get_filter_cl', comment_add='iscale outside range')
      end if

      ! Check filters constructed.
      if(.not. filter%sky_status) then
        call s2fil_error(S2FIL_ERROR_INIT, 's2fil_filter_get_filter_cl', &
        comment_add='Attempt to get filter cl before computed')
      end if

      ! Compute and initialise a pl object from the filter(iscale) sky.
      filter_cl = s2_sky_get_cl(filter%sky(iscale))

    end function s2fil_filter_get_filter_cl


    !--------------------------------------------------------------------------
    ! s2fil_filter_get_tmpl_cl
    !
    !!  Compute and return tmpl cl spectrum for the dilation corresponding 
    !!  to the iscale index specified.
    !!
    !! Notes:
    !!   - Initialises a new tmpl_cl (s2_pl) object.  The returned 
    !!     tmpl_cl must be freed by the calling routine.
    !!   - Requires dilation of template to specified scale, then alms
    !!     of dilated template must be computed, before computing cls.
    !!
    !! Variables:
    !!   - filter: Filter object to get the variable of.
    !!   - iscale: Index of scale to get.
    !!   - tmpl_cl: Object tmpl cl spectrum corresponding to dilation specified
    !!     by iscale.
    !
    !! @author J. D. McEwen
    !! @version 0.1 January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function s2fil_filter_get_tmpl_cl(filter, iscale) result(tmpl_cl)

      type(s2fil_filter), intent(in) :: filter
      integer, intent(in) :: iscale
      type(s2_pl) :: tmpl_cl

      type(s2_sky) :: temp_tmpl

      ! Check object initialised.
      if(.not. filter%init) then
        call s2fil_error(S2FIL_ERROR_NOT_INIT, 's2fil_filter_get_tmpl_spec')
      end if 

      ! Check iscale within valid range.
      if(iscale < 1 .or. iscale > filter%n_scale) then
          call s2fil_error(S2FIL_ERROR_SIZE_INVALID, &
            's2fil_filter_get_tmpl_cl', comment_add='iscale outside range')
      end if

      ! Get temporary template and dilate for scale considered.

      ! Get temporary template and dilate for scale considered without 
      ! altering original template.
      temp_tmpl = s2_sky_init(filter%tmpl)

      ! Dilate by scale corresponding to iscale.
      call s2_sky_dilate(temp_tmpl, filter%scale(iscale,1), &
           filter%scale(iscale,2), filter%norm_pres_dil)

      ! Compute alms of dilated template.
      ! Won't be computed since don't save dilated templates used for
      ! constructing optimal filter.
      call s2_sky_compute_alm(temp_tmpl, filter%lmax, filter%mmax)

      ! Convolve filter at current scale with beam if present.
      if(filter%beam_status) then
         call s2_sky_conv(temp_tmpl, filter%beam)
      end if

      ! Compute and initialise a pl object from the tmpl dilated by scale 
      ! corresponding to iscale.
      tmpl_cl = s2_sky_get_cl(temp_tmpl)
      
      ! Free temporary template used for dilation.
      call s2_sky_free(temp_tmpl)

    end function s2fil_filter_get_tmpl_cl


end module s2fil_filter_mod
