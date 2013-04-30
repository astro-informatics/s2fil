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
!!   - [-file_type_in file_type_in]: Type of input file (map; alm; sky).
!!   - [-filter_mask filename_mask]: Name of optional mask.
!!   - [-filter_data filename_filter_data]: Name of filter data file.
!!   - [-theta_filter_adj theta_filter_adj (degrees)]: Maximum angular
!!     separation between filters that are considered to be adjacent.
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
  use s2_error_mod
  use s2_sky_mod
  use s2_vect_mod
  use pix_tools, only: ang2pix_ring, ang2pix_nest

  implicit none

  interface 
    function abs_pix(x) result (val)
      use s2_types_mod
      real(s2_sp), intent(in) :: x
      real(s2_sp) :: val
    end function abs_pix
  end interface

  character(len=*), parameter ::  MAP_FILE = 'map'
  character(len=*), parameter ::  ALM_FILE = 'alm'
  character(len=*), parameter ::  SKY_FILE = 'sky'
  integer :: file_type_in = S2_SKY_FILE_TYPE_MAP
  character(len=S2_STRING_LEN) :: file_type_in_str = MAP_FILE

  character(len=S2_STRING_LEN) :: filename_inp
  character(len=S2_STRING_LEN) :: filename_mask
  logical :: apply_mask = .false.
  character(len=S2_STRING_LEN) :: filename_filter_data
  character(len=S2_STRING_LEN) :: filename_out_prefix

  character(len=1), parameter :: COMMENT_CHAR = '#'
  integer :: fileid = 21, iostat, fail = 0
  character(len=S2_STRING_LEN) :: line, line2
  integer :: verbosity = 5
  integer :: nside, lmax, nfil, ifil, iadj, ireg, iregadj
  type(s2_sky) :: sky
  type(s2_sky) :: sky_mask, sky_temp
  type(s2_sky), allocatable :: filter(:)
  type(s2_sky), allocatable :: mean(:)
  type(s2_sky), allocatable :: std(:)
  type(s2_sky), allocatable :: filtered(:)
  type(s2_sky), allocatable :: sig(:)
  type(s2_sky), allocatable :: mask(:)
  type(s2_sky) :: tmp

  real(s2_sp), allocatable :: filter_data_theta(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_filter(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_mean(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_std(:)
  real(s2_sp), allocatable :: filter_data_nstd(:)

  integer, parameter :: NCENTRES_MAX = 500
  integer, allocatable :: ncentres(:)
  real(s2_dp), allocatable :: centres_theta_tmp(:)
  real(s2_dp), allocatable :: centres_phi_tmp(:)
  real(s2_dp), allocatable :: centres_radius_tmp(:)
  real(s2_dp), allocatable :: centres_theta(:,:)
  real(s2_dp), allocatable :: centres_phi(:,:)
  real(s2_dp), allocatable :: centres_radius(:,:)
  real(s2_dp) :: peak_radius, peak_radius_adj

  type(s2_vect) :: vec0, vec1
  real(s2_dp) :: dot
  real(s2_sp) :: reg_sep_ang
  logical :: discard
  integer :: isource, nsource, iamp
  real(s2_sp) :: theta_filter_adj = 5.0
  real(s2_sp) :: max_amp, max_amp_nearby, amp, sig_max
  logical, allocatable :: adj(:)
  real(s2_sp), allocatable :: regions_amp(:)
  real(s2_sp), allocatable :: regions_size(:)
  real(s2_dp), allocatable :: regions_theta(:)
  real(s2_dp), allocatable :: regions_phi(:)
  real(s2_dp), allocatable :: regions_sig(:)
  real(s2_dp), allocatable :: regions_sig_radius(:)
  real(s2_dp), allocatable :: regions_mask_overlap(:)
  real(s2_dp), allocatable :: regions_prior_size_up(:)
  real(s2_dp), allocatable :: regions_prior_size_down(:)
  integer :: prior_size_up_ifil, prior_size_down_ifil, ifil_p
  integer :: masked
  real(s2_sp) :: size_prior_radius_factor = 0.10


  !----------------------------------------------------------------------------
  ! Read command line options
  !----------------------------------------------------------------------------

  call parse_options()


  !----------------------------------------------------------------------------
  ! Read data
  !----------------------------------------------------------------------------

  ! Set input sky file type.
  select case (trim(file_type_in_str))
    case (MAP_FILE)
       file_type_in = S2_SKY_FILE_TYPE_MAP
    case (ALM_FILE)
       file_type_in = S2_SKY_FILE_TYPE_ALM
    case (SKY_FILE)
       file_type_in = S2_SKY_FILE_TYPE_SKY
    case default
       call s2_error(S2_ERROR_SKY_FILE_INVALID, 's2_axiconv', &
         comment_add='Invalid file type option')
  end select

  ! Read input sky.
  sky = s2_sky_init(filename_inp, file_type_in)

  ! Read filter data file.
  open(fileid, file=filename_filter_data, form='formatted', status='old')

  ! Count number of filters.
  nfil = 0
  do 
     read(fileid,'(a)',iostat=iostat) line
     if (iostat < 0) exit
     if (line(1:1) /= COMMENT_CHAR) nfil = nfil + 1
  end do
  rewind(fileid)

  ! Allocate space for filter data.
  allocate(filter_data_theta(0:nfil-1), stat=fail)
  allocate(filter_data_filename_filter(0:nfil-1), stat=fail)
  allocate(filter_data_filename_mean(0:nfil-1), stat=fail)
  allocate(filter_data_filename_std(0:nfil-1), stat=fail)
  allocate(filter_data_nstd(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if

  ! Read filter data.
  ifil = 0
  do 
     read(fileid,'(a)',iostat=iostat) line
     if (iostat < 0) exit
     if (line(1:1) /= COMMENT_CHAR) then
        read(line,'(f4.1, 2X, f4.2, 2X, a)') &
             filter_data_theta(ifil), &
             filter_data_nstd(ifil), &
             line2
        filter_data_filename_filter(ifil) = &
             trim(line2(1:index(line2,',')-1));
        filter_data_filename_mean(ifil) = &
             trim(line2(index(line2,',')+2:index(line2,',',back=.true.)-1));
        filter_data_filename_std(ifil) = &
             trim(line2(index(line2,',',back=.true.)+2:S2_STRING_LEN));
        ifil = ifil + 1
     end if
  end do

  ! Write filter data.
  if (verbosity >= 5) then
     do ifil = 0,nfil-1
        write(*,'(a,i2,a,f4.1)') 'filter_data_theta(', ifil, ')           = ', &
             filter_data_theta(ifil)
        write(*,'(a,i2,a,f4.1)') 'filter_data_nstd(', ifil, ')            = ', &
             filter_data_nstd(ifil)
        write(*,'(a,i2,a,a)') 'filter_data_filename_filter(', ifil, ') = ', &
             trim(filter_data_filename_filter(ifil))
        write(*,'(a,i2,a,a)') 'filter_data_filename_mean(', ifil, ')   = ', &
             trim(filter_data_filename_mean(ifil))
        write(*,'(a,i2,a,a)') 'filter_data_filename_std(', ifil, ')    = ', &
             trim(filter_data_filename_std(ifil))
     end do
  end if
  close(fileid)

  ! Read filters.
  allocate(filter(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1
     filter(ifil) = s2_sky_init(trim(filter_data_filename_filter(ifil)), &
          S2_SKY_FILE_TYPE_SKY)
  end do

  ! Read mean maps.
  allocate(mean(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1
     mean(ifil) = s2_sky_init(trim(filter_data_filename_mean(ifil)), &
          S2_SKY_FILE_TYPE_MAP)
  end do

  ! Read std maps.
  allocate(std(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1
     std(ifil) = s2_sky_init(trim(filter_data_filename_std(ifil)), &
          S2_SKY_FILE_TYPE_MAP)
  end do

  ! Read input mask.
  sky_mask = s2_sky_init(filename_mask, S2_SKY_FILE_TYPE_MAP)


  !----------------------------------------------------------------------------
  ! Apply matched filters
  !----------------------------------------------------------------------------

  ! Compute harmonic coefficients of input map.
  call s2_sky_compute_alm(sky, lmax, lmax)
  
  ! Filter maps.
  allocate(filtered(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1
     filtered(ifil) = s2_sky_axiconv(sky, filter(ifil), compute_map=.true.)
  end do

  ! Mask filtered maps.
  if (apply_mask) then

     ! Apply mask to each filtered map.
     do ifil = 0,nfil-1
        sky_temp = s2_sky_product(filtered(ifil), sky_mask)
        call s2_sky_free(filtered(ifil))
        filtered(ifil) = s2_sky_init(sky_temp)
        call s2_sky_free(sky_temp)
     end do

  end if

  ! Save filtered maps.
  if (verbosity >= 1) then
     do ifil = 0,nfil-1
        write(line,'(a,a,i2.2,a)') trim(filename_out_prefix), &
             '_filtered_ifil', ifil, '.fits'
        call s2_sky_write_file(filtered(ifil), trim(line), S2_SKY_FILE_TYPE_MAP) 
     end do
  end if


  !----------------------------------------------------------------------------
  ! Compute significance maps and threshold
  !----------------------------------------------------------------------------

  ! Compute sig = (filtered - mean) / std maps.
  allocate(sig(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1  
     tmp = s2_sky_add(filtered(ifil), mean(ifil), subtract=.true.)
     call s2_sky_fun(tmp, abs_pix)
     sig(ifil) = s2_sky_product(tmp, std(ifil), divide=.true.)
     call s2_sky_free(tmp)
  end do

  ! Save sig maps.
  if (verbosity >= 1) then
     do ifil = 0,nfil-1
        write(line,'(a,a,i2.2,a)') trim(filename_out_prefix), &
             '_sig_ifil', ifil, '.fits'
        call s2_sky_write_file(sig(ifil), trim(line), S2_SKY_FILE_TYPE_MAP) 
     end do
  end if

  ! Threshold sig maps.
  do ifil = 0,nfil-1  
     call s2_sky_thres_abs(sig(ifil), filter_data_nstd(ifil))
  end do

  ! Save thresholded sig maps.
  if (verbosity >= 1) then
     do ifil = 0,nfil-1
        write(line,'(a,a,i2.2,f3.2,a,i2.2,a)') trim(filename_out_prefix), &
             '_sigthres_nstd', &
             floor(filter_data_nstd(ifil)), &
             filter_data_nstd(ifil)-real(floor(filter_data_nstd(ifil)),s2_sp), &
             '_ifil', ifil, '.fits'
        call s2_sky_write_file(sig(ifil), trim(line), S2_SKY_FILE_TYPE_MAP) 
     end do
  end if


  !----------------------------------------------------------------------------
  ! Find localised regions at each scale
  !----------------------------------------------------------------------------

  ! Allocate space for masks.
  allocate(mask(0:nfil-1), stat=fail)
  allocate(ncentres(0:nfil-1), stat=fail)
  allocate(centres_theta(0:nfil-1, 0:NCENTRES_MAX-1), stat=fail)
  allocate(centres_phi(0:nfil-1, 0:NCENTRES_MAX-1), stat=fail)
  allocate(centres_radius(0:nfil-1, 0:NCENTRES_MAX-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if

  ! Localise peak regions for each filter scale.
  do ifil = 0,nfil-1  

     ! Find peak regions.
     peak_radius = filter_data_theta(ifil) / 180.0 * PI
     call s2_sky_thres_peaks(sig(ifil), peak_radius, mask(ifil), &
          ncentres(ifil), centres_theta_tmp, centres_phi_tmp, &
          centres_radius_tmp)
     
     ! Save peak region data.
     if (ncentres(ifil) > NCENTRES_MAX) then
        call s2_error(S2_ERROR_MEM_OUT_OF_BOUNDS, 's2fil_axiloc', &
             comment_add='Number of centres exceeds limit')
     end if    
     centres_theta(ifil, 0:ncentres(ifil)-1) = &
          centres_theta_tmp(0:ncentres(ifil)-1)
     centres_phi(ifil, 0:ncentres(ifil)-1) = &
          centres_phi_tmp(0:ncentres(ifil)-1)
     centres_radius(ifil, 0:ncentres(ifil)-1) = &
          centres_radius_tmp(0:ncentres(ifil)-1)
     deallocate(centres_theta_tmp, centres_phi_tmp, centres_radius_tmp)

  end do

  ! Save region mask maps.
  if (verbosity >= 1) then
     do ifil = 0,nfil-1
        write(line,'(a,a,i2.2,a)') trim(filename_out_prefix), &
             '_regionmask_ifil', ifil, '.fits'
        call s2_sky_write_file(mask(ifil), trim(line), S2_SKY_FILE_TYPE_MAP) 
     end do
  end if

  ! Save source parameters for each filter scale to output files.
  if (verbosity >= 5) then
     do ifil = 0,nfil-1  

        write(line,'(a,a,i2.2,a)') trim(filename_out_prefix), &
             '_sources_ifil', ifil, '.txt'
        open(unit=fileid, file=line, status='replace', action='write')
        write(fileid,'(a,a)') COMMENT_CHAR, ' Localised source positions'
        write(fileid,'(a,a)') COMMENT_CHAR, ' Written by s2fil_axiloc'
        write(fileid,'(a)') COMMENT_CHAR
        write(fileid,'(a,i20)') 'n_sources= ', ncentres(ifil)
        do ireg = 0, ncentres(ifil)-1

           ! Get amplitude of filtered field and significance level at source position.
           if(s2_sky_get_pix_scheme(filtered(ifil)) == S2_SKY_RING) then
              call ang2pix_ring(s2_sky_get_nside(filtered(ifil)), &
                   centres_theta(ifil,ireg), centres_phi(ifil,ireg), iamp)
           else if(s2_sky_get_pix_scheme(filtered(ifil)) == S2_SKY_NEST) then
              call ang2pix_nest(s2_sky_get_nside(filtered(ifil)), &
                   centres_theta(ifil,ireg), centres_phi(ifil,ireg), iamp)
           else
              call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2fil_axiloc')
           end if
           amp = s2_sky_get_map_pix(filtered(ifil), iamp)
           sig_max = s2_sky_get_map_pix(sig(ifil), iamp)

           ! Write detected source.
           write(fileid, '(a)') COMMENT_CHAR
           write(fileid,'(a,e20.10)') 'amplitude= ', amp
           write(fileid,'(a,e24.10)') 'alpha= ', centres_phi(ifil,ireg)
           write(fileid,'(a,e24.10)') 'beta=  ', centres_theta(ifil,ireg)
           write(fileid,'(a,e24.10)') 'gamma= ', 0.0
           write(fileid,'(a,e24.10)') 'size=  ', filter_data_theta(ifil) / 180 * PI
           write(fileid,'(a,e24.10)') 'sig=   ', sig_max
           write(fileid,'(a,e19.10)') 'sig_radius= ', centres_radius(ifil,ireg)
        end do
        close(fileid)

     end do
  end if


  !----------------------------------------------------------------------------
  ! Locate final regions by looking across scales
  !----------------------------------------------------------------------------

  ! Allocate memory to store sources found.
  allocate(adj(0:nfil-1), stat=fail)
  allocate(regions_amp(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_size(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_theta(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_phi(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_sig(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_sig_radius(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_mask_overlap(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_prior_size_up(0:NCENTRES_MAX-1), stat=fail)
  allocate(regions_prior_size_down(0:NCENTRES_MAX-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if

  ! Locate sources.
  isource = 0
  do ifil = 0, nfil-1

     ! Find adjacent scales.
     adj(0:nfil-1) = .false.
     do iadj = 0, nfil-1
        if (iadj /= ifil) then
           if ( abs(filter_data_theta(ifil) - filter_data_theta(iadj)) &
                <= theta_filter_adj) then
              adj(iadj) = .true.
           end if
        end if
     end do

     ! Write adjacent filter list.
     if (verbosity >= 5) then
        write(*,*)
        write(*,'(a,i4)') 'ifil              = ', ifil
        write(*,'(a,f4.1)') 'filter_data_theta = ', filter_data_theta(ifil)
        write(*,'(a,a)') '  iadj adj(iadj)'
        do iadj = 0, nfil-1
           write(*,'(i6,l10)') iadj, adj(iadj)
        end do
     end if

     ! Set peak radius for current filter.
     peak_radius = filter_data_theta(ifil) / 180.0 * PI

     ! Consider each region for the current filter.
     do ireg = 0, ncentres(ifil)-1

        ! Find max absolute amplitude of filtered field in region of
        ! current filter.
        max_amp = s2_sky_region_max(filtered(ifil), peak_radius, &
             centres_theta(ifil,ireg), centres_phi(ifil,ireg))

        ! Get amplitude of filtered field at source position.
        if(s2_sky_get_pix_scheme(filtered(ifil)) == S2_SKY_RING) then
           call ang2pix_ring(s2_sky_get_nside(filtered(ifil)), &
                centres_theta(ifil,ireg), centres_phi(ifil,ireg), iamp)
        else if(s2_sky_get_pix_scheme(filtered(ifil)) == S2_SKY_NEST) then
           call ang2pix_nest(s2_sky_get_nside(filtered(ifil)), &
                centres_theta(ifil,ireg), centres_phi(ifil,ireg), iamp)
        else
           call s2_error(S2_ERROR_SKY_PIX_INVALID, 's2fil_axiloc')
        end if
        amp = s2_sky_get_map_pix(filtered(ifil), iamp)
        sig_max = s2_sky_get_map_pix(sig(ifil), iamp)

        ! Write current region.
        if (verbosity >= 5) then
           write(*,'(a)') 'Candidate region:'
           write(*,'(a,i10)') '  ireg       = ', ireg
           write(*,'(a,f10.1)') '  size       = ', filter_data_theta(ifil)
           write(*,'(a,f10.1)') '  theta      = ', centres_theta(ifil,ireg) / PI * 180
           write(*,'(a,f10.1)') '  phi        = ', centres_phi(ifil,ireg) / PI * 180
           write(*,'(a,e10.4)') '  max_amp    = ', max_amp
           write(*,'(a,e10.4)') '  amp        = ', amp
           write(*,'(a,e10.4)') '  sig        = ', sig_max
           write(*,'(a,f10.1)') '  sig_radius = ', centres_radius(ifil,ireg) / PI * 180
        end if

        ! Find nearby regions of adjacent scales.
        discard = .false.
        do iadj = 0, nfil-1
           if (iadj /= ifil .and. adj(iadj) .and. (.not. discard)) then

              ! Set peak radius for adjacent filter.
              peak_radius_adj = filter_data_theta(iadj) / 180.0 * PI

              ! Consider each region for the adjacent filters.
              do iregadj = 0, ncentres(iadj)-1
                 
                 ! Compute separation between centres of regions.
                 vec0 = s2_vect_init(1.0, &
                      real(centres_theta(ifil,ireg),s2_sp), &
                      real(centres_phi(ifil,ireg),s2_sp))                 
                 vec1 = s2_vect_init(1.0, &
                      real(centres_theta(iadj,iregadj),s2_sp), &
                      real(centres_phi(iadj,iregadj),s2_sp))
                 dot = s2_vect_dot(vec0, vec1)         
                 call s2_vect_free(vec0)
                 call s2_vect_free(vec1)
                 if (dot > 1d0) dot = 1d0  ! Remove numerical noise that
                                           ! could cause acos to fail.
                 reg_sep_ang = acos(dot)
                 reg_sep_ang = reg_sep_ang/ PI * 180.0

                 ! If regions nearby then find peak in adjacent filtered field.
                 if (reg_sep_ang < &
                      (filter_data_theta(ifil) + filter_data_theta(iadj))) then

                    ! Find max amplitude of filtered field in region of 
                    ! adjacent filter.
                    max_amp_nearby = &
                         s2_sky_region_max(filtered(iadj), peak_radius_adj, &
                         centres_theta(iadj,iregadj), centres_phi(iadj,iregadj))

                    ! If max amplitude greater than max amplitude of region 
                    ! for current filter, and amplitudes have the same sign, 
                    ! then discard.
                    if (abs(max_amp) < abs(max_amp_nearby) .and. ((max_amp * max_amp_nearby) > 0)) then
                       discard = .true.
                       ! Write region that rejected candidate.
                       if (verbosity >= 5) then
                          write(*,'(a)') 'Candidate region rejected by:'
                          write(*,'(a,i10)') '  iadj       = ', iadj
                          write(*,'(a,f10.1)') '  size       = ', filter_data_theta(iadj)
                          write(*,'(a,i10)') '  iregadj    = ', iregadj
                          write(*,'(a,f10.1)') '  theta      = ', centres_theta(iadj,iregadj) / PI * 180
                          write(*,'(a,f10.1)') '  phi        = ', centres_phi(iadj,iregadj) / PI * 180
                          write(*,'(a,e10.4)') '  max_amp    = ', max_amp_nearby
                          write(*,'(a,f10.1)') '  sig_radius = ', centres_radius(iadj,iregadj) / PI * 180
                       end if
                       exit
                    end if

                 end if

              end do

           end if
        end do

        ! Save region at this scale if not discarded.
        if (.not. discard) then

           ! Check do not overflow buffer.
           if (isource >= NCENTRES_MAX) then
              call s2_error(S2_ERROR_MEM_OUT_OF_BOUNDS, 's2fil_axiloc', &
                   comment_add='Number of sources found exceeds limit')
           end if

           ! Compute theta prior edges as first adjacent scales where
           ! don't see candidate at similar location.

           ! Look up...
           if (ifil == nfil-1) then
              prior_size_up_ifil = ifil
           else
              do ifil_p = ifil+1, nfil-1
                 masked = nint(s2_sky_region_max(mask(ifil_p),  &
                      size_prior_radius_factor*peak_radius, &
                      centres_theta(ifil,ireg), centres_phi(ifil,ireg)))
                 prior_size_up_ifil = ifil_p
                 if (masked /= 1) exit
              end do
           end if

           ! Look down...
           if (ifil == 0) then
              prior_size_down_ifil = ifil
           else
              do ifil_p = ifil-1, 0, -1
                 masked = nint(s2_sky_region_max(mask(ifil_p), &
                      size_prior_radius_factor*peak_radius, &
                      centres_theta(ifil,ireg), centres_phi(ifil,ireg)))
                 prior_size_down_ifil = ifil_p
                 if (masked /= 1) exit
              end do
           end if

           ! Save region.
           regions_amp(isource) = amp
           regions_size(isource) = filter_data_theta(ifil) / 180 * PI
           regions_theta(isource) = centres_theta(ifil,ireg)
           regions_phi(isource) = centres_phi(ifil,ireg)
           regions_sig(isource) = sig_max
           regions_sig_radius(isource) = centres_radius(ifil,ireg) 
           regions_prior_size_up(isource) = filter_data_theta(prior_size_up_ifil) / 180.0 * PI
           regions_prior_size_down(isource) = filter_data_theta(prior_size_down_ifil) / 180.0 * PI
           isource = isource + 1

        end if

     end do

  end do
  nsource = isource

  ! Calculate overlap of each region with mask.
  do isource = 0, nsource-1
     regions_mask_overlap(isource) = s2_sky_mask_overlap(sky_mask, &
          regions_theta(isource), regions_phi(isource), &
          real(regions_size(isource),s2_dp))
  end do

  ! Write source parameters.
  if (verbosity >= 5) then
     write(*,*) 
     write(*,'(a)') 'Detected sources:'
     do isource = 0, nsource-1
        write(*,*)      
        write(*,'(a,i10)') '  isource      = ', isource
        write(*,'(a,f10.1)') '  size         = ', regions_size(isource) / PI * 180
        write(*,'(a,f10.1)') '  theta        = ', regions_theta(isource) / PI * 180
        write(*,'(a,f10.1)') '  phi          = ', regions_phi(isource) / PI * 180
        write(*,'(a,e10.4)') '  amp          = ', regions_amp(isource)
        write(*,'(a,e10.4)') '  sig          = ', regions_sig(isource)
        write(*,'(a,e10.4)') '  sig_radius   = ', regions_sig_radius(isource) / PI * 180
        write(*,'(a,e10.4)') '  mask_overlap = ', regions_mask_overlap(isource)
        write(*,'(a,e10.4)') '  size_lo      = ', regions_prior_size_down(isource) / PI * 180
        write(*,'(a,e10.4)') '  size_hi      = ', regions_prior_size_up(isource) / PI * 180

     end do
  end if

  ! Save source parameters to output file.
  write(line,'(a,a)') trim(filename_out_prefix), &
       '_sources.txt'
  open(unit=fileid, file=line, status='replace', action='write')
  write(fileid,'(a,a)') COMMENT_CHAR, ' Localised source positions'
  write(fileid,'(a,a)') COMMENT_CHAR, ' Written by s2fil_axiloc'
  write(fileid,'(a)') COMMENT_CHAR
  write(fileid,'(a,i20)') 'n_sources= ', nsource
  do isource = 0, nsource-1
     write(fileid, '(a)') COMMENT_CHAR
     write(fileid,'(a,e20.10)') 'amplitude= ', regions_amp(isource)
     write(fileid,'(a,e24.10)') 'alpha= ', regions_phi(isource)
     write(fileid,'(a,e24.10)') 'beta=  ', regions_theta(isource)
     write(fileid,'(a,e24.10)') 'gamma= ', 0.0
     write(fileid,'(a,e24.10)') 'size=  ', regions_size(isource)
     write(fileid,'(a,e24.10)') 'sig=   ', regions_sig(isource)
     write(fileid,'(a,e19.10)') 'sig_radius= ', regions_sig_radius(isource)
     write(fileid,'(a,e17.10)') 'mask_overlap= ', regions_mask_overlap(isource)
     write(fileid,'(a,e22.10)') 'size_lo= ', regions_prior_size_down(isource)
     write(fileid,'(a,e22.10)') 'size_hi= ', regions_prior_size_up(isource)       
  end do
  close(fileid)


  !----------------------------------------------------------------------------
  ! Free memory
  !----------------------------------------------------------------------------

  call s2_sky_free(sky)
  call s2_sky_free(sky_mask)
  deallocate(filter_data_theta)
  deallocate(filter_data_nstd)
  deallocate(filter_data_filename_filter)
  deallocate(filter_data_filename_mean)
  deallocate(filter_data_filename_std)
  deallocate(ncentres, centres_theta, centres_phi, centres_radius)
  deallocate(adj)
  deallocate(regions_amp, regions_sig, regions_sig_radius)
  deallocate(regions_mask_overlap)
  deallocate(regions_size, regions_theta, regions_phi)
  deallocate(regions_prior_size_down, regions_prior_size_up)
  do ifil = 0, nfil-1
     call s2_sky_free(filter(ifil))
     call s2_sky_free(mean(ifil))
     call s2_sky_free(std(ifil))
     call s2_sky_free(filtered(ifil))
     call s2_sky_free(sig(ifil))
     call s2_sky_free(mask(ifil))
  end do
  deallocate(filter)
  deallocate(mean)
  deallocate(std)
  deallocate(filtered)
  deallocate(sig)
  deallocate(mask)


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
              '[-file_type_in file_type_in (map; alm; sky)]'
            write(*,'(a,a)') '                    ', &
              '[-mask filename_mask]'
            write(*,'(a,a)') '                    ', &
              '[-apply_mask apply_mask (optional)]'
            write(*,'(a,a)') '                    ', &
              '[-filter_data filename_filter_data]'
            write(*,'(a,a)') '                    ', &
              '[-theta_filter_adj theta_filter_adj (degrees)]'
            write(*,'(a,a)') '                    ', &
              '[-size_prior_radius_factor size_prior_radius_factor]'
!            write(*,'(a,a)') '                    ', &
!              '[-min_peak_area min_peak_area (steradians)]'
            write(*,'(a,a)') '                    ', &
              '[-nside nside]'
            write(*,'(a,a)') '                    ', &
              '[-lmax lmax]'
            write(*,'(a,a)') '                    ', &
              '[-out out_prefix]'
            stop

          case ('-inp')
            filename_inp = trim(arg) 

          case ('-file_type_in')
            file_type_in_str = trim(arg) 

          case ('-mask')
            filename_mask = trim(arg)

          case ('-apply_mask')
            read(arg,*) apply_mask

          case ('-filter_data')
            filename_filter_data = trim(arg)

          case ('-theta_filter_adj')
            read(arg,*) theta_filter_adj

          case ('-size_prior_radius_factor')
            read(arg,*) size_prior_radius_factor

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


!---------------------------------------------------------------------
! abs_pix
!
!! Function to take absolute value of a map pixel value.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   December 2012 - Written by Jason McEwen 
!---------------------------------------------------------------------

function abs_pix(x) result (val)
  
  use s2_types_mod
  
  implicit none
  
  real(s2_sp), intent(in) :: x
  real(s2_sp) :: val
  
  val = abs(x)
  
end function abs_pix
