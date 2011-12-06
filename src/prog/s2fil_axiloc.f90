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
  use s2_error_mod
  use s2_sky_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_inp
  character(len=S2_STRING_LEN) :: filename_filter_data
  character(len=S2_STRING_LEN) :: filename_out_prefix

  character(len=1), parameter :: COMMENT_CHAR = '#'
  integer :: fileid = 21, iostat, fail = 0
  character(len=S2_STRING_LEN) :: line, line2
  integer :: verbosity = 5
  integer :: nside, lmax, nfil, ifil
  type(s2_sky) :: sky
  type(s2_sky), allocatable :: filter(:)
  type(s2_sky), allocatable :: mean(:)
  type(s2_sky), allocatable :: std(:)
  type(s2_sky), allocatable :: filtered(:)
  type(s2_sky), allocatable :: sig(:)
  type(s2_sky) :: tmp

  real(s2_sp), allocatable :: filter_data_theta(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_filter(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_mean(:)
  character(len=S2_STRING_LEN), allocatable :: filter_data_filename_std(:)
  real(s2_sp), allocatable :: filter_data_nstd(:)


  !----------------------------------------------------------------------------
  ! Read command line options
  !----------------------------------------------------------------------------

  call parse_options()


  !----------------------------------------------------------------------------
  ! Read data
  !----------------------------------------------------------------------------

  ! Read input sky.
  sky = s2_sky_init(filename_inp, S2_SKY_FILE_TYPE_MAP)

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

  ! Save filtered maps.
  if (verbosity > 1) then
     do ifil = 0,nfil-1
        write(line,'(a,a,i2.2,a)') trim(filename_out_prefix), &
             '_filtered_ifil', ifil, '.fits'
        call s2_sky_write_file(filtered(ifil), trim(line), S2_SKY_FILE_TYPE_MAP) 
     end do
  end if


  !----------------------------------------------------------------------------
  ! Threshold filtered maps
  !----------------------------------------------------------------------------

  ! Compute sig = (filtered - mean) / std maps.
  allocate(sig(0:nfil-1), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 's2fil_axiloc')
  end if
  do ifil = 0,nfil-1  
     tmp = s2_sky_add(filtered(ifil), mean(ifil), subtract=.true.)
     sig(ifil) = s2_sky_product(tmp, std(ifil), divide=.true.)
     call s2_sky_free(tmp)
  end do

  ! Save sig maps.
  if (verbosity > 1) then
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
  if (verbosity > 1) then
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
  ! Free memory
  !----------------------------------------------------------------------------

  call s2_sky_free(sky)
  deallocate(filter_data_theta)
  deallocate(filter_data_nstd)
  deallocate(filter_data_filename_filter)
  deallocate(filter_data_filename_mean)
  deallocate(filter_data_filename_std)
  do ifil = 0, nfil-1
     call s2_sky_free(filter(ifil))
     call s2_sky_free(mean(ifil))
     call s2_sky_free(std(ifil))
     call s2_sky_free(filtered(ifil))
     call s2_sky_free(sig(ifil))
  end do
  deallocate(filter)
  deallocate(mean)
  deallocate(std)
  deallocate(filtered)
  deallocate(sig)


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
