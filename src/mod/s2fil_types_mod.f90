!------------------------------------------------------------------------------
! s2fil_types_mod -- S2FIL library types class
!
!! Definition of intrinsic types and constants used in the s2fil library.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 January 2005
!
! Revisions:
!   January 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2fil_types_mod

  implicit none

  private


  ! --------------------------------------
  ! Intrinsic type definitions
  ! --------------------------------------

  !! Size of second dimension of scale(/dilation).
  integer, public, parameter :: S2FIL_SCALE_DIM2_SIZE = 2

  !! Filter type: Matched filter (symmetric)
!  integer, public, parameter :: S2FIL_FILTER_TYPE_MFS = 1

  !! Filter type: Matched filter (directional)
!  integer, public, parameter :: S2FIL_FILTER_TYPE_MFD = 2

  !! Filter type: Scale adaptive filter (symmetric)
!  integer, public, parameter :: S2FIL_FILTER_TYPE_SAFS = 3

  !! Filter type: Scale adaptive filter (directional)
!  integer, public, parameter :: S2FIL_FILTER_TYPE_SAFD = 4

  !! Default norm preserving dilation status.
  logical, public, parameter :: S2FIL_NORM_PRES_DIL_DEFAULT = .false.


end module s2fil_types_mod
