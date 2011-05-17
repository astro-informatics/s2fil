!------------------------------------------------------------------------------
! s2fil_error_mod -- S2FIL library error class
!
!! Functionality to handle errors that may occur in the S2FIL library.  Public
!! S2FIL error codes are defined, with corresponding private error comments
!! and default halt execution status.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - January 2005
!
! Revisions:
!   January 2005 - Written by Jason McEwen
!------------------------------------------------------------------------------

module s2fil_error_mod

  use s2_types_mod, only: S2_STRING_LEN

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: s2fil_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: S2FIL_ERROR_NUM = 16

  integer, public, parameter :: &
    S2FIL_ERROR_NONE = 0, &
    S2FIL_ERROR_INIT = 1, &
    S2FIL_ERROR_NOT_INIT = 2, &
    S2FIL_ERROR_INIT_FAIL = 3, &
    S2FIL_ERROR_SIZE_INVALID = 4, &
    S2FIL_ERROR_MEM_ALLOC_FAIL = 5, &
    S2FIL_ERROR_FILTER_TYPE_INVALID = 6, &
    S2FIL_ERROR_SCALE_TYPE_INVALID = 7, &
    S2FIL_ERROR_DIV_BY_ZERO = 8, &
    S2FIL_ERROR_PROG_INPUT_INVALID = 9, &
    S2FIL_ERROR_BEAM_NOT_DEF = 10, &
    S2FIL_ERROR_FIL_FILE_EXISTS = 11, &
    S2FIL_ERROR_FIL_FILE_INVALID = 12,&
    S2FIL_ERROR_FIELD_TR_DEF = 13, &
    S2FIL_ERROR_FIELD_FILE_EXISTS = 14, &
    S2FIL_ERROR_FIELD_FILE_INVALID = 15


  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.  
  !! Comment associated with each error type.
  character(len=S2_STRING_LEN), parameter :: &
    error_comment(S2FIL_ERROR_NUM) = &
      (/ & 
      'No error                                                                 ', &
      'Object already initialised                                               ', &
      'Object not initialised                                                   ', &
      'Object initialisation failed                                             ', &
      'Sizes invalid                                                            ', &
      'Memory allocation failed                                                 ', &
      'Invalid filter type                                                      ', &
      'Invalid scale type                                                       ', &
      'Divide by zero avoided                                                   ', &
      'Invalid program command line option                                      ', &
      'Beam not defined                                                         ', &
      'Filter file already exists                                               ', &
      'Filter file invalid                                                      ', &
      'Filtered field already defined                                           ', &
      'Field file already exists                                                ', &
      'Field file invalid                                                       ' &
      /) 

  !! Default program halt status of each error type.
  logical, parameter :: &
    halt_default(S2FIL_ERROR_NUM) = &
      (/ &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.  /)
  
  
  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! s2fil_error
    ! 
    !! Display error message corresponding to error_code and halt program 
    !! execution if required.
    !!
    !! Variables:
    !!   - error_code:
    !!   - [procedure]: Procedure name where s2_error called from.  Displayed 
    !!     when error message printed to screen.
    !!   - [comment_add]: If present, additional comment to append to default 
    !!     error comment.
    !!   - [comment_out]: If present the error comment is copied to comment_out
    !!     on output.
    !!   - [halt_in]: If present overrides default halt value.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine s2fil_error(error_code, procedure, comment_add, &
      comment_out, halt_in)

      integer, intent(in) :: error_code
      character(len=*), intent(in), optional :: procedure, comment_add
      character(len=*), intent(inout), optional :: comment_out
      logical, intent(in), optional :: halt_in

      logical :: halt
      character(len=*), parameter :: comment_prefix = 'S2FIL_ERROR: '

      !---------------------------------------
      ! Display error message
      !---------------------------------------

      if(present(procedure)) then

        if(present(comment_add)) then
          write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            '''', &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            ''''
        end if
 
     else

        if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
            trim(error_comment(error_code+1)), &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
        end if

      end if

      ! Copy error comment if comment_out present.
      if(present(comment_out)) comment_out = error_comment(error_code+1)

      !---------------------------------------
      ! Halt program execution if required
      !---------------------------------------
      
      if( present(halt_in) ) then
        halt = halt_in
      else
        halt = halt_default(error_code+1)
      end if

      if( halt ) then
        write(*,'(a,a,a,a,a)') comment_prefix, &
          '  Halting program execution ', &
          'due to error ''', trim(error_comment(error_code+1)), ''''
        stop
      end if

    end subroutine s2fil_error


end module s2fil_error_mod
