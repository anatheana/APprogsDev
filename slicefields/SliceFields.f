!! Field slicer
!! -
!! Command line tool for slicing and combining electric fields from
!! ADDA internal field output files. Slicing will divide the 3D voxel file into
!! 2D slices in given dimensions. The code can also apply the same slicing to several
!! fields and combine the results (combine X- and Y-polarizations, and substract incident field)
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, ADDA internal field slicing tool (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2013).
!! -
!! Antti Penttilä
!! 2013
!! Department of Physics, University of Helsinki
PROGRAM SliceFields
  !! Double precision FP_KIND type parameter from PaulV
  USE Type_Kinds
  !! File handling routines from PaulV
  USE File_Utility
  !! String routines from PaulV
  USE String_Utility
  !! My own math- and printing-related routines
  USE AP_utils
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: str_len = 128, line_len = 512, &
    fn_length = 128, dbl = FP_KIND, mem_block = 262144 ! 64^3
  
  INTEGER :: n_field_files, i, j, k, co_dir, inp_ch, out_ch, astat, &
    block_n = 0, array_length = 0, val_n = 0
  
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: co_vals, unsrt
  
  CHARACTER(LEN=fn_length) :: field_fn_X, field_fn_Y, empty_fn_X, empty_fn_Y
  
  WRITE(*,*) ""
  WRITE(*,*) "Starting SliceFields..."
  CALL handle_input()
  
  WRITE(*,*) ""
  WRITE(*,'(A49,I0)') " First pass, find distinct values of coordinate #", co_dir
  CALL first_pass()
  
  WRITE(*,*) ""
  WRITE(*,*) "Ready"
  
CONTAINS

SUBROUTINE first_pass()
  IMPLICIT NONE
  REAL(dbl), DIMENSION(3) :: temp_three
  CHARACTER(LEN=line_len) :: line
  
  inp_ch = Get_Lun()
  OPEN(inp_ch, FILE=TRIM(field_fn_X), ACTION='read', STATUS='old', IOSTAT=astat)
  IF(astat /= 0) THEN
    WRITE(*,'(A,A,A)') " Error: cannot open file '", TRIM(field_fn_X), "' for reading"
  END IF
  
  ! Allocate space for coordinate values
  CALL allocate_array()
  
  ! Read out first line headers
  READ(inp_ch, *) line
  READ(inp_ch, '(A)', IOSTAT=astat) line
  DO WHILE(astat == 0)
    val_n = val_n+1
    IF(val_n > array_length) CALL allocate_array()
    READ(line, *) temp_three
    unsrt(val_n) = temp_three(co_dir)
    READ(inp_ch, '(A)', IOSTAT=astat) line
  END DO
  CLOSE(inp_ch)
  
  WRITE(*,'(A6,I0,A8)') " read ", val_n, " records"
  CALL list_union(unsrt, val_n) ! TEE TÄMÄ
  ALLOCATE(co_vals(val_n), STAT=astat)
  IF(astat /= 0) THEN
    WRITE(*,*) "Memory allocation error"
    STOP
  END IF
  co_vals(1:val_n) = unsrt(1:val_n)
  DEALLOCATE(unsrt)
  
  WRITE(*,*) co_vals
  
END SUBROUTINE first_pass


SUBROUTINE allocate_array()
  IMPLICIT NONE
  INTEGER :: bn, ex_length
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: temp_array
  
  bn = block_n+1
  ex_length = array_length
  array_length = bn * mem_block
  
  ! First time
  IF(bn == 1) THEN
    ALLOCATE(unsrt(array_length), STAT=astat)
    IF(astat /= 0) THEN
      WRITE(*,*) "Memory allocation error"
      STOP
    END IF
  ELSE
    ALLOCATE(temp_array(ex_length), STAT=astat)
    IF(astat /= 0) THEN
      WRITE(*,*) "Memory allocation error"
      STOP
    END IF
    temp_array(1:ex_length) = unsrt(1:ex_length)
    DEALLOCATE(unsrt)
    ALLOCATE(unsrt(array_length), STAT=astat)
    IF(astat /= 0) THEN
      WRITE(*,*) "Memory allocation error"
      STOP
    END IF
    unsrt(1:ex_length) = temp_array(1:ex_length)
    DEALLOCATE(temp_array)
  END IF

END SUBROUTINE allocate_array


! alpha-version
SUBROUTINE handle_input()
  IMPLICIT NONE
  CHARACTER(LEN=64) :: arg
  
  n_field_files = COMMAND_ARGUMENT_COUNT()
  
  IF(n_field_files < 2 .OR. 5 < n_field_files) THEN
    WRITE(*,*) " Error: User must supply slicing direction (number 1-3) and 1-4 filenames to process"
    STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1,arg)
  READ(arg, *) co_dir
  CALL GET_COMMAND_ARGUMENT(2,field_fn_X)
  IF(n_field_files > 2) CALL GET_COMMAND_ARGUMENT(3,field_fn_Y)
  IF(n_field_files > 3) CALL GET_COMMAND_ARGUMENT(4,empty_fn_X)
  IF(n_field_files > 4) CALL GET_COMMAND_ARGUMENT(5,empty_fn_Y)
  
END SUBROUTINE handle_input


SUBROUTINE open_files()
  IMPLICIT NONE

  ! Open files and channels
  ! fp_in = Get_Lun()
  ! OPEN(fp_in, FILE=fn_in, STATUS='old', ACTION='read', IOSTAT=stats)
  ! IF(stats /= 0) THEN
    ! WRITE(*,*) "Error in opening file '", TRIM(fn_in), "' for reading"
    ! STOP
  ! END IF
  ! fp_out = Get_Lun()
  ! OPEN(fp_out, FILE=fn_out, STATUS='replace', ACTION='write', IOSTAT=stats)
  ! IF(stats /= 0) THEN
    ! WRITE(*,*) "Error in opening file '", TRIM(fn_out), "' for writing"
    ! STOP
  ! END IF
  
END SUBROUTINE open_files



END PROGRAM SliceFields
