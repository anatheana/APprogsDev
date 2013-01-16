!! AP_Quicksort
!! Implements recursive in-place quicksort and list partition algorithms,
!! as they are presented in Wikipedia (June 2012). Implements them to integer
!! and real arrays as well as for integer and real tables.
!! This implementation is quite fast, perhaps due to the design
!! where generic interface wrappers are called but the computation
!! is done inside simpler internal functions.
!! -
!! Can be used with Type_Kinds module by PaulV, default is without.
!! -
!! Two public interfaces -- AP_Qsort and AP_partition.
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Fortran 95 implementation of the Quicksort algorithm (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2012).
!! -
!! Antti Penttilä
!! 2012
!! Department of Physics, University of Helsinki
MODULE AP_quicksort
  ! Uncomment the next line out if you have Type_Kinds module
!	USE Type_Kinds
  
	IMPLICIT NONE

	PUBLIC
	
  !! Public interface AP_Qsort implements the quicksort algorithm to
  !! integer arrays and tables, and to real array and tables.
  !! Different versions of the interface are:
  !! -
  !! subroutine AP_Qsort(vec, LOW, HIGH), where
  !! vec is either real or integer array to be sorted.
  !! Optional arguments LOW and HIGH can be used to give the start index in vec array
  !! and the end index in vec array. By default, the whole array is sorted.
  !! -
  !! OR
  !! -
  !! subroutine AP_Qsort(mat, LOW, HIGH, COLUMNS), where
  !! mat is either real or integer table (i.e. matrix) to be sorted.
  !! Extra optional argument COLUMNS gives the indices of columns that are used in sorting.
  !! Default is to sort using all the columns, and in order of first to last column.
	INTERFACE AP_Qsort
	  MODULE PROCEDURE Qsort_real, Qsort_int, Qsort_real_table, Qsort_int_table
	END INTERFACE AP_Qsort
  
  !! Public interface AP_partition implements partition of array or table
  !! using a given pivot element. Versions of the interface are:
  !! -
  !! function AP_partition(vec, pivot, LOW, HIGH), where
  !! vec is either real or integer array to be partitioned.
  !! pivot is the index of the pivot value in vec by which the array will be partitioned.
  !! Optional arguments LOW and HIGH can be used to give the start index in vec array
  !! and the end index in vec array. By default, the whole array is partitioned.
  !! Function will return an index of new pivot value. Array values from indices LOW to pivot will
  !! be smaller than or equal to value in pivot, and values from indices pivot+1 to HIGH will be larger.
  !! -
  !! OR
  !! -
  !! function AP_partition(mat, pivot, LOW, HIGH, COLUMNS), where
  !! mat is either real or integer table (i.e. matrix) to be partitioned.
  !! Extra optional argument COLUMNS gives the indices of columns that are used in partitioning.
  !! Default is to partition using all the columns, and in order of first to last column.
  INTERFACE AP_partition
	  MODULE PROCEDURE partition_real, partition_int, partition_int_table, partition_real_table
	END INTERFACE AP_partition
  
  ! Private interface
  INTERFACE is_vec_smaller
    MODULE PROCEDURE is_vec_smaller_int, is_vec_smaller_real
  END INTERFACE is_vec_smaller
  
  !! Use if no Type_Kinds module
  !! define single precision reals
!  INTEGER, PARAMETER :: internal_double_type =  SELECTED_REAL_KIND(6)
  !! Use if no Type_Kinds module
  !! define double precision reals
  INTEGER, PARAMETER :: internal_double_type = SELECTED_REAL_KIND(15)
  !! Use with Type_Kinds module
!  INTEGER, PARAMETER :: internal_double_type = FP_Kind
  
  ! Private variables
  INTEGER :: i, temp, pv_value, cols_n
  REAL(internal_double_type) :: rtemp, rpv_value
  INTEGER, DIMENSION(:), ALLOCATABLE :: tempvec, pv_valuevec, colsvec
  REAL(internal_double_type), DIMENSION(:), ALLOCATABLE :: rtempvec, rpv_valuevec

	PRIVATE :: Qsort_int_rec, Qsort_real_rec, partition_int_privat, partition_real_privat, &
    i, temp, pv_value, rtemp, rpv_value, is_vec_smaller, is_vec_smaller_int, is_vec_smaller_real, &
    tempvec, pv_valuevec, colsvec, rtempvec, rpv_valuevec

CONTAINS

! Public routines

! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for integer quicksort
SUBROUTINE Qsort_int(vec, LOW, HIGH)
  IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER :: n, left, right
  
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(vec)
  END IF

  IF(left >= right) RETURN
  
  CALL Qsort_int_rec(vec, left, right)
  
END SUBROUTINE Qsort_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for real quicksort
SUBROUTINE Qsort_real(vec, LOW, HIGH)
  IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER :: n, left, right
  
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(vec)
  END IF

  IF(left >= right) RETURN
  
  CALL Qsort_real_rec(vec, left, right)
  
END SUBROUTINE Qsort_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for integer table quicksort
SUBROUTINE Qsort_int_table(mat, LOW, HIGH, COLUMNS)
  IMPLICIT NONE
	INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: COLUMNS
  INTEGER :: n, left, right
  
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(mat,1)
  END IF

  IF(left >= right) RETURN
  
  n = SIZE(mat,2)
  ALLOCATE(tempvec(n), pv_valuevec(n))
  
  IF(PRESENT(COLUMNS)) THEN
    cols_n = SIZE(COLUMNS)
    ALLOCATE(colsvec(cols_n))
    colsvec(:) = COLUMNS(:)
  ELSE
    cols_n = n
    ALLOCATE(colsvec(cols_n))
    colsvec(:) = (/ (i, i=1,cols_n) /)
  END IF
  
  CALL Qsort_int_table_rec(mat, left, right)
  
END SUBROUTINE Qsort_int_table


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for real table quicksort
SUBROUTINE Qsort_real_table(mat, LOW, HIGH, COLUMNS)
  IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: COLUMNS
  INTEGER :: n, left, right
  
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(mat,1)
  END IF

  IF(left >= right) RETURN
  
  n = SIZE(mat,2)
  ALLOCATE(rtempvec(n), rpv_valuevec(n))
  
  IF(PRESENT(COLUMNS)) THEN
    cols_n = SIZE(COLUMNS)
    ALLOCATE(colsvec(cols_n))
    colsvec(:) = COLUMNS(:)
  ELSE
    cols_n = n
    ALLOCATE(colsvec(cols_n))
    colsvec(:) = (/ (i, i=1,cols_n) /)
  END IF
  
  CALL Qsort_real_table_rec(mat, left, right)
  
END SUBROUTINE Qsort_real_table


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for partition integer array
FUNCTION partition_int(vec, pivot, LOW, HIGH) RESULT(new_pivot)
	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: pivot
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER :: new_pivot
	INTEGER :: temp, left, right, pv_value, i
	
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(vec)
  END IF
  
  IF(left >= right) RETURN
  
  pv_value = vec(pivot)
  vec(pivot) = vec(right)
  vec(right) = pv_value

  new_pivot = left
  DO i=left,right-1
    IF(vec(i) < pv_value) THEN
      temp = vec(i)
      vec(i) = vec(new_pivot)
      vec(new_pivot) = temp
      new_pivot = new_pivot+1
    END IF
  END DO
  temp = vec(new_pivot)
  vec(new_pivot) = vec(right)
  vec(right) = temp
  
END FUNCTION partition_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for partition real array
FUNCTION partition_real(vec, pivot, LOW, HIGH) RESULT(new_pivot)
	IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: pivot
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER :: new_pivot
	INTEGER :: left, right, i
  REAL(internal_double_type) :: temp, pv_value
	
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(vec)
  END IF
  
  IF(left >= right) RETURN
  
  pv_value = vec(pivot)
  vec(pivot) = vec(right)
  vec(right) = pv_value

  new_pivot = left
  DO i=left,right-1
    IF(vec(i) < pv_value) THEN
      temp = vec(i)
      vec(i) = vec(new_pivot)
      vec(new_pivot) = temp
      new_pivot = new_pivot+1
    END IF
  END DO
  temp = vec(new_pivot)
  vec(new_pivot) = vec(right)
  vec(right) = temp
  
END FUNCTION partition_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for partition integer table
! Partition using COLUMNS, in given order. Default is to use all columns, from 1 to end
FUNCTION partition_int_table(mat, pivot, LOW, HIGH, COLUMNS) RESULT(new_pivot)
	IMPLICIT NONE
	INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: pivot
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: COLUMNS
  INTEGER :: new_pivot
	INTEGER :: left, right, i, n, k
  INTEGER, DIMENSION(:), ALLOCATABLE :: temp, pv_value, cols
  
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(mat,1)
  END IF
  
  IF(left >= right) RETURN
  
  n = SIZE(mat,2)
  ALLOCATE(temp(n), pv_value(n))
  
  IF(PRESENT(COLUMNS)) THEN
    k = SIZE(COLUMNS)
    ALLOCATE(cols(k))
    cols(:) = COLUMNS(:)
  ELSE
    k = n
    ALLOCATE(cols(k))
    cols(:) = (/ (i, i=1,n) /)
  END IF
  
  pv_value(:) = mat(pivot,:)
  mat(pivot,:) = mat(right,:)
  mat(right,:) = pv_value(:)

  new_pivot = left
  DO i=left,right-1
    IF( is_vec_smaller(mat(i,:), pv_value(:), cols) ) THEN
      temp(:) = mat(i,:)
      mat(i,:) = mat(new_pivot,:)
      mat(new_pivot,:) = temp(:)
      new_pivot = new_pivot+1
    END IF
  END DO
  temp(:) = mat(new_pivot,:)
  mat(new_pivot,:) = mat(right,:)
  mat(right,:) = temp(:)
  
END FUNCTION partition_int_table


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface for partition real table
! Partition using COLUMNS, in given order. Default is to use all columns, from 1 to end
FUNCTION partition_real_table(mat, pivot, LOW, HIGH, COLUMNS) RESULT(new_pivot)
	IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: pivot
  INTEGER, INTENT(IN), OPTIONAL :: LOW, HIGH
  INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: COLUMNS
  INTEGER :: new_pivot
	INTEGER :: left, right, i, n, k
  REAL(internal_double_type), DIMENSION(:), ALLOCATABLE :: temp, pv_value
  INTEGER, DIMENSION(:), ALLOCATABLE ::  cols
	
	IF(PRESENT(LOW)) THEN
    left = LOW
  ELSE
    left = 1
  END IF
  IF(PRESENT(HIGH)) THEN
    right = HIGH
  ELSE
    right = SIZE(mat,1)
  END IF
  
  IF(left >= right) RETURN
  
  n = SIZE(mat,2)
  ALLOCATE(temp(n), pv_value(n))
  
  IF(PRESENT(COLUMNS)) THEN
    k = SIZE(COLUMNS)
    ALLOCATE(cols(k))
    cols(:) = COLUMNS(:)
  ELSE
    k = n
    ALLOCATE(cols(k))
    cols(:) = (/ (i, i=1,n) /)
  END IF
  
  pv_value(:) = mat(pivot,:)
  mat(pivot,:) = mat(right,:)
  mat(right,:) = pv_value(:)

  new_pivot = left
  DO i=left,right-1
    IF( is_vec_smaller(mat(i,:), pv_value(:), cols) ) THEN
      temp(:) = mat(i,:)
      mat(i,:) = mat(new_pivot,:)
      mat(new_pivot,:) = temp(:)
      new_pivot = new_pivot+1
    END IF
  END DO
  temp(:) = mat(new_pivot,:)
  mat(new_pivot,:) = mat(right,:)
  mat(right,:) = temp(:)
  
END FUNCTION partition_real_table


! Internal, private routines


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation, do not call with size 1 array
RECURSIVE SUBROUTINE Qsort_int_rec(vec, left, right)
  IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: n, pivot
  
  n = right-left+1
  IF(MOD(n,2) == 0) THEN
    pivot = left-1 + n/2
  ELSE
    pivot = left-1 + (n+1)/2
  END IF
  
  pivot = partition_int_privat(vec, pivot, left, right)
   
  IF(left < pivot-1) CALL Qsort_int_rec(vec, left, pivot-1)
  IF(right > pivot+1) CALL Qsort_int_rec(vec, pivot+1, right)

END SUBROUTINE Qsort_int_rec


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation, do not call with size 1 array
RECURSIVE SUBROUTINE Qsort_real_rec(vec, left, right)
  IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: n, pivot
  
  n = right-left+1
  IF(MOD(n,2) == 0) THEN
    pivot = left-1 + n/2
  ELSE
    pivot = left-1 + (n+1)/2
  END IF
  
  pivot = partition_real_privat(vec, pivot, left, right)
   
  IF(left < pivot-1) CALL Qsort_real_rec(vec, left, pivot-1)
  IF(right > pivot+1) CALL Qsort_real_rec(vec, pivot+1, right)

END SUBROUTINE Qsort_real_rec


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation, do not call with size 1 array
RECURSIVE SUBROUTINE Qsort_int_table_rec(mat, left, right)
  IMPLICIT NONE
	INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: n, pivot
  
  n = right-left+1
  IF(MOD(n,2) == 0) THEN
    pivot = left-1 + n/2
  ELSE
    pivot = left-1 + (n+1)/2
  END IF
  
  pivot = partition_int_table_privat(mat, pivot, left, right)
   
  IF(left < pivot-1) CALL Qsort_int_table_rec(mat, left, pivot-1)
  IF(right > pivot+1) CALL Qsort_int_table_rec(mat, pivot+1, right)

END SUBROUTINE Qsort_int_table_rec


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation, do not call with size 1 array
RECURSIVE SUBROUTINE Qsort_real_table_rec(mat, left, right)
  IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: n, pivot
  
  n = right-left+1
  IF(MOD(n,2) == 0) THEN
    pivot = left-1 + n/2
  ELSE
    pivot = left-1 + (n+1)/2
  END IF
  
  pivot = partition_real_table_privat(mat, pivot, left, right)
   
  IF(left < pivot-1) CALL Qsort_real_table_rec(mat, left, pivot-1)
  IF(right > pivot+1) CALL Qsort_real_table_rec(mat, pivot+1, right)

END SUBROUTINE Qsort_real_table_rec


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation of partition integer array,
! do not call with left >= right
FUNCTION partition_int_privat(vec, pivot, left, right) RESULT(new_pivot)
	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: pivot, left, right
  INTEGER :: new_pivot
  
  pv_value = vec(pivot)
  vec(pivot) = vec(right)
  vec(right) = pv_value

  new_pivot = left
  DO i=left,right-1
    IF(vec(i) < pv_value) THEN
      temp = vec(i)
      vec(i) = vec(new_pivot)
      vec(new_pivot) = temp
      new_pivot = new_pivot+1
    END IF
  END DO
  temp = vec(new_pivot)
  vec(new_pivot) = vec(right)
  vec(right) = temp

END FUNCTION partition_int_privat


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation of partition real array,
! do not call with left >= right
FUNCTION partition_real_privat(vec, pivot, left, right) RESULT(new_pivot)
	IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:), INTENT(INOUT) :: vec
  INTEGER, INTENT(IN) :: left, right
  INTEGER :: new_pivot, pivot
  
  rpv_value = vec(pivot)
  vec(pivot) = vec(right)
  vec(right) = rpv_value

  new_pivot = left
  DO i=left,right-1
    IF(vec(i) < rpv_value) THEN
      rtemp = vec(i)
      vec(i) = vec(new_pivot)
      vec(new_pivot) = rtemp
      new_pivot = new_pivot+1
    END IF
  END DO
  rtemp = vec(new_pivot)
  vec(new_pivot) = vec(right)
  vec(right) = rtemp

END FUNCTION partition_real_privat


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation of partition integer table,
! do not call with left >= right or SIZE(cols) > SIZE(mat,2)
FUNCTION partition_int_table_privat(mat, pivot, left, right) RESULT(new_pivot)
	IMPLICIT NONE
	INTEGER, DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: pivot, left, right
  INTEGER :: new_pivot
  
  pv_valuevec(:) = mat(pivot,:)
  mat(pivot,:) = mat(right,:)
  mat(right,:) = pv_valuevec(:)

  new_pivot = left
  DO i=left,right-1
    IF( is_vec_smaller_int(mat(i,:), pv_valuevec(:), colsvec) ) THEN
      tempvec(:) = mat(i,:)
      mat(i,:) = mat(new_pivot,:)
      mat(new_pivot,:) = tempvec(:)
      new_pivot = new_pivot+1
    END IF
  END DO
  tempvec(:) = mat(new_pivot,:)
  mat(new_pivot,:) = mat(right,:)
  mat(right,:) = tempvec(:)

END FUNCTION partition_int_table_privat


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal implementation of partition real table,
! do not call with left >= right or SIZE(cols) > SIZE(mat,2)
FUNCTION partition_real_table_privat(mat, pivot, left, right) RESULT(new_pivot)
	IMPLICIT NONE
	REAL(internal_double_type), DIMENSION(:,:), INTENT(INOUT) :: mat
  INTEGER, INTENT(IN) :: pivot, left, right
  INTEGER :: new_pivot
  
  rpv_valuevec(:) = mat(pivot,:)
  mat(pivot,:) = mat(right,:)
  mat(right,:) = rpv_valuevec(:)

  new_pivot = left
  DO i=left,right-1
    IF( is_vec_smaller_real(mat(i,:), rpv_valuevec(:), colsvec) ) THEN
      rtempvec(:) = mat(i,:)
      mat(i,:) = mat(new_pivot,:)
      mat(new_pivot,:) = rtempvec(:)
      new_pivot = new_pivot+1
    END IF
  END DO
  rtempvec(:) = mat(new_pivot,:)
  mat(new_pivot,:) = mat(right,:)
  mat(right,:) = rtempvec(:)

END FUNCTION partition_real_table_privat


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal interface function
! checks if array 1 is smaller than or equal to array 2, elementwise
! Make sure that SIZE(vec1) and SIZE(vec2) > SIZE(cols)
FUNCTION is_vec_smaller_int(vec1, vec2, cols) RESULT(smaller)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(IN) :: vec1, vec2, cols
  LOGICAL :: smaller
  INTEGER :: i, k, diff

  k = SIZE(cols)
  smaller = .TRUE.
  DO i=1,k
    diff = vec2(cols(i)) - vec1(cols(i))
    IF(diff > 0) RETURN
    IF(diff == 0) CYCLE
    smaller = .FALSE.
    RETURN
  END DO
  
END FUNCTION is_vec_smaller_int


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal interface function
! checks if array 1 is smaller than or equal to array 2, elementwise
! Make sure that SIZE(vec1) and SIZE(vec2) > SIZE(cols)
FUNCTION is_vec_smaller_real(vec1, vec2, cols) RESULT(smaller)
  IMPLICIT NONE
  REAL(internal_double_type), DIMENSION(:), INTENT(IN) :: vec1, vec2
  INTEGER, DIMENSION(:), INTENT(IN) :: cols
  LOGICAL :: smaller
  INTEGER :: i, k
  REAL(internal_double_type) :: diff

  k = SIZE(cols)
  smaller = .TRUE.
  DO i=1,k
    diff = vec2(cols(i)) - vec1(cols(i))
    IF(diff > 0) RETURN
    IF(diff == 0.0_internal_double_type) CYCLE
    smaller = .FALSE.
    RETURN
  END DO
  
END FUNCTION is_vec_smaller_real

END MODULE AP_quicksort