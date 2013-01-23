!! Module for collected general-use math and related algorithms
!! -
!! Antti Penttilä
!! 2012
!! Department of Physics, University of Helsinki
MODULE AP_utils

  !! Needs PaulV Type_Kinds module to define 'FP_KIND' precision type
	USE Type_Kinds
  !! Needs PaulV Fundamental_Constants module for various constants
	USE Fundamental_Constants
  !! Needs AP_quicksort for quicksort implementation
  USE AP_quicksort
  
	IMPLICIT NONE

! Comment out with gfortran
!	EXTERNAL :: dpotrf, dpotri, dgetrf, dgetri
	
	PUBLIC

  !! Public interface for a subroutine, writes string to given unit or
  !! to a list of units. Make sure that units are open for writing.
  !! Call write_to(ch,str), where
  !! str is the string that needs to be written, and
  !! ch is either one unit number (integer), or a list of
  !! unit numbers. Use number 0 to write to std_out.
	INTERFACE write_to
		MODULE PROCEDURE write_to_one_ch, write_to_ch_list
	END INTERFACE write_to
	
  !! Public interface for a subroutine that randomly sorts
  !! the elements in a real or integer array. Call with
  !! random_sort(array), the array is returned randomly mixed.
	INTERFACE random_sort
		MODULE PROCEDURE random_sort_real, random_sort_int
	END INTERFACE random_sort

  !! Public interface for a function that will test if
  !! given integer/real array is sorted in ascending order.
  !! Call as 'sorted = is_sorted(array)', the resultint logical value
  !! 'sorted' will have a value .TRUE. is the array is sorted.
  INTERFACE is_sorted
    MODULE PROCEDURE is_sorted_int, is_sorted_real
  END INTERFACE is_sorted

  LOGICAL :: rng_inited = .FALSE.
  
  !! Combined data type for van der Corput series
	TYPE corput_series
    !! base of the series
  	INTEGER :: base
    !! sigma-parameters
  	INTEGER, DIMENSION(:), ALLOCATABLE :: sigma
    !! values of the series
  	REAL(FP_KIND), DIMENSION(:), ALLOCATABLE :: vdc
    !! current length of the series
    INTEGER :: cur_k
  	INTEGER :: alloc_n
	END TYPE corput_series

	! Private stuff
	
	INTEGER :: tictoc1 = 0, tictoc2 = 0, tictoc3 = 0

  ! For Mersenne twister 
  ! Default seed
  INTEGER, PARAMETER :: defaultsd = 4357
  ! Period parameters
  INTEGER, PARAMETER :: N = 624, N1 = N + 1
  ! the array for the state vector
  INTEGER, SAVE, DIMENSION(0:N-1) :: mt
  INTEGER, SAVE :: mti = N1
  ! Overload procedures for saving and getting mt state
  INTERFACE mtsave
    MODULE PROCEDURE mtsavef
    MODULE PROCEDURE mtsaveu
  END INTERFACE
  INTERFACE mtget
    MODULE PROCEDURE mtgetf
    MODULE PROCEDURE mtgetu
  END INTERFACE

	! What is private
	PRIVATE :: defaultsd, N, N1, mt, mti, mtsave, mtget, tictoc1, tictoc2, tictoc3
 
CONTAINS

!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Start to record time or reset the timer
SUBROUTINE tic()
	IMPLICIT NONE
	
	CALL SYSTEM_CLOCK(tictoc1, tictoc2, tictoc3)

END SUBROUTINE tic


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get the time that has elapsed from calling tic()
!! in seconds. Be sure to call tic() at least once before
!! toc() is used. If optional argument LAP is .TRUE., will
!! not reset the timer. By default will reset the timer.
FUNCTION toc(LAP) RESULT(time)
	IMPLICIT NONE
  !! Optional input
  LOGICAL, OPTIONAL :: LAP
  !! Output, elapsed time in seconds
	REAL(FP_KIND) :: time
  INTEGER :: t
	
	CALL SYSTEM_CLOCK(t)
	IF(t < tictoc1) THEN
    t = t + tictoc3-tictoc1
  ELSE
    t = t - tictoc1
  END IF
  time = (1.0_FP_KIND * t)/tictoc2
  
  IF(PRESENT(LAP)) THEN
    IF(.NOT. LAP) tictoc1 = t
  END IF

END FUNCTION toc


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Is sorted, integer?
FUNCTION is_sorted_int(vec) RESULT(sorted)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(IN) :: vec
  LOGICAL :: sorted
  INTEGER :: i, n
  
  n = SIZE(vec)
  sorted = .TRUE.
  DO i=1,n-1
    IF(vec(i) > vec(i+1)) THEN
      sorted = .FALSE.
      RETURN
    END IF
  END DO

END FUNCTION is_sorted_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Is sorted, real?
FUNCTION is_sorted_real(vec) RESULT(sorted)
  IMPLICIT NONE
  REAL(FP_KIND), DIMENSION(:), INTENT(IN) :: vec
  LOGICAL :: sorted
  INTEGER :: i, n
  
  n = SIZE(vec)
  sorted = .TRUE.
  DO i=1,n-1
    IF(vec(i) > vec(i+1)) THEN
      sorted = .FALSE.
      RETURN
    END IF
  END DO

END FUNCTION is_sorted_real


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Union of list of integers. Call with
!! union_integer(ivec, n, REMOVE), where
!! ivec is the list of integers, that will, on return, contain a
!! sorted union of input array.
!! When calling, argument n is the number of elements in the input array
!! from which the union will be formed. Use 0 if all teh elements are used.
!! On return, n will be the number of elements in the union, and the union elements
!! will be in places ivec(1:n).
!! If optional integer 'REMOVE' is given, that value is not included in the union
SUBROUTINE union_integer(ivec, n, REMOVE)
	IMPLICIT NONE
  !! Input - integer array / output - sorted union of elements
	INTEGER, DIMENSION(:), INTENT(INOUT) :: ivec
  !! Input - length of array to be handled or 0 for all values.
  !! Output - number of distinct elements in ivec
	INTEGER, INTENT(INOUT) :: n
  !! Optional integer value to be removed from the results
	INTEGER, OPTIONAL :: REMOVE

	INTEGER :: i, ex, nn
	
	IF(n == 0) THEN
		n = SIZE(ivec, 1)
	END IF
	
	CALL AP_Qsort(ivec, HIGH=n)
	
	! More code, less CPU
	IF(PRESENT(REMOVE)) THEN
		nn = 0
		ex = REMOVE
		DO i=1,n
			IF(ex == ivec(i) .OR. REMOVE == ivec(i)) CYCLE
			nn = nn+1
			ex = ivec(i)
			ivec(nn) = ex
		END DO
	ELSE
		nn = 1
		ex = ivec(1)
		DO i=2,n
			IF(ex == ivec(i)) CYCLE
			nn = nn+1
			ex = ivec(i)
			ivec(nn) = ex
		END DO
	END IF
	n = nn
	
END SUBROUTINE union_integer


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Norm of (vec_a - vec_b), where both are 3-dimensional
!! real vectors
FUNCTION vec3_norm_of_diff(vec_a, vec_b) RESULT(d)
	IMPLICIT NONE
  !! Input vectors
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: vec_a, vec_b
  !! Output real
	REAL(FP_KIND) :: d

	d = SQRT((vec_a(1)-vec_b(1))**2 + (vec_a(2)-vec_b(2))**2 + (vec_a(3)-vec_b(3))**2)

END FUNCTION vec3_norm_of_diff


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random sort
SUBROUTINE random_sort_real(d)
	IMPLICIT NONE
	REAL(FP_KIND), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER :: i,j,k
	REAL(FP_KIND) :: x,l
	
	k = SIZE(d,1)
	DO i=k,2,-1
		CALL RANDOM_NUMBER(x)
		j = FLOOR(x*i)+1
		l = d(i)
		d(i) = d(j)
		d(j) = l
	END DO
	
END SUBROUTINE random_sort_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random sort
SUBROUTINE random_sort_int(d)
	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(INOUT) :: d
	INTEGER :: i,j,k,l
	REAL(FP_KIND) :: x
	
	k = SIZE(d,1)
	DO i=k,2,-1
		CALL RANDOM_NUMBER(x)
		j = FLOOR(x*i)+1
		l = d(i)
		d(i) = d(j)
		d(j) = l
	END DO
	
END SUBROUTINE random_sort_int


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Adds 3-dimensional vec2 to vec1, and optionally scales to length sca
SUBROUTINE vec3_add_and_scale(vec1, vec2, sca)
	IMPLICIT NONE
  !! Input - vec1, output - added and scaled vector
	REAL(FP_KIND), DIMENSION(3), INTENT(INOUT) :: vec1
  !! Second input vector
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: vec2
  !! Optional input, resulting length of the vector.
  !! If not given, scaled to 1.
	REAL(FP_KIND), OPTIONAL :: sca
	REAL(FP_KIND) :: r
	
	vec1(1) = vec1(1)+vec2(1)
	vec1(2) = vec1(2)+vec2(2)
	vec1(3) = vec1(3)+vec2(3)
	IF(PRESENT(sca)) THEN
    r = SQRT(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
    vec1(1) = sca*vec1(1)/r
    vec1(2) = sca*vec1(2)/r
    vec1(3) = sca*vec1(3)/r
  END IF
	
END SUBROUTINE vec3_add_and_scale


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Scales 3-dimensional vec to length sca
SUBROUTINE vec3_scale(vec, sca)
	IMPLICIT NONE
  !! Input and output
	REAL(FP_KIND), DIMENSION(3), INTENT(INOUT) :: vec
  !! Input, length of vec on return
	REAL(FP_KIND) :: sca
	REAL(FP_KIND) :: r
	
	r = SQRT(vec(1)**2 + vec(2)**2 + vec(3)**2)
	vec(1) = sca*vec(1)/r
	vec(2) = sca*vec(2)/r
	vec(3) = sca*vec(3)/r
	
END SUBROUTINE vec3_scale


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!
!! Adds two 3-dimensional vectors
SUBROUTINE vec3_add(vec1, vec2)
	IMPLICIT NONE
  !! Added vector will be in vec1 on return
	REAL(FP_KIND), DIMENSION(3), INTENT(INOUT) :: vec1
  !! Input
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: vec2

	vec1(1) = vec1(1)+vec2(1)
	vec1(2) = vec1(2)+vec2(2)
	vec1(3) = vec1(3)+vec2(3)
	
END SUBROUTINE vec3_add


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3-dimensional vector norm
FUNCTION vec3_norm(vec) RESULT(x)
	IMPLICIT NONE
  !! Input
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: vec
  !! Result
  REAL(FP_KIND) :: x

	x = SQRT(vec(1)**2 + vec(2)**2 + vec(3)**2)
	
END FUNCTION vec3_norm


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculates the angle between 3-dimensional vec1 and vec2 in radians
FUNCTION vec3_angle_between(vec1, vec2) RESULT(ang)
	IMPLICIT NONE
  !! Inputs
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: vec1, vec2
  !! Result
	REAL(FP_KIND) :: ang

	ang = (vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)) / &
		SQRT((vec1(1)**2 + vec1(2)**2 + vec1(3)**2)* &
		(vec2(1)**2 + vec2(2)**2 + vec2(3)**2))
	ang = ACOS(ang)

END FUNCTION vec3_angle_between


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! On random number from normal distribution.
!! Initialize random number generator yourself or
!! call init_rng() once before calling any random number function.
FUNCTION normal_rnd(MEAN, STD) RESULT(x)
	IMPLICIT NONE
  !! Optional input arguments, defaults are 0.0 and 1.0
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: MEAN, STD
  !! Result
	REAL(FP_KIND) :: x
  REAL(FP_KIND) :: u, v
	
	CALL RANDOM_NUMBER(u)
	CALL RANDOM_NUMBER(v)
	x = SQRT(-2.0_FP_KIND*LOG(u)) * COS(2*PI*v)
	IF(PRESENT(STD)) x = x*STD
	IF(PRESENT(MEAN)) x = x+MEAN
	
END FUNCTION normal_rnd


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Vector of random numbers from normal distribution.
!! Initialize random number generator yourself or
!! call init_rng() once before calling any random number function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE normal_rnd_vec(vec, MEAN, STD)
	IMPLICIT NONE
  !! Output vector for random numbers
	REAL(FP_KIND), DIMENSION(:), INTENT(OUT) :: vec
  !! Optional input arguments, defaults are 0.0 and 1.0
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: MEAN, STD
	REAL(FP_KIND) :: m, s, u, v
	INTEGER :: n, i
	LOGICAL :: odd
	
	IF(PRESENT(STD)) THEN
		s = STD
	ELSE
		s = 1.0_FP_KIND
	END IF
	IF(PRESENT(MEAN)) THEN
		m = MEAN
	ELSE
		m = 0.0_FP_KIND
	END IF
	
	n = SIZE(vec, 1)

	DO i=1,n-1,2	
		CALL RANDOM_NUMBER(u)
		CALL RANDOM_NUMBER(v)
		u = SQRT(-2.0_FP_KIND*LOG(u))
		v = 2*PI*v
		vec(i) = s * u * COS(v) + m
		vec(i+1) = s * u * SIN(v) + m
	END DO
	IF(is_odd(n)) THEN
		vec(i) = s * u * SIN(v) + m
	END IF
	
END SUBROUTINE normal_rnd_vec


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! One random number from log-normal distribution.
!! Two optional parametrizations; with MEAN and STD (standard version) or
!! with REP_MEAN and REP_STD where these will be the mean and std
!! of the resulting log-normal distribution.
!! Initialize random number generator yourself or
!! call init_rng() once before calling any random number function.
FUNCTION log_normal_rnd(MEAN, STD, REP_MEAN, REP_STD) RESULT(x)
	IMPLICIT NONE
  !! Optional inputs. Default is standard parametrization with underlying
  !! normal distribution having mean 0 and std 1. If REP_MEAN is given
  !! the reparametrized version will be used, where default REP_STD = 1
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: MEAN, STD, REP_MEAN, REP_STD
  !! Result
	REAL(FP_KIND) :: x
  REAL(FP_KIND) :: m, s, y
	
	! Re-parametrized version
	IF(PRESENT(REP_MEAN)) THEN
		m = REP_MEAN
		IF(PRESENT(REP_STD)) THEN
			s = REP_STD
		ELSE
			s = 1.0_FP_KIND
		END IF
		y = LOG(s**2/m**2 + 1.0_FP_KIND)
		m = LOG(m**2/SQRT(m**2 + s**2))
		s = SQRT(y)
	! Standard version
	ELSE IF(PRESENT(MEAN)) THEN
		m = MEAN
		IF(PRESENT(STD)) THEN
			s = STD
		ELSE
			s = 1.0_FP_KIND
		END IF	
	ELSE
		m = 0.0_FP_KIND
		IF(PRESENT(STD)) THEN
			s = STD
		ELSE
			s = 1.0_FP_KIND
		END IF	
	END IF
	
	x = EXP(normal_rnd(m,s))
	
END FUNCTION log_normal_rnd


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Vector of random numbers from log-normal distribution.
!! Two optional parametrizations; with MEAN and STD (standard version) or
!! with REP_MEAN and REP_STD where these will be the mean and std
!! of the resulting log-normal distribution.
!! Initialize random number generator yourself or
!! call init_rng() once before calling any random number function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE log_normal_rnd_vec(vec, MEAN, STD, REP_MEAN, REP_STD)
	IMPLICIT NONE
  !! Output vector of random numbers
	REAL(FP_KIND), DIMENSION(:), INTENT(OUT) :: vec
  !! Optional inputs. Default is standard parametrization with underlying
  !! normal distribution having mean 0 and std 1. If REP_MEAN is given
  !! the reparametrized version will be used, where default REP_STD = 1
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: MEAN, STD, REP_MEAN, REP_STD
	REAL(FP_KIND) :: m, s, y
	INTEGER :: i, n
	
	! Re-parametrized version
	IF(PRESENT(REP_MEAN)) THEN
		m = REP_MEAN
		IF(PRESENT(REP_STD)) THEN
			s = REP_STD
		ELSE
			s = 1.0_FP_KIND
		END IF
		y = LOG(s**2/m**2 + 1.0_FP_KIND)
		m = LOG(m) - 0.5_FP_KIND*y
		s = SQRT(y)
	! Standard version
	ELSE IF(PRESENT(MEAN)) THEN
		m = MEAN
		IF(PRESENT(STD)) THEN
			s = STD
		ELSE
			s = 1.0_FP_KIND
		END IF	
	ELSE
		m = 0.0_FP_KIND
		IF(PRESENT(STD)) THEN
			s = STD
		ELSE
			s = 1.0_FP_KIND
		END IF	
	END IF
	
	CALL normal_rnd_vec(vec,m,s)
	n = SIZE(vec, 1)
	DO i=1,n
		vec(i) = EXP(vec(i))
	END DO
	
END SUBROUTINE log_normal_rnd_vec


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Linear interpolation. Interpolate value y at x, based on the
!! linear interpolation between nearest points in table tab
FUNCTION linear_interpolation(x, tab, ind, STAT, PREV_I) RESULT(y)
	IMPLICIT NONE
  !! Input, where to interpolate
	REAL(FP_KIND), INTENT(IN) :: x
  !! Input data matrix. Must have (at least) two columns where the 
  !! x and y -values of the function are stored, and rows must be
  !! sorted in ascending order in x.
	REAL(FP_KIND), DIMENSION(:,:), INTENT(IN) :: tab
  !! Input, column numbers of x and y -values in data matrix.
	INTEGER, DIMENSION(2), INTENT(IN) :: ind
  !! Optional status code, output. If the x is below the value range
  !! in data, -1 will be returned, and -2 if x is above the range.
  !! In normal operation 0 will be returned.
	INTEGER, OPTIONAL, INTENT(OUT) :: STAT
  !! Optional input. If you know beforehand the row numbers
  !! that enclose the given x value, you can give the lower limit row number
  !! to fasten the computation.
	INTEGER, OPTIONAL, INTENT(INOUT) :: PREV_I
  !! Result
	REAL(FP_KIND) :: y
	REAL(FP_KIND) :: c
	INTEGER :: n, i, j, ind1, ind2
	
	ind1 = ind(1)
	ind2 = ind(2)
	
	IF(PRESENT(PREV_I)) THEN
		IF(PREV_I > 0) THEN
			i = PREV_I
			j = i+1
			GOTO 99
		END IF
	END IF

	n = SIZE(tab,1)
	! Check that value is inside range
	IF(x < tab(1,ind1)) THEN
		IF(PRESENT(STAT)) STAT = -1
		RETURN
	ELSE IF(x > tab(n,ind1)) THEN
		IF(PRESENT(STAT)) STAT = -2
		RETURN
	END IF
	
	i = FLOOR(n/2.0_FP_KIND)
	IF(tab(i,ind1) > x) THEN
	! Search downwards
		i = i-1
		DO WHILE(tab(i,ind1) > x)
			i = i-1
		END DO
		j = i+1
	ELSE
	! Search upwards
		i = i+1
		DO WHILE(tab(i,ind1) < x)
			i = i+1
		END DO
		j = i
		i = i-1
	END IF
	
99 CONTINUE
	c = (x-tab(i,ind1))/(tab(j,ind1)-tab(i,ind1))
	y = tab(i,ind2) + c*(tab(j,ind2)-tab(i,ind2))
	
	IF(PRESENT(STAT)) STAT = 0
	IF(PRESENT(PREV_I)) PREV_I = i
	
END FUNCTION linear_interpolation


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inits the random number generator from system clock.
! Call once before using the random number generator.
! Works at least with Absoft Fortran, comment out for others
! SUBROUTINE init_rng()
	! IMPLICIT NONE
	! INTEGER :: s_size
	! INTEGER, DIMENSION(8) :: t
	! INTEGER, DIMENSION(1) :: seed
	! REAL(Double) :: x

	! CALL RANDOM_SEED(SIZE=s_size)
	! IF(s_size /= 1) THEN
		! WRITE(*,*) "AnttiUtils|init_rng: RNG needs multiple seeds"
		! STOP
	! END IF
	! CALL DATE_AND_TIME(VALUES=t)
	! seed = 100*t(7) + t(8)/10
	! CALL RANDOM_SEED(PUT=seed)
	! CALL RANDOM_NUMBER(x) ! First number is not good
  
 ! rng_inited = .TRUE.

! END SUBROUTINE init_rng


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Inits the random number generator from system clock.
!! Call once before using the random number generator.
!! Works at least with gfortran
SUBROUTINE init_rng()
	IMPLICIT NONE
	INTEGER :: ssize
	INTEGER, DIMENSION(8) :: t
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed
	REAL(FP_KIND) :: x

	CALL RANDOM_SEED(SIZE=ssize)
	ALLOCATE(seed(ssize))
	CALL DATE_AND_TIME(VALUES=t)
	seed = 100*t(7) + t(8)/10
	CALL RANDOM_SEED(PUT=seed)
	CALL RANDOM_NUMBER(x) ! First number is not good

  rng_inited = .TRUE.

END SUBROUTINE init_rng


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Real-valued random number in a given interval.
FUNCTION give_rnr(LOW, UP) RESULT(rn)
	IMPLICIT NONE
  !! Optional inputs, lower and upper limits for the
  !! random number. Default is [0,1[. Result is always
  !! smaller than the upper limit.
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: LOW, UP
  !! Result
	REAL(FP_KIND) :: rn
	REAL(FP_KIND) :: ll, ul, lw, x
	
	ll = 0.0_FP_KIND
	ul = 1.0_FP_KIND
	IF(PRESENT(LOW)) ll = LOW
	IF(PRESENT(UP)) ul = UP

	CALL RANDOM_NUMBER(x)
	lw = ul-ll
	rn = ll + lw*x
	
END FUNCTION give_rnr


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Integer-valued random number in a given interval.
FUNCTION give_rnr_int(LOW, UP) RESULT(rn)
	IMPLICIT NONE
  !! Optional inputs, lower and upper limits for the
  !! random number. Default is [0,1]. Both lower and upper
  !! limits are possible to attain.
	INTEGER, INTENT(IN), OPTIONAL :: LOW, UP
  !! Result
	INTEGER :: rn
	INTEGER :: ll, ul, lw
	REAL(FP_KIND) :: x

	ll = 0
	ul = 2
	IF(PRESENT(LOW)) ll = LOW
	IF(PRESENT(UP)) ul = UP+1

	CALL RANDOM_NUMBER(x)
	lw = ul-ll
	rn = ll + FLOOR(lw*x)
	
END FUNCTION give_rnr_int


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Give rotation matrix R, which rotates
!! vector (0,0,1) to (a,b,c) when applied
!! as (R . v). The (a,b,c) and v should have unit-lenght
SUBROUTINE give_rotation_matrix(rvec, rmat)
	IMPLICIT NONE
  !! Input, where to rotate
	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: rvec
  !! Output, rotation matrix
	REAL(FP_KIND), DIMENSION(3,3), INTENT(INOUT) :: rmat
	
	rmat(1,1) = (rvec(2)**2 + rvec(1)**2*rvec(3))/(rvec(1)**2 + rvec(2)**2)
	rmat(1,2) = (rvec(1)*rvec(2)*(rvec(3)-1))/(rvec(1)**2 + rvec(2)**2)
	rmat(1,3) = rvec(1)
	rmat(2,1) = (rvec(1)*rvec(2)*(rvec(3)-1))/(rvec(1)**2 + rvec(2)**2)
	rmat(2,2) = (rvec(1)**2 + rvec(2)**2*rvec(3))/(rvec(1)**2 + rvec(2)**2)
	rmat(2,3) = rvec(2)
	rmat(3,1) = -rvec(1)
	rmat(3,2) = -rvec(2)
	rmat(3,3) = rvec(3)

END SUBROUTINE give_rotation_matrix


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Random direction / vector on unit sphere
SUBROUTINE vec3_rnd_dir(dir_v)
	IMPLICIT NONE
  !! Output, random vector on unit sphere
	REAL(FP_KIND), DIMENSION(3), INTENT(OUT) :: dir_v
	REAL(FP_KIND) :: theta,phi,x
	
	x = give_rnr(-1.0_FP_KIND, 1.0_FP_KIND)
	theta = ACOS(x)
	phi = give_rnr(0.0_FP_KIND, 2*PI)
	
	dir_v(1) = COS(phi)*SIN(theta)
	dir_v(2) = SIN(phi)*SIN(theta)
	dir_v(3) = COS(theta)
	
END SUBROUTINE vec3_rnd_dir


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Random direction and length from Henyey-Greenstsein
!! direction distribution and reparametrized log-normal
!! radius direction. The radius distribution is a reparametrized log-normal
!! where tau is both the mean and the standard deviation of the distribution.
!! The tau is computed as sca*exp(1-round^(1/D)). Parameters round and D
!! are related to the simulated annealing temperature factor.
SUBROUTINE rnd_HGR(res_v, g, round, D, sca)
	IMPLICIT NONE
  !! Output, vector with direction from Henyey-Greenstsein and
  !! length from lognormal distributions
	REAL(FP_KIND), DIMENSION(3), INTENT(OUT) :: res_v
  !! Input, asymmetry factor, parameter for H-G
	REAL(FP_KIND), INTENT(IN) :: g
  !! Input, round parameter for simulated annealing scaling. Use 1 for no scaling.
	INTEGER, INTENT(IN) :: round
  !! Input, D parameter for simulated annealing scaling.
	INTEGER, INTENT(IN) :: D
  !! Optional extra scaling factor for radius.
	REAL(FP_KIND), INTENT(IN), OPTIONAL :: sca
	REAL(FP_KIND) :: sca_c, stheta, ctheta, phi, cphi, sphi, tau, r, x
	
	IF(PRESENT(sca)) THEN
		sca_c = sca
	ELSE
		sca_c = 1.0_FP_KIND
	END IF
	
	x = give_rnr()
	IF(g == 0) THEN
		ctheta = 2*x-1
	ELSE
		ctheta = (1 + g**2 - ((1-g**2)/(1+g*(2*x-1)))**2) / (2*g)
	END IF
	stheta = SQRT(1-ctheta**2)
	phi = give_rnr(0.0_FP_KIND, 2*PI)
	cphi = COS(phi)
	sphi = SIN(phi)
	tau = sca_c*EXP(1 - round**(1.0_FP_KIND/D))
	r = log_normal_rnd(REP_MEAN=tau, REP_STD=tau)
	
	res_v(1) = r*stheta*cphi
	res_v(2) = r*stheta*sphi
	res_v(3) = r*ctheta
	
END SUBROUTINE rnd_HGR


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Random rotation for 3-dimensional vector
!! Rotation angle is from exponential distribution in range 0-180 deg
SUBROUTINE vec3_rnd_rotation(dir_v, p)
	IMPLICIT NONE
  !! Input - vector to be rotated, output - rotated vector
	REAL(FP_KIND), DIMENSION(3), INTENT(INOUT) :: dir_v
  !! Input, parameter for exponential distribution
	REAL(FP_KIND), INTENT(IN) :: p

	REAL(FP_KIND) :: theta,phi
	REAL(FP_KIND), DIMENSION(3) :: ov,rv
	
	! First random theta from exponential
	theta = -p * LOG(give_rnr())
	IF(theta > 180) theta = 179.9_FP_KIND
	! Then random phi
	phi = give_rnr(0.0_FP_KIND, 2*PI)
	
	rv(1) = COS(phi)*SIN(theta)
	rv(2) = SIN(phi)*SIN(theta)
	rv(3) = COS(theta)
	ov(:) = dir_v(:)
	
	! Rotation
	dir_v(1) = (ov(2)**2*rv(1) + ov(1)**2*ov(3)*rv(1) + ov(1)*ov(2)*(-1 + ov(3))*rv(2)) / &
		(ov(1)**2 + ov(2)**2) + ov(1)*rv(3)
	dir_v(2) = (ov(1)*ov(2)*(-1 + ov(3))*rv(1) + ov(1)**2*rv(2) + ov(2)**2*ov(3)*rv(2)) / &
		(ov(1)**2 + ov(2)**2) + ov(2)*rv(3)
	dir_v(3) = -(ov(1)*rv(1)) - ov(2)*rv(2) + ov(3)*rv(3)
	
END SUBROUTINE vec3_rnd_rotation


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Print matrix row-vise to a given chanel.
!! Use channel 0 for std_out, and make sure
!! that the channel is open for writing.
SUBROUTINE print_matrix(m, CHANNEL)
	IMPLICIT NONE
  !! Input, real matrix to be printed
	REAL(FP_KIND), DIMENSION(:,:), INTENT(IN) :: m
  !! Optional, channel where to print. 0 (default) for std_out
	INTEGER, OPTIONAL, INTENT(IN) :: CHANNEL
	INTEGER :: lch, i, ms
	
	ms = SIZE(m,1)
	IF(PRESENT(CHANNEL)) THEN
		lch = CHANNEL
	ELSE
		lch = 0
	END IF
	
	IF(lch == 0) THEN
		DO i=1,ms
			WRITE(*,*) m(i,:)
		END DO
	ELSE
		DO i=1,ms
			WRITE(lch,*) m(i,:)
		END DO
	END IF
	
END SUBROUTINE print_matrix


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writes to one channel
SUBROUTINE write_to_one_ch(ch, str)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: ch
	CHARACTER(LEN=*), INTENT(IN) :: str
	
	IF(ch == 0) THEN
		WRITE(*,*) str
	ELSE IF(ch > 0) THEN
		WRITE(ch,*) str
	END IF
	
END SUBROUTINE write_to_one_ch


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writes to multiple channels
SUBROUTINE write_to_ch_list(chl, str)
	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(IN) :: chl
	CHARACTER(LEN=*), INTENT(IN) :: str
	INTEGER :: i, n
	
	n = SIZE(chl, 1)
	DO i=1,n
		IF(chl(i) == 0) THEN
			WRITE(*,*) str
		ELSE IF(chl(i) > 0) THEN
			WRITE(chl(i),*) str
		END IF
	END DO
	
END SUBROUTINE write_to_ch_list


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!
!! is integer odd?
FUNCTION is_odd(n) RESULT(odd)
  IMPLICIT NONE
  !! Input
  INTEGER, INTENT(IN) :: n
  !! Result, .TRUE. if n is odd
  LOGICAL :: odd
  
	IF(MOD(n,2) == 0) THEN
		odd = .FALSE.
	ELSE
		odd = .TRUE.
	END IF

END FUNCTION is_odd


! Below are subroutines etc. for Mersenne twister RNG and
! some utilities, plus van der Corput chain
!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
!! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!! Fortran version rewritten as an F90 module and mt state saving and getting
!! subroutines added by Richard Woloshyn. June 30, 1999
!! -
!! Inits van der Corput chain
FUNCTION init_corput(base, RANDOMINIT, SER_LEN) RESULT(vdcs)
  IMPLICIT NONE
  !! Input, base
  INTEGER, INTENT(IN) :: base
  !! Optional input, set .TRUE. if you need to init
  !! the random number generator
  LOGICAL, INTENT(IN), OPTIONAL :: RANDOMINIT
  !! Optional length of series, default is 50
  INTEGER, INTENT(IN), OPTIONAL :: SER_LEN
  !! Result, van der Corput series
  TYPE(corput_series) :: vdcs
  
  INTEGER, PARAMETER :: def_ser_len = 50
  INTEGER :: i, sl, astat
  
  vdcs%base = base
  ALLOCATE(vdcs%sigma(0:base-1))
  DO i=0,base-1
    vdcs%sigma(i) = i
  END DO
  IF(PRESENT(RANDOMINIT)) THEN
    IF(RANDOMINIT) THEN
      IF(.NOT. rng_inited) THEN
        CALL init_rng()
      END IF
      CALL random_sort(vdcs%sigma)
    END IF
  END IF
 
 IF(PRESENT(SER_LEN)) THEN
    sl = SER_LEN
  ELSE
    sl = def_ser_len
  END IF
  
  ALLOCATE(vdcs%vdc(0:sl), STAT=astat)
  IF(astat /= 0) THEN
    WRITE(*,*) "init_corput: error in allocation."
    STOP
  END IF
  vdcs%alloc_n = sl
  
  vdcs%vdc(0) = REAL(vdcs%sigma(0),FP_KIND) / (vdcs%base-1)
  vdcs%cur_k = 0
  
END FUNCTION init_corput

!  PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Computes van der Corput series up to order k.
!! series (vdcs) must be inited first.
SUBROUTINE corput(vdcs, k)
  IMPLICIT NONE
  !! Input - exisiting series object
  !! Output - same series but with k elemenents computed.
  TYPE(corput_series), INTENT(INOUT) :: vdcs
  !! Input, order up to which compute elements
  INTEGER, INTENT(IN) :: k
  
  INTEGER :: astat, i, new_alloc, n, r
  REAL(FP_KIND), DIMENSION(:), ALLOCATABLE :: temp_vdc

  IF(k < vdcs%cur_k) RETURN
    
  IF(k > vdcs%alloc_n) THEN
    ALLOCATE(temp_vdc(0:vdcs%alloc_n), STAT=astat)
    IF(astat /= 0) THEN
      WRITE(*,*) "corput: allocation error(1)"
      STOP
    END IF
    temp_vdc(0:vdcs%alloc_n) = vdcs%vdc(0:vdcs%alloc_n)
    DEALLOCATE(vdcs%vdc)
    new_alloc = 2*k
    ALLOCATE(vdcs%vdc(0:new_alloc), STAT=astat)
    IF(astat /= 0) THEN
      WRITE(*,*) "corput: allocation error(2)"
      STOP
    END IF
    vdcs%vdc(0:vdcs%alloc_n) = temp_vdc(0:vdcs%alloc_n)
    vdcs%alloc_n = new_alloc
  END IF
  
  DO i=vdcs%cur_k+1, k
    n = INT(i/vdcs%base)
    r = MOD(i, vdcs%base)
    vdcs%vdc(i) = (vdcs%vdc(n)+vdcs%sigma(r))/vdcs%base
  END DO
  vdcs%cur_k = k
 
END SUBROUTINE corput

! PRIVATE!!!!!!!!!!!!!!!!!!!!
! Initialization subroutine
SUBROUTINE sgrnd(seed)
  IMPLICIT NONE
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
  INTEGER, INTENT(IN), OPTIONAL :: seed
  
  ! My additions for system clock init
  INTEGER, DIMENSION(8) :: t
  INTEGER :: lseed
  
  IF(.NOT. PRESENT(seed)) THEN
    CALL DATE_AND_TIME(VALUES=t)
    lseed = 100*t(7) + t(8)/10
  ELSE
    lseed = seed
  END IF

  mt(0) = IAND(lseed,-1)
  DO mti=1,N-1
    mt(mti) = IAND(69069 * mt(mti-1),-1)
  END DO

END SUBROUTINE sgrnd


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mersene-twister random number generator. Call sgrnd first.
FUNCTION grnd() RESULT(res)
  IMPLICIT NONE

  REAL(FP_KIND) :: res
  
! Period parameters
  INTEGER, PARAMETER :: M = 397, MATA  = -1727483681
!                                    constant vector a
  INTEGER, PARAMETER :: LMASK =  2147483647
!                                    least signIFicant r bits
  INTEGER, PARAMETER :: UMASK = -LMASK - 1
!                                    most signIFicant w-r bits
! Tempering parameters
  INTEGER, PARAMETER :: TMASKB= -1658038656, TMASKC= -272236544

  INTEGER, DIMENSION(0:1), SAVE ::  mag01 = (/ 0, MATA /)
! mag01(x) = x * MATA for x=0,1

  INTEGER :: y, kk

  IF(mti >= N) THEN
!                       generate N words at one time
    IF(mti == N+1) THEN
!                            IF sgrnd() has not been called,
      CALL sgrnd( defaultsd )
!                              a default initial seed is used
    END IF

    DO kk=0,N-M-1
      y=ior(IAND(mt(kk),UMASK),IAND(mt(kk+1),LMASK))
      mt(kk)=IEOR(IEOR(mt(kk+M),ISHFT(y,-1)),mag01(IAND(y,1)))
    END DO
    DO kk=N-M,N-2
      y=ior(IAND(mt(kk),UMASK),IAND(mt(kk+1),LMASK))
      mt(kk)=IEOR(IEOR(mt(kk+(M-N)),ISHFT(y,-1)),mag01(IAND(y,1)))
    END DO
    y=ior(IAND(mt(N-1),UMASK),IAND(mt(0),LMASK))
    mt(N-1)=ieor(IEOR(mt(M-1),ISHFT(y,-1)),mag01(IAND(y,1)))
    mti = 0
  END IF

  y=mt(mti)
  mti = mti + 1 
  y=IEOR(y,ISHFT(y,-11))
  y=IEOR(y,IAND(ISHFT(y,7),TMASKB))
  y=IEOR(y,IAND(ISHFT(y,15),TMASKC))
  y=IEOR(y,ISHFT(y,-18))

  IF(y < 0) THEN
    res=(DBLE(y)+2.0_FP_KIND**32)/(2.0_FP_KIND**32-1.0_FP_KIND)
  ELSE
    res=DBLE(y)/(2.0_FP_KIND**32-1.0_FP_KIND)
  END IF

end FUNCTION grnd

! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! State saving SUBROUTINEs.
! Usage:  call mtsave( file_name, format_character )
!    or   call mtsave( unit_number, format_character )
! where   format_character = 'u' or 'U' will save in unformatted form, otherwise
!         state information will be written in formatted form.
SUBROUTINE mtsavef( fname, forma )

!NOTE: This SUBROUTINE APPENDS to the end of the file "fname".
  IMPLICIT NONE
  character(*), intent(in) :: fname
  character, intent(in)    :: forma

  select case (forma)
    case('u','U')
     open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED', &
          position='APPEND')
     write(10)mti
     write(10)mt

    case default
     open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED', &
          position='APPEND')
     write(10,*)mti
     write(10,*)mt

  end select
  close(10)

END SUBROUTINE mtsavef

SUBROUTINE mtsaveu( unum, forma )
  IMPLICIT NONE
  INTEGER, intent(in)    :: unum
  character, intent(in)  :: forma

  select case (forma)
    case('u','U')
     write(unum)mti
     write(unum)mt

    case default
     write(unum,*)mti
     write(unum,*)mt

  end select

end SUBROUTINE mtsaveu

! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! State getting SUBROUTINEs.
! Usage:  call mtget( file_name, format_character )
!    or   call mtget( unit_number, format_character )
! where   format_character = 'u' or 'U' will read in unformatted form, otherwise
!         state information will be read in formatted form.
SUBROUTINE mtgetf( fname, forma )
  IMPLICIT NONE
  character(*), intent(in) :: fname
  character, intent(in)    :: forma

  select case (forma)
    case('u','U')
     open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
     read(10)mti
     read(10)mt

    case default
     open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
     read(10,*)mti
     read(10,*)mt

  end select
  close(10)

end SUBROUTINE mtgetf

! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mtgetu( unum, forma )
  IMPLICIT NONE
  INTEGER, intent(in)    :: unum
  character, intent(in)  :: forma

  select case (forma)
    case('u','U')
     read(unum)mti
     read(unum)mt

    case default
     read(unum,*)mti
     read(unum,*)mt

  end select

END SUBROUTINE mtgetu


END MODULE AP_utils