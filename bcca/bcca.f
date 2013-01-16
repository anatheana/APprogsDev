!! Ballistic cluster-cluster aggregation
!! -
!! Command line tool for simulating ballistic cluster-cluster aggregation
!! using homogeneous spheres
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Ballistic cluster-cluster aggregation (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2013).
!! -
!! Antti Penttilä
!! 2013
!! Department of Physics, University of Helsinki

PROGRAM bcca
  !! Double precision FP_KIND type parameter from PaulV
  USE Type_Kinds
  !! File handling routines from PaulV
  USE File_Utility
  !! String routines from PaulV
  USE String_Utility
  !! My own math- and printing-related routines
  USE AP_utils
  
  INTEGER, PARAMETER :: fp = FP_KIND, fn_max_char = 256, arg_max_char = 256
  REAL(fp), PARAMETER :: sphVolC = 4.18879020478639_fp, & ! 4/3*pi
    r = 0.5_fp, r2 = 1.0_fp ! fix diameter = 1
  CHARACTER(LEN=fn_max_char) :: def_fn_out = "bcca.geom"
  
  INTEGER :: n, i, j, k, stata, fp_out, nPow
  
  REAL(fp) :: x, y, z, enc1, enc2, volMat, volEnc, pD, encR
  REAL(fp), DIMENSION(3) :: co1, co2, cent1, cent2, dir1, dir2
  REAL(fp), DIMENSION(3,3) :: rmat
  REAL(fp), DIMENSION(:), ALLOCATABLE :: encl
  REAL(fp), DIMENSION(:,:), ALLOCATABLE :: sph, tsph
  
  CHARACTER(LEN=fn_max_char) :: fn_out

  WRITE(*,*) ""
  WRITE(*,*) "*****"
  WRITE(*,*) "Starting BCCA..."

 ! Init
  WRITE(*,*) ""
  WRITE(*,*) "initializing..."
  CALL init_rng()
  WRITE(*,*) "OK"
  
  ! Input
  WRITE(*,*) ""
  WRITE(*,*) "handling input..."
  CALL handle_input()
  WRITE(*,'(A30,I0,A11)') " building BCCA aggregate with ", n, " spheres..."
  
  ! Do BCCA
  CALL first_round()
  DO i=2,nPow
    k = 2**i
    DO j=1,n,k
      CALL aggregate(j,k)
    END DO
  END DO
  
  ! Analyzing
  x = 0.0_fp
  DO i=1,n
    y = vec3_norm(sph(i,:))
    IF(x < y) x = y
  END DO
  encR = x + r
  volEnc = sphVolC * encR**3
  volMat = n * sphVolC * r**3
  pD = volMat/volEnc
  WRITE(*,*) ""
  WRITE(*,'(A28,I0,A24,E12.6,A1)') " Ready, BCCA aggregate with ", n, " uniform spheres with r=", r, ","
  WRITE(*,'(A17,F8.4,A31,E12.6)') " Packing density ", 100*pD, " %, radius of enclosing sphere ", encR

  ! Writing output
  WRITE(*,*) ""
  WRITE(*,*) "writing output..."
  fp_out = Get_Lun()
  OPEN(fp_out, FILE=fn_out, ACTION='write', STATUS='replace', IOSTAT=stata)
  DO i=1,n
    WRITE(fp_out, '(E13.6,1X,E13.6,1X,E13.6,1X,E13.6)') r, sph(i,:)
  END DO
  CLOSE(fp_out)
  
CONTAINS


SUBROUTINE aggregate(ind, num)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind, num
  INTEGER :: li, lj, lk, ind2, ind3, old_n
  REAL(fp) :: d, max_d
  LOGICAL :: success
  
  old_n = num/2
  ind2 = ind+old_n
  ind3 = ind+num
  enc1 = encl(ind)
  enc2 = encl(ind2)
  
  success = .FALSE.
  ! Loop until hit
  hitloop: DO WHILE(.NOT. success)
    tsph(1:old_n,:) = sph(ind2:ind3-1,:)
    ! random direction
    CALL vec3_rnd_dir(dir1)
    ! Random place in unit circle in x-y-plane
    x = give_rnr(UP=2*PI)
    dir2 = (/ SIN(x), COS(x), 0.0_fp /)
    ! Rotation matrix from (0,0,1) to dir1
    CALL give_rotation_matrix(dir1, rmat)
    dir2 = MATMUL(rmat,dir2)
    ! After scaling, dir2 is inside enc1+enc2-circle with normal dir1
    ! so, place where aggregate2 is approaching from direction -dir1
    y = SQRT(give_rnr())*(enc1+enc2)
    CALL vec3_scale(dir2, y)
    ! Now move second aggregate to new center cent2 = dir2 + z*dir1
    z = SQRT((enc1+enc2)**2 - y**2)
    cent2 = dir2 + z*dir1
    DO li=1,old_n
      CALL vec3_add_and_scale(tsph(li,:), cent2)
    END DO
    
    ! Find possible contact
    max_d = HUGE(max_d)
    DO li=1,old_n
      DO lj=ind,ind2-1
        ! Distance from line
        d = SQRT((dir1(2)*(tsph(li,1) - sph(lj,1)) + dir1(1)*(-tsph(li,2) + &
          sph(lj,2)))**2 + (dir1(3)*(tsph(li,1) - sph(lj,1)) + dir1(1)*(-tsph(li,3) + &
          sph(lj,3)))**2 + (dir1(3)*(tsph(li,2) - sph(lj,2)) + dir1(2)*(-tsph(li,3) + &
          sph(lj,3)))**2)
        IF(d > r2) CYCLE
        ! Contact, possible new place and distance to move along -dir1
        y = SQRT(r2**2 - d**2)
        z = vec3_norm_of_diff(tsph(li,:), sph(lj,:))
        x = SQRT(z**2 - d**2) - y
        IF(x < max_d) max_d = x
        success = .TRUE.
      END DO
    END DO
  
  END DO hitloop

  sph(ind2:ind3-1,:) = tsph(1:old_n,:)
  ! Move to correct place, find combinded origin and move
  co1 = -max_d*dir1
  DO li=ind2,ind3-1
    CALL vec3_add_and_scale(sph(li,:), co1)
  END DO
  cent1(1) = SUM(sph(ind:ind3-1,1))/num
  cent1(2) = SUM(sph(ind:ind3-1,2))/num
  cent1(3) = SUM(sph(ind:ind3-1,3))/num
  DO li=ind,ind3-1
    CALL vec3_add_and_scale(sph(li,:), -cent1)
  END DO
  ! Find new enclosing sphere
  enc1 = 0.0_fp
  DO li=ind,ind3-1
    enc2 = vec3_norm(sph(li,:))
    IF(enc2 > enc1) enc1 = enc2
  END DO
  encl(ind) = enc1+r
  
END SUBROUTINE aggregate


!! First round of BCCA, build two-sphere aggregates
SUBROUTINE first_round()
  IMPLICIT NONE

  DO i=1,n,2
    CALL vec3_rnd_dir(co1)
    sph(i,:) = -r * co1(:)
    sph(i+1,:) = r * co1(:)
    encl(i) = r2
  END DO
  
END SUBROUTINE first_round


!! alpha-version
SUBROUTINE handle_input()
  IMPLICIT NONE
  CHARACTER(LEN=arg_max_char) :: arg
  
  i = COMMAND_ARGUMENT_COUNT()
  SELECT CASE(i)
  CASE(2)
    CALL GET_COMMAND_ARGUMENT(1, arg)
    READ(arg, *) nPow
    CALL GET_COMMAND_ARGUMENT(2, arg)
    WRITE(fn_out, '(A)') TRIM(arg)    
  CASE(1)
    CALL GET_COMMAND_ARGUMENT(1, arg)
    READ(arg, *) nPow
    WRITE(fn_out, '(A9)') "bcca.geom"
  CASE(0)
    nPow = 5
    WRITE(fn_out, '(A9)') "bcca.geom"
  CASE DEFAULT
    WRITE(*,*) "Wrong number of arguments, 0-2 allowed"
    STOP
  END SELECT
  n = 2**nPow
  
  ALLOCATE(sph(n,3), tsph(n/2,3), encl(n), STAT=stata)
  IF(stata /= 0) THEN
    WRITE(*,*) "Error in memory allocation"
    STOP
  END IF
  encl(:) = 0.0_fp
  
END SUBROUTINE handle_input


END PROGRAM bcca