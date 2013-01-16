!! Packing properties
!! -
!! Command line tool for simulating different properties and their distributions
!! in three-dimensional sphere packings, including nearest-neigbour.
!! Vessel including the spheres is a minimum-bounding-box and can be periodic.
!! 
!! Spheres can have radii or they can be points (zero volume).
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Properties of sphere packings (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2013).
!! -
!! Antti Penttilä
!! 2013
!! Department of Physics, University of Helsinki

PROGRAM PackingProperties
  !! Double precision FP_KIND type parameter from PaulV
  USE Type_Kinds
  !! File handling routines from PaulV
  USE File_Utility
  !! String routines from PaulV
  USE String_Utility
  !! My own math- and printing-related routines
	USE AP_utils
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: fp = FP_KIND, pack_report = 50, pack_max_try = 10000, &
    fn_max_char = 256
  REAL(fp), PARAMETER :: max_pD = 0.641_fp, &
    sphVolC = 4.18879020478639_fp, & ! 4/3*pi
    r = 0.5_fp ! fix diameter = 1
  
  INTEGER :: n, i, j, k, stata, fp_out
  
  REAL(fp) :: pD, x, y, z, matVol, boxVol, sphVol, box_x, box_y, box_z
  REAL(fp), DIMENSION(:), ALLOCATABLE :: near_neigh, near_cToc
  REAL(fp), DIMENSION(:,:), ALLOCATABLE :: sph
  
  CHARACTER(LEN=fn_max_char) :: fn_out, fn_geom

  WRITE(*,*) ""
  WRITE(*,*) "*****"
  WRITE(*,*) "Starting PackingProperties..."

 ! Init
  WRITE(*,*) ""
  WRITE(*,*) "initializing..."
  CALL init_rng()
  WRITE(*,*) "OK"
  
  ! Input
  WRITE(*,*) ""
  WRITE(*,*) "handling input..."
  CALL handle_input()
  IF(pD > max_pD) THEN
    WRITE(*,*) "Too high packing density requested (", pD, &
      "), must be less than ", max_pD
    RETURN
  END IF
  WRITE(*,'(A23,I0,A13,F0.3,A3,F0.3,A3,F0.3,A18,F0.0)') " input read, requested ", &
    n, " spheres in (", box_x, " x ", box_y,     " x ", box_z, ") periodic cuboid,"
  WRITE(*,'(A31,F5.3)') "  resulting packing density of ", pD
  
  ! Pack vessel
  WRITE(*,*) ""
  WRITE(*,*) "packing spheres..."
  CALL do_pack()
  
  ! Compute properties
  WRITE(*,*) ""
  WRITE(*,*) "computing packing properties..."
  CALL comp_prop()
  
  ! Writing output
  WRITE(*,*) ""
  WRITE(*,*) "writing output..."
  fp_out = Get_Lun()
  ! outproperties
  OPEN(fp_out, FILE=fn_out, ACTION='write', STATUS='replace', IOSTAT=stata)
  WRITE(fp_out, '(A4,I0)') "n = ", n
  WRITE(fp_out, '(A5,F6.4)') "pD = ", pD
  WRITE(fp_out, '(A23,E12.6)') "mean-nearest-neighbour ", SUM(near_neigh)/n
  WRITE(fp_out, '(A30,E12.6)') "mean-nearest-center-to-center ", SUM(near_cToc)/n
  WRITE(fp_out, '(A19)') "nearest-neighbours:"
  DO i=1,n
    WRITE(fp_out, '(E13.6)'), near_neigh(i)
  END DO
  WRITE(fp_out, '(A25)') "nearest-center-to-center:"
  DO i=1,n
    WRITE(fp_out, '(E13.6)'), near_cToc(i)
  END DO
  CLOSE(fp_out)
  ! geometry
  OPEN(fp_out, FILE=fn_geom, ACTION='write', STATUS='replace', IOSTAT=stata)
  CALL print_matrix(sph, CHANNEL=fp_out)
  CLOSE(fp_out)
  WRITE(*,*) "OK"
  
  
CONTAINS


!! Compute packing properties
SUBROUTINE comp_prop()
  IMPLICIT NONE
  INTEGER :: px, py, pz
  REAL(fp) :: d_min, d, cx, cy, cz, c_min, c, cr
  
  ! nearest neighbour
  DO i=1,n
    d_min = HUGE(d_min)
    c_min = HUGE(c_min)
    x = sph(i,1)
    y = sph(i,2)
    z = sph(i,3)
    cr = sph(i,4)
    DO j=1,n
      IF(j == i) CYCLE
      ! Periodic boxes
      DO  px=-1,1
        cx = x + px*box_x
        DO py=-1,1
          cy = y + py*box_y
          DO pz=-1,1
            cz = z + pz*box_z
            c = SQRT((sph(j,1)-cx)**2 + (sph(j,2)-cy)**2 + (sph(j,3)-cz)**2)
            IF(c < c_min) c_min = c
            d = c - (r+sph(j,4))
            IF(d < d_min) d_min = d
          END DO
        END DO
      END DO
    END DO
    near_neigh(i) = d_min
    near_cToc(i) = c_min
  END DO
  
END SUBROUTINE comp_prop


!! Pack spheres
SUBROUTINE do_pack()
  IMPLICIT NONE
  INTEGER :: px,py,pz
  REAL(fp) :: cx,cy,cz
  LOGICAL :: collide
  
  CALL rnd_xyz()
  sph(1,:) = (/ x, y, z, r /)
  
  WRITE(*,'(A24,I0,A3)', ADVANCE='no') "  packing particle # of ", n, " : "
  DO i=2,n
    IF(MOD(i,pack_report) == 0) THEN
      WRITE(*,'(I0,1X)', ADVANCE='no') i
    END IF
    collide = .TRUE.
    j = 1
    DO WHILE(collide)
      CALL rnd_xyz()
      collide = .FALSE.
      check_others: DO k=1,i-1
        ! Check periodic boxes
        DO px=-1,1
          cx = x + px*box_x
          DO py=-1,1
            cy = y + py*box_y
            DO pz=-1,1
              cz = z + pz*box_z
              IF((sph(k,1)-cx)**2 + (sph(k,2)-cy)**2 + (sph(k,3)-cz)**2 < &
                (sph(k,4)+r)**2) THEN
                COLLIDE = .TRUE.
                EXIT check_others
              END IF
            END DO
          END DO
        END DO
      END DO check_others
      j = j+1
      IF (j > pack_max_try) THEN
      WRITE(*,*) ""
        WRITE(*,'(A34,I0,A33)') "  ERROR: maximum number of tries (", pack_max_try, &
          ") reached in packing, stopping..."
        STOP
      END IF
    END DO
    sph(i,:) = (/ x, y, z, r /)
  END DO
  
END SUBROUTINE do_pack


!! alpha-version
SUBROUTINE handle_input()
  IMPLICIT NONE
  
  box_x = 10.0_fp
  box_y = 10.0_fp
  box_z = 10.0_fp
  
  WRITE(fn_geom, '(A9)'), "pack.geom"
  WRITE(fn_out, '(A8)'), "pack.out"
  
  sphVol = sphVolC * r**3
  boxVol = box_x*box_y*box_z

  ! give n
  ! n = 100
  ! matVol = n*sphVol
  ! pD = matVol/boxVol
  
  ! give pD
  pD = 0.1_fp
  n = NINT(pD*boxVol/sphVol)
  matVol = n*sphVol
  pD = matVol/boxVol
 
  ALLOCATE(sph(n,4),near_neigh(n),near_cToc(n), STAT=stata)
  IF(stata /= 0) THEN
    WRITE(*,*) "Error in memory allocation"
    STOP
  END IF

END SUBROUTINE handle_input


! random coordinates
SUBROUTINE rnd_xyz()
  IMPLICIT NONE
  
  x = give_rnr(0.0_fp,box_x)
  y = give_rnr(0.0_fp,box_y)
  z = give_rnr(0.0_fp,box_z)
  
END SUBROUTINE rnd_xyz


END PROGRAM PackingProperties
