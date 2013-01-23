!! Mesh converter
!! -
!! Command line tool for converting different polygon mesh and voxel file formats.
!! Currently, possible input formats are:
!! .OBJ, open format polygonal mesh file
!! .MRT, Macke Ray Tracing file format
!! and possible output formats are
!! .POV, POV-Ray raytracing format
!! .ADD, Adda format for rectangular grid volume-based voxel format
!! .OBJ, open format polygonal mesh file
!! .NFO, general information about the structure defined in the input file
!! -
!! Usage: MeshConvert infile informat outfile outformat [reduce_vertices].
!! Note that input and output file formats are concluded from the command line argument
!! informat (obj/mrt) and outformat (pov/add/obj/nfo), and not from the file suffix of the input
!! or output file name. Optionally, keyword 'reduce_vertices' can be given. In that case the code
!! tries to look for dublicate vertex points and reduce the size of the polygonal mesh. This can be slow.
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Mesh converting software (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2012).
!! -
!! Antti Penttilä
!! 2012
!! Department of Physics, University of Helsinki
PROGRAM MeshConvert
  !! Double precision FP_KIND type parameter from PaulV
	USE Type_Kinds
  !! File handling routines from PaulV
	USE File_Utility
  !! String routines from PaulV
	USE String_Utility
  !! My own math- and printing-related routines
	USE AP_utils
  !! experimental
  USE meshVolume
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: str_len = 128, line_len = 512
	INTEGER, PARAMETER :: in_format_n = 2
	INTEGER, PARAMETER :: out_format_n = 4
	INTEGER, PARAMETER :: vertex_block = 1000
	INTEGER, PARAMETER :: face_block = 2000
	INTEGER, PARAMETER :: report = 1000
	INTEGER, PARAMETER :: iform_obj = 1, iform_mrt = 2, oform_pov = 1, oform_add = 2, oform_obj = 3, &
	  oform_nfo = 4
	
	CHARACTER(LEN=3), PARAMETER, DIMENSION(in_format_n) :: in_formats = (/ "obj", "mrt" /)
	CHARACTER(LEN=3), PARAMETER, DIMENSION(out_format_n) :: out_formats = (/ "pov", "add", "obj", "nfo" /)
	
	INTEGER :: i, j, k, argc, fp_in, fp_out, stats, form_in_i, form_out_i, &
		vertex_list_bn, vertex_list_n, face_list_bn, face_list_n
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: face_list, i3temp
	
	REAL(FP_KIND) :: x, y, z
	REAL(FP_KIND), DIMENSION(:,:), ALLOCATABLE :: vertex_list, f3temp, normal_list
	
	CHARACTER :: c
	CHARACTER(LEN=3) :: suf_in, suf_out
	CHARACTER(LEN=str_len) :: fn_in, fn_out, inp_str
	CHARACTER(LEN=line_len) :: line
	
	LOGICAL :: normals_on
	
	! External function, comment out with gfortran
	INTEGER, EXTERNAL :: IARGC
	
	! Command line input
	argc = IARGC()
	IF(argc /= 4 .AND. argc /= 5) THEN
		WRITE(*,*) "You must give four or more arguments, at least input file name and format, and ", &
			"output file name and format"
		WRITE(*,*) "  'shape.in obj shape.out pov'"
		STOP
	END IF
	CALL GETARG(1,fn_in)
	CALL GETARG(2,suf_in)
	CALL GETARG(3, fn_out)
	CALL GETARG(4,suf_out)
	
	CALL check_formats()
	
	CALL open_files()

	! Read input
	SELECT CASE(form_in_i)
	CASE(iform_obj) ! OBJ
		CALL read_obj()
	CASE(iform_mrt) ! MRT - Macke ray-tracing format
		CALL read_mrt()
	END SELECT
	CLOSE(fp_in)
	
	! Reduce vertices?
	IF(argc == 5) THEN
		CALL GETARG(5,inp_str)
		IF(StrLowCase(TRIM(inp_str)) == "reduce_vertices") THEN
			CALL reduce_vertices()
		END IF
	END IF

	! Write output
	SELECT CASE(form_out_i)
	CASE(oform_pov) ! POV
		CALL write_pov()
	CASE(oform_add) ! ADDA voxel format
		CALL write_add()		
	CASE(oform_obj) ! Wavefront OBJ format
		CALL write_obj()
	CASE(oform_nfo) ! General information
		CALL write_nfo()		
	END SELECT
	CLOSE(fp_out)	
	
	WRITE(*,*) ""
	WRITE(*,*) "Ready"
	
CONTAINS


! Tries to reduce the number of vertices by searching dublicates
SUBROUTINE reduce_vertices()
	IMPLICIT NONE
	
	INTEGER :: li, lj, lk, rvc, i_is, i_was
	INTEGER, DIMENSION(:), ALLOCATABLE :: vertex_reduce_list, union_list
	REAL(FP_KIND), DIMENSION(:,:), ALLOCATABLE :: vertex_list_temp, normal_list_temp
	
  WRITE(*,*) ""
  WRITE(*,*) "Starting vertex reduction..."

	ALLOCATE(vertex_reduce_list(vertex_list_n),union_list(vertex_list_n), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
	
	DO li=1,vertex_list_n
		vertex_reduce_list(li) = li
	END DO
	
  WRITE(*,'(A,I0,A)', ADVANCE='no') "Searching for duplicate vertices, handling # of ", vertex_list_n, ":"
	rvc = 1
	mainloop: DO li=2,vertex_list_n
	  IF(MOD(li,report) == 0) THEN
  	  WRITE(*,'(A1,I0)', ADVANCE='no') " ", li
  	END IF
		searchloop: DO lj=1,li-1
			IF(vertex_list(lj,1) == vertex_list(li,1) .AND. &
				vertex_list(lj,2) == vertex_list(li,2) .AND. &
				vertex_list(lj,3) == vertex_list(li,3)) THEN
				vertex_reduce_list(li) = lj
				CYCLE mainloop
			END IF
		END DO searchloop
		rvc = rvc+1
	END DO mainloop
	
	IF(normals_on) THEN
		ALLOCATE(vertex_list_temp(rvc,3),normal_list_temp(rvc,3), STAT=stats)
	ELSE
		ALLOCATE(vertex_list_temp(rvc,3), STAT=stats)
	END IF
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
	
  WRITE(*,*) ""
  WRITE(*,*) "Computing sorted union of vertex list"
	li = 0
	union_list(:) = vertex_reduce_list(:)
	CALL union_integer(union_list, li)
	
  WRITE(*,'(A,I0,A)', ADVANCE='no') "Updating vertex list, handling # of ", rvc, " :"
	DO li=1,rvc
	  IF(MOD(li,report) == 0) THEN
  	  WRITE(*,'(A1,I0)', ADVANCE='no') " ", li
  	END IF
		i_was = union_list(li)
		i_is = li
		vertex_list_temp(i_is,:) = vertex_list(i_was,:)
		IF(normals_on) normal_list_temp(i_is,:) = normal_list(i_was,:)
		DO lj=1,vertex_list_n
			IF(vertex_reduce_list(lj) == i_was) THEN
				vertex_reduce_list(lj) = i_is
			END IF
		END DO
	END DO
	DEALLOCATE(vertex_list)
	IF(normals_on) THEN
		DEALLOCATE(normal_list)
		ALLOCATE(vertex_list(rvc,3),normal_list(rvc,3), STAT=stats)
	ELSE
		ALLOCATE(vertex_list(rvc,3), STAT=stats)
	END IF
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
	vertex_list(:,:) = vertex_list_temp(1:rvc,:)
	DEALLOCATE(vertex_list_temp)
	IF(normals_on) THEN
		normal_list(:,:) = normal_list_temp(1:rvc,:)
		DEALLOCATE(normal_list_temp)
	END IF

	WRITE(*,*) ""
	WRITE(*,'(A,I0,A)', ADVANCE='no') "Updating face list, handling # of ", face_list_n, " :"
	DO li=1,face_list_n
	  IF(MOD(li,report) == 0) THEN
  	  WRITE(*,'(A1,I0)', ADVANCE='no') " ", li
  	END IF
		DO lj=1,3
			i_was = face_list(li,lj)
			i_is = vertex_reduce_list(i_was)
			face_list(li,lj) = i_is
		END DO
	END DO
	DEALLOCATE(vertex_reduce_list)
	
	WRITE(*,*) ""
	WRITE(*,*) "Vertex reduction OK."
	WRITE(*,*) "Compacting from ", vertex_list_n, " to ", rvc, " vertices"
	
	
	vertex_list_n = rvc
	vertex_list_bn = vertex_block*(vertex_list_n/vertex_block+1)

END SUBROUTINE reduce_vertices


! Writes general information about the target
SUBROUTINE write_nfo()
	IMPLICIT NONE

	INTEGER :: li, lj, fa, fb, fc
	REAL(FP_KIND) :: min_x, max_x, min_y, max_y, min_z, max_z, triA, max_triA, min_triA, triA_temp, vol
	REAL(FP_KIND), DIMENSION(3) :: av, bv, cv
	
	min_x = MINVAL(vertex_list(1:vertex_list_n,1))
	max_x = MAXVAL(vertex_list(1:vertex_list_n,1))
	min_y = MINVAL(vertex_list(1:vertex_list_n,2))
	max_y = MAXVAL(vertex_list(1:vertex_list_n,2))
	min_z = MINVAL(vertex_list(1:vertex_list_n,3))
	max_z = MAXVAL(vertex_list(1:vertex_list_n,3))
	
	min_triA = HUGE(min_triA)
	max_triA = 0.0_FP_KIND
  triA = 0.0_FP_KIND
  DO li=1,face_list_n
    fa = face_list(li,1)
    fb = face_list(li,2)
    fc = face_list(li,3)
    av = vertex_list(fa,:)
    bv = vertex_list(fb,:)
    cv = vertex_list(fc,:)
    triA_temp = sqrt((av(2)*(bv(1) - cv(1)) + bv(2)*cv(1) - bv(1)*cv(2) + av(1)*(-bv(2) + cv(2)))**2 + &
      (av(3)*(bv(1) - cv(1)) + bv(3)*cv(1) - bv(1)*cv(3) + av(1)*(-bv(3) + cv(3)))**2 + &
      (av(3)*(bv(2) - cv(2)) + bv(3)*cv(2) - bv(2)*cv(3) + av(2)*(-bv(3) + cv(3)))**2)/2.0_FP_KIND
    triA = triA + triA_temp
    IF(triA_temp < min_triA) min_triA = triA_temp
    IF(triA_temp > max_triA) max_triA = triA_temp
   END DO

	
	WRITE(fp_out, '(A,A)') "Object information, read from file ", TRIM(fn_in)
	WRITE(fp_out,*) "Contains ", face_list_n, " triangles using ", vertex_list_n, " vertices."
	WRITE(fp_out,*) "Object range is ", min_x, " to ", max_x, " (in x)"
	WRITE(fp_out,*) "  ", min_y, " to ", max_y, " (in y)"
	WRITE(fp_out,*) "  ", min_z, " to ", max_z, " (in z)"
	WRITE(fp_out,*) "Total face area is ", triA, ", and"
	WRITE(fp_out,*) "  min. area ", min_triA
	WRITE(fp_out,*) "  mean area ", triA/face_list_n
	WRITE(fp_out,*) "  max. area ", max_triA

	WRITE(*,*) ""
	WRITE(*,*) "NFO format write OK."
  
  vol = calculate_volume(vertex_list(1:vertex_list_n,:), face_list(1:face_list_n,:))
  WRITE(*,*) "Experimental: volume from mesh is ", vol
  WRITE(fp_out,*) "Experimental: volume from mesh is ", vol
	
END SUBROUTINE write_nfo


! Write ADDA voxel geometry format
SUBROUTINE write_add()
	IMPLICIT NONE

	INTEGER :: li, lj, lk, xn, yn, zn, cxn, cc
	REAL(FP_KIND) :: min_x, max_x, min_y, max_y, min_z, max_z, b, x, y, z, hb, c1, c2, c3, t1, rv, &
		jy, jz, vol1, vol2
	REAL(FP_KIND), DIMENSION(3) :: u1, u2, u3
	LOGICAL :: does_cross
	LOGICAL, DIMENSION(:), ALLOCATABLE :: scanline
	
	IF(argc == 4) THEN
		b = 1.0_FP_KIND
	ELSE
		CALL GETARG(5, inp_str)
		READ(inp_str, *) b
	END IF
	hb = 0.5_FP_KIND * b
	
	CALL init_rnd()
	
	min_x = MINVAL(vertex_list(1:vertex_list_n,1))
	max_x = MAXVAL(vertex_list(1:vertex_list_n,1))
	min_y = MINVAL(vertex_list(1:vertex_list_n,2))
	max_y = MAXVAL(vertex_list(1:vertex_list_n,2))
	min_z = MINVAL(vertex_list(1:vertex_list_n,3))
	max_z = MAXVAL(vertex_list(1:vertex_list_n,3))
	
	xn = CEILING((max_x-min_x)/b)+2
	yn = CEILING((max_y-min_y)/b)+2
	zn = CEILING((max_z-min_z)/b)+2

	WRITE(fp_out, '(A,A)') "# ADDA geometry file from ", TRIM(fn_in)
	WRITE(fp_out, '(A,ES14.8)') "# cellsize ", b
	WRITE(fp_out, '(A,I0,A,I0,A,I0)') "# cells: x-dir ", xn, ", y-dir ", yn, ", z-dir ", zn

	ALLOCATE(scanline(xn), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF

	cc = 0
	x = min_x - hb
	WRITE(*,*) ""
	WRITE(*,'(A9,I0,A17)',ADVANCE='no') "Scanning ", yn*zn, " lines, now ready"
	yloop: DO lj=1,yn
	  WRITE(*,'(A1,F6.2,A1)',ADVANCE='no') " ", 100.0_FP_KIND*(1+(lj-1)*zn)/(yn*zn), "%"
		CALL RANDOM_NUMBER(rv)
		y = min_y - hb + b*(lj-1)
		jy = y + (rv-0.5_FP_KIND)*0.01_FP_KIND*b
		zloop: DO lk=1,zn
			CALL RANDOM_NUMBER(rv)
			z = min_z - hb + b*(lk-1)
			jz = z + (rv-0.5_FP_KIND)*0.01_FP_KIND*b
			scanline(:) = .FALSE.			
			scanloop: DO li=1,face_list_n
				u1(:) = vertex_list(face_list(li,1),:)
				u2(:) = vertex_list(face_list(li,2),:)
				u3(:) = vertex_list(face_list(li,3),:)
				does_cross = poly_crossing(u1,u2,u3,x,jy,jz,c1,c2,c3)
				IF(does_cross) THEN
					t1 = -(c3-hb)/b
					cxn = INT(t1)
					IF(cxn /= t1) cxn=cxn+1
					scanline(cxn+1:xn) = (scanline(cxn+1:xn) .EQV. .FALSE.)
				END IF
			END DO scanloop
			DO li=1,xn
				IF(scanline(li)) THEN
					WRITE(fp_out,'(I0,1X,I0,1X,I0)') li, lj, lk
					cc = cc+1
				END IF
			END DO
		END DO zloop
	END DO yloop
	WRITE(*,*) ""

	WRITE(*,*) ""
	WRITE(*,*) "ADDA format write OK."
	WRITE(*,*) "Wrote ", cc, " cells in ", xn, "x", yn, "x", zn, " grid with cellsize ", b
	vol1 = cc*b**3
	vol2 = (max_x-min_x)*(max_y-min_y)*(max_z-min_z)
	WRITE(*,'(A,ES12.5,A,F7.3,A)') " Cell volume is ", cc*b**3, " and porosity inside surrounding box is ", vol1/vol2, "%"
	WRITE(*,*) "Object range was ", min_x, " to ", max_x, " (in x)"
	WRITE(*,*) "  ", min_y, " to ", max_y, " (in y)"
	WRITE(*,*) "  ", min_z, " to ", max_z, " (in z)"

END SUBROUTINE write_add


! Write POV-Ray file format
SUBROUTINE write_pov()
	IMPLICIT NONE

	INTEGER :: li, lj
	REAL(FP_KIND) :: mr
	REAL(FP_KIND), DIMENSION(3) :: big, small, center
	
	big = MAXVAL(vertex_list)
	small = MINVAL(vertex_list)
	center = SUM(vertex_list)/vertex_list_n
	mr = MAX(big(1)-center(1),center(1)-small(1),big(2)-center(2),center(2)-small(2), &
		big(3)-center(3),center(3)-small(3))
	
	WRITE(fp_out, *) "#declare box_max = <", big(1), ",", big(2), ",", big(3), ">;"
	WRITE(fp_out, *) "#declare box_center = <", center(1), ",", center(2), ",", center(3), ">;"
	WRITE(fp_out, *) "#declare box_min = <", small(1), ",", small(2), ",", small(3), ">;"
	WRITE(fp_out, *) "#declare max_r = ", mr, ";"
	
	WRITE(fp_out, *) "mesh2 {"
	WRITE(fp_out, *) " vertex_vectors {"
	WRITE(fp_out, *) "  ", vertex_list_n, ","
	DO li=1,vertex_list_n-1
		WRITE(fp_out, *) "  <", vertex_list(li,1), ",", vertex_list(li,2), ",", vertex_list(li,3), ">,"
	END DO
	WRITE(fp_out, *) "  <", vertex_list(li,1), ",", vertex_list(li,2), ",", vertex_list(li,3), ">"
	WRITE(fp_out, *) " }"
	IF(normals_on) THEN
		WRITE(fp_out, *) " normal_vectors {"
		WRITE(fp_out, *) "  ", vertex_list_n, ","
		DO li=1,vertex_list_n-1
			WRITE(fp_out, *) "  <", normal_list(li,1), ",", normal_list(li,2), ",", normal_list(li,3), ">,"
		END DO
		WRITE(fp_out, *) "  <", normal_list(li,1), ",", normal_list(li,2), ",", normal_list(li,3), ">"
		WRITE(fp_out, *) " }"
	END IF
	WRITE(fp_out, *) " face_indices {"
	WRITE(fp_out, *) "  ", face_list_n, ","
	DO li=1,face_list_n-1
		WRITE(fp_out, *) "  <", face_list(li,1)-1, ",", face_list(li,2)-1, ",", &
			face_list(li,3)-1, ">,"
	END DO
	WRITE(fp_out, *) "  <", face_list(li,1)-1, ",", face_list(li,2)-1, ",", &
		face_list(li,3)-1, ">"
	WRITE(fp_out, *) " }"
	WRITE(fp_out, *) " material1()"
	WRITE(fp_out, *) "}"

	WRITE(*,*) ""
	WRITE(*,*) "POV-Ray format write OK."
	
END SUBROUTINE write_pov


! Write Wavefront OBJ file format
SUBROUTINE write_obj()
	IMPLICIT NONE

	INTEGER :: li, lj
	
	WRITE(fp_out, '(A,A)') "# Converted from file ", TRIM(fn_in)
	WRITE(fp_out, '(A)') "g default"
	WRITE(fp_out, '(A)') "s off"
	
	DO li=1,vertex_list_n
		WRITE(fp_out, '(A2,F12.8,A1,F12.8,A1,F12.8)') "v ", vertex_list(li,1), " ", vertex_list(li,2), " ", &
			vertex_list(li,3)
	END DO
	IF(normals_on) THEN
		DO li=1,vertex_list_n
			WRITE(fp_out, '(A3,F12.8,A1,F12.8,A1,F12.8)') "vn ", normal_list(li,1), " ", normal_list(li,2), " ", &
				normal_list(li,3)
		END DO
	END IF
	WRITE(fp_out, '(A)') "o object0"
	IF(normals_on) THEN
		DO li=1,face_list_n
			WRITE(fp_out, '(A2,I0,A2,I0,A1,I0,A2,I0,A1,I0,A2,I0)') "f ", face_list(li,1), "//", &
				face_list(li,1), " ", face_list(li,2), "//", face_list(li,2), " ", face_list(li,3), &
				"//", face_list(li,3)
		END DO
	ELSE
		DO li=1,face_list_n
			WRITE(fp_out, '(A2,I0,A1,I0,A1,I0)') "f ", face_list(li,1), " ", face_list(li,2), " ", &
				face_list(li,3)
		END DO
	END IF
	
	WRITE(*,*) ""
	WRITE(*,*) "OBJ format write OK."
	
END SUBROUTINE write_obj


! Read Macke's ray-tracing file format
SUBROUTINE read_mrt()
  IMPLICIT NONE
  INTEGER :: li, lj, lk, pn, fn
  INTEGER, DIMENSION(:), ALLOCATABLE :: vn, two
  
  READ(fp_in, *) pn
  ALLOCATE(vn(pn),two(pn), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
  DO li=1,pn
    READ(fp_in, *, IOSTAT=stats) vn(li)
    IF(stats /= 0) THEN
      WRITE(*,*) "Error (", stats, ") in reading 'mrt' format file ", TRIM(fn_in)
      WRITE(*,*) "File ends or other I/O error while reading ", pn, " polygons at line ", li+1
      STOP
    END IF
    SELECT CASE(vn(li))
    CASE(3)
    	two(li) = 1
    CASE(4)
    	two(li) = 2
    CASE DEFAULT
			WRITE(*,*) "Error: this program supports only polygons with 3 or 4 vertices."
			STOP
		END SELECT
  END DO
  face_list_n = SUM(two)
  vertex_list_n = SUM(vn)
  
  face_list_bn = face_block*(face_list_n/face_block+1)
  vertex_list_bn = vertex_block*(vertex_list_n/vertex_block+1)
  
  ALLOCATE(vertex_list(vertex_list_bn,3), face_list(face_list_bn,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
	
	DO li=1,vertex_list_n
		READ(fp_in, *, IOSTAT=stats) vertex_list(li,:)
    IF(stats /= 0) THEN
      WRITE(*,*) "Error (", stats, ") in reading 'mrt' format file ", TRIM(fn_in)
      WRITE(*,*) "File ends or other I/O error while reading vertex ", li, &
        " of total ", vertex_list_n, " vertices at line ", li+1+pn
      STOP
    END IF
	END DO

	lj = 0
	lk = 0
	DO li=1,pn
		lj = lj+1
		face_list(lj,1) = lk+1
		face_list(lj,2) = lk+2
		face_list(lj,3) = lk+3
		IF(two(li) == 1) THEN
			lk = lk+3
		ELSE
			lj = lj+1
			face_list(lj,1) = lk+1
			face_list(lj,2) = lk+3
			face_list(lj,3) = lk+4
			lk = lk+4
		END IF
	END DO

  DEALLOCATE(vn)

	WRITE(*,*) ""
	WRITE(*,*) "Macke ray-tracing format read OK"
	WRITE(*,*) "Read ", face_list_n, " triangles using ", vertex_list_n, " vertices."
	
END SUBROUTINE read_mrt


! Read Wavefront OBJ file format
SUBROUTINE read_obj()
	IMPLICIT NONE
	INTEGER :: li, lj, lk, fi1, fi2, away
	CHARACTER(LEN=2) :: c2
	CHARACTER :: c1

	ALLOCATE(vertex_list(vertex_block,3), face_list(face_block,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation failed"
		STOP
	END IF
	vertex_list_bn = vertex_block
	face_list_bn = face_block
	
	li = 0
	lj = 0
	lk = 0
	DO
		IF(li == vertex_list_bn) CALL allocate_vertex()
		IF(lj == face_list_bn) CALL allocate_face()
		READ(fp_in, '(A)', IOSTAT=stats) line
		IF(stats /= 0) EXIT

		c = line(1:1)
		SELECT CASE(c)
		CASE('v')
			IF(line(2:2) == 'n') THEN
				IF(.NOT. normals_on) THEN ! First time normals
					normals_on = .TRUE.
					ALLOCATE(normal_list(vertex_list_bn,3), STAT=stats)
					IF(stats /= 0) THEN
						WRITE(*,*) "Memory allocation failed"
						STOP
					END IF
				END IF
				lk = lk+1
				IF(lk > vertex_list_bn) THEN
					WRITE(*,*) "Number of normals exceeding number of vertices"
					STOP
				END IF
				READ(line(3:),*) normal_list(lk,:)
			ELSE
				li = li+1
				READ(line(2:),*) vertex_list(li,:)
			END IF
		CASE('f')
			lj = lj+1
			fi2 = INDEX(line, "//")
			fi1 = INDEX(line, "/")
			IF(fi2 > 0) THEN
				CALL mytrim(line, '/')
				READ(line(2:),*) face_list(lj,1), away, face_list(lj,2), away, face_list(lj,3), away
			ELSE IF(fi1 > 0) THEN
				CALL mytrim(line, '/')
				READ(line(2:),*) face_list(lj,1), away,away, face_list(lj,2), away,away, &
					face_list(lj,3), away,away
			ELSE
				READ(line(2:),*) face_list(lj,:)
			END IF
		END SELECT
	END DO
	
	vertex_list_n = li
	face_list_n = lj
	IF(normals_on .AND. li /= lk) THEN
		WRITE(*,*) "Number of normals does not equal number of vertices, which is not supported here."
		STOP
	END IF

	WRITE(*,*) ""
	WRITE(*,*) "OBJ format read OK."
	WRITE(*,*) "Read ", face_list_n, " triangles using ", vertex_list_n, " vertices."

END SUBROUTINE read_obj


SUBROUTINE allocate_vertex()
	IMPLICIT NONE
	
	INTEGER :: bn
	
	bn = vertex_list_bn+vertex_block
	ALLOCATE(f3temp(vertex_list_bn,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation error"
		STOP
	END IF
	f3temp(:,:) = vertex_list(:,:)
	DEALLOCATE(vertex_list, STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory deallocation error"
		STOP
	END IF
	ALLOCATE(vertex_list(bn,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation error"
		STOP
	END IF
	vertex_list(1:vertex_list_bn,:) = f3temp(:,:)
	
	! Normals, if needed
	IF(normals_on) THEN
		f3temp(:,:) = normal_list(:,:)
		DEALLOCATE(normal_list, STAT=stats)
		IF(stats /= 0) THEN
			WRITE(*,*) "Memory deallocation error"
			STOP
		END IF
		ALLOCATE(normal_list(bn,3), STAT=stats)
		IF(stats /= 0) THEN
			WRITE(*,*) "Memory allocation error"
			STOP
		END IF
		normal_list(1:vertex_list_bn,:) = f3temp(:,:)
	END IF
	
	DEALLOCATE(f3temp, STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory deallocation error"
		STOP
	END IF
	vertex_list_bn = bn

END SUBROUTINE allocate_vertex


SUBROUTINE allocate_face()
	IMPLICIT NONE
	
	INTEGER :: bn
	
	bn = face_list_bn+face_block
	ALLOCATE(i3temp(face_list_bn,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation error"
		STOP
	END IF
	i3temp(:,:) = face_list(:,:)
	DEALLOCATE(face_list, STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory deallocation error"
		STOP
	END IF
	ALLOCATE(face_list(bn,3), STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory allocation error"
		STOP
	END IF
	face_list(1:face_list_bn,:) = i3temp(:,:)
	DEALLOCATE(i3temp, STAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Memory deallocation error"
		STOP
	END IF
	face_list_bn = bn
	
END SUBROUTINE allocate_face

SUBROUTINE open_files()
	IMPLICIT NONE

	! Open files and channels
	fp_in = Get_Lun()
	OPEN(fp_in, FILE=fn_in, STATUS='old', ACTION='read', IOSTAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Error in opening file '", TRIM(fn_in), "' for reading"
		STOP
	END IF
	fp_out = Get_Lun()
	OPEN(fp_out, FILE=fn_out, STATUS='replace', ACTION='write', IOSTAT=stats)
	IF(stats /= 0) THEN
		WRITE(*,*) "Error in opening file '", TRIM(fn_out), "' for writing"
		STOP
	END IF
	
END SUBROUTINE open_files


SUBROUTINE check_formats()
	IMPLICIT NONE
	INTEGER :: li

	WRITE(*,*) ""

	! Check format codes
	form_in_i = 0
	DO li=1,in_format_n
		IF(StrLowCase(in_formats(li)) == StrLowCase(suf_in)) THEN
			form_in_i = li
			EXIT
		END IF
	END DO
	IF(form_in_i == 0) THEN
		WRITE(*,*) "Unknown input format '", suf_in, "'"
		STOP
	ELSE
		WRITE(*,*) "Input file format: '", in_formats(form_in_i), "'"
	END IF
	form_out_i = 0
	DO li=1,out_format_n
		IF(StrLowCase(out_formats(li)) == StrLowCase(suf_out)) THEN
			form_out_i = li
			EXIT
		END IF
	END DO
	IF(form_out_i == 0) THEN
		WRITE(*,*) "Unknown output format '", suf_out, "'"
		STOP
	ELSE
		WRITE(*,*) "Output file format: '", out_formats(form_out_i), "'"
	END IF
	
END SUBROUTINE check_formats


SUBROUTINE mytrim(str, ch, REP)
	IMPLICIT NONE
	CHARACTER(LEN=*), INTENT(INOUT) :: str
	CHARACTER, INTENT(IN) :: ch
	CHARACTER, INTENT(IN), OPTIONAL :: REP
	INTEGER :: li, sl
	CHARACTER :: lrep
	
	IF(PRESENT(REP)) THEN
		lrep = REP
	ELSE
		lrep = ' '
	END IF
	
	sl = LEN_TRIM(str)
	DO li=1,sl
		IF(str(li:li) == ch) WRITE(str(li:li), '(A1)') lrep
	END DO

END SUBROUTINE mytrim


FUNCTION poly_crossing(u1,u2,u3,x,y,z,c1,c2,c3) RESULT(does_cross)
	IMPLICIT NONE

	REAL(FP_KIND), DIMENSION(3), INTENT(IN) :: u1,u2,u3
	REAL(FP_KIND), INTENT(IN) :: x,y,z
	REAL(FP_KIND), INTENT(OUT) :: c1,c2,c3
	LOGICAL :: does_cross
	REAL(FP_KIND) :: div

	c1 = -1.0_FP_KIND
	c2 = -1.0_FP_KIND
	c3 = 0.0_FP_KIND
	does_cross = .FALSE.
	
	div = (u1(3)*(u2(2) - u3(2)) + u2(3)*u3(2) - u2(2)*u3(3) + u1(2)*(-u2(3) + u3(3)))
	IF(div == 0) RETURN
	c1 = (-(z*u1(2)) + y*u1(3) + z*u3(2) - u1(3)*u3(2) + (-y + u1(2))*u3(3)) / div
	IF(c1 < 0 .OR. 1 < c1) RETURN

	div = (-((z - u2(3))*(u1(2) - u3(2))) + (y - u2(2))*(u1(3) - u3(3)))
	IF(div == 0) RETURN
	c2 = (-(z*u1(2)) + y*u1(3) + z*u2(2) - u1(3)*u2(2) + (-y + u1(2))*u2(3)) / div
	IF(c2 < 0 .OR. 1 < c2) RETURN

	div = (u1(3)*(u2(2) - u3(2)) + u2(3)*u3(2) - u2(2)*u3(3) + u1(2)*(-u2(3) + u3(3)))
	IF(div == 0) RETURN
	c3 = (-(y*u1(3)*u2(1)) + x*u1(3)*u2(2) + y*u1(1)*u2(3) - x*u1(2)*u2(3) + y*u1(3)*u3(1) - &
		u1(3)*u2(2)*u3(1) - y*u2(3)*u3(1) + u1(2)*u2(3)*u3(1) - x*u1(3)*u3(2) + &
		u1(3)*u2(1)*u3(2) + x*u2(3)*u3(2) - u1(1)*u2(3)*u3(2) + &
		z*(u1(2)*(u2(1) - u3(1)) + u2(2)*u3(1) - u2(1)*u3(2) + u1(1)*(-u2(2) + u3(2))) + & 
		(-(y*u1(1)) + x*u1(2) + y*u2(1) - u1(2)*u2(1) + (-x + u1(1))*u2(2))*u3(3)) / div
	does_cross = .TRUE.

END FUNCTION poly_crossing 


SUBROUTINE init_rnd
	IMPLICIT NONE

	INTEGER :: rsize
	INTEGER, DIMENSION(8) :: t
	INTEGER, DIMENSION(:), ALLOCATABLE :: rseed

	CALL RANDOM_SEED(size=rsize)
	ALLOCATE(rseed(rsize))

	CALL DATE_AND_TIME(VALUES=t)
	rseed = 100*t(7) + t(8)/10
	CALL RANDOM_SEED(PUT=rseed)

END SUBROUTINE init_rnd


END PROGRAM MeshConvert
