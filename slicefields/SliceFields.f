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
  !! experimental
  USE meshVolume
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: str_len = 128, line_len = 512

	
	WRITE(*,*) ""
	WRITE(*,*) "Ready"
	
CONTAINS


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



END PROGRAM SliceFields
