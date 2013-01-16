!! Packing properties
!! -
!! Command line tool for simulating different properties and their distributions
!! in three-dimensional sphere packings, including nearest-neigbour.
!! Vessel including the spheres is a periodic minimum-bounding-box.
!! 
!! Spheres can have radii or they can be points (zero volume).
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Properties of Sphere packings (computer code),
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
  IMPLICIT NONE
  
  INTEGER :: i, j, k, sph_n
  
  REAL(FP_KIND) :: box_x, box_y, box_z, x, y, z, r, &
    pD, box_vol, sph_vol
  REAL(FP_KIND), DIMENSION(:,:), ALLOCATABLE :: sph
  
  LOGICAL :: point_process


  
  
CONTAINS

! First version
SUBROUTINE handle_input()
  IMPLICIT NONE
  
  point_process = .FALSE.
  pD = 0.1_FP_KIND
  box_x = 10.0_FP_KIND
  box_y = 10.0_FP_KIND
  box_z = 10.0_FP_KIND
  
  

END SUBROUTINE handle_input


END PROGRAM PackingProperties
