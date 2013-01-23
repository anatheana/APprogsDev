!! 
!! Mesh volume
!! -
!! Module for computing the volume of closed, non-self-intersecting
!! triangular mesh. The algorithm is direct conversion of C algorithm
!! based on O' Rourke's Computational Geometry in C (2nd ed.).
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Mesh volume (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2013).
!! -
!! Antti Penttilä
!! 2013
!! Department of Physics, University of Helsinki

MODULE meshVolume

  USE Type_Kinds
  IMPLICIT NONE
  
  PUBLIC

CONTAINS

!! Calculate the volume of the tetrahedron given as input. The
!! calculation is carried out in the way described in Code 4.16 of
!! O'Rourke's Computational Geometry in C (2nd ed.), that is, by
!! translating the tetrahedron so that the vertex d is placed at the
!! origin.
FUNCTION tetrahedron_volume(tetra) RESULT(vol)
  IMPLICIT NONE
  REAL(FP_KIND), DIMENSION(4,3), INTENT(IN) :: tetra
  REAL(FP_KIND) :: vol
  
  REAL(FP_KIND) :: ax, ay, az, bx, by, bz, cx, cy, cz
  
  ax = tetra(1,1) - tetra(4,1)
  ay = tetra(1,2) - tetra(4,2)
  az = tetra(1,3) - tetra(4,3)
  bx = tetra(2,1) - tetra(4,1)
  by = tetra(2,2) - tetra(4,2)
  bz = tetra(2,3) - tetra(4,3)
  cx = tetra(3,1) - tetra(4,1)
  cy = tetra(3,2) - tetra(4,2)
  cz = tetra(3,3) - tetra(4,3)

  vol = (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + &
    az * (bx * cy - by * cx)) / 6.0_FP_KIND  

END FUNCTION tetrahedron_volume


!! An implementation of polyhedral volume calculation following the
!! divergence theorem and O' Rourke's Computational Geometry in C (2nd ed.).
!! The function arguments are as in exercise 4.7.7:
!! 
!! num_vertices is the number of vertices of the polyhedron
!! vectors is an array of 3D vectors of size num_vertices
!! num_trfaces is the number of triangle faces of the polyhedron.
!! vectors_indices is an array of vector indices, such that vector_indices(i,1),
!! vector_indices(i,2) and vector_indices(i,3) are the vectors composing
!! the triangular face i.
!!
!! It is assumed that the vertices are given in counter clockwise order.
FUNCTION calculate_volume(vectors, vector_indices) RESULT(volume)
  IMPLICIT NONE
  REAL(FP_KIND), DIMENSION(:,:), INTENT(IN) :: vectors
  INTEGER, DIMENSION(:,:), INTENT(IN) :: vector_indices
  REAL(FP_KIND) :: volume
  INTEGER :: i, num_trfaces
  REAL(FP_KIND), DIMENSION(4,3) :: tetra

  num_trfaces = SIZE(vector_indices, 1)
  volume = 0.0_FP_KIND
  tetra(4,:) = vectors(vector_indices(1,1),:)
  DO i=1,num_trfaces
    tetra(1,:) = vectors(vector_indices(i,1),:)
    tetra(2,:) = vectors(vector_indices(i,2),:)
    tetra(3,:) = vectors(vector_indices(i,3),:)
    volume = volume + tetrahedron_volume(tetra)
  END DO

END FUNCTION calculate_volume
  
END MODULE meshVolume
