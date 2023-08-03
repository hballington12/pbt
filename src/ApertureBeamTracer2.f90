! ApertureBeamTracer1.f90
! main program for beam tracing code developed by Harry Ballington

program main

use input_mod
use misc_submod
use beam_loop_mod
use types_mod
use diff_mod
use omp_lib
use outputs_mod

implicit none

! to do:
! add vector re-normalisation checks after rotations (ie. sr rotate_into_aperture_system)

! shared
real(8) start, finish ! cpu timing variables

! input
character(len=*), parameter :: ifn = 'input.txt' ! input filename

! sr SDATIN
character(100) cfn ! crystal filename
character(100) afn ! apertures filename
real(8) la ! wavelength
real(8) rbi ! real part of the refractive index
real(8) ibi ! imaginary part of the refractive index
integer(8) rec ! max number of internal beam recursions
character(100) rot_method ! rotation method
logical is_multithreaded ! whether or not code should use multithreading

! sr PDAL2
integer(8) num_vert ! number of unique vertices
integer(8) num_face !  number of faces
integer(8) num_norm ! number of face normals
integer(8), dimension(:,:), allocatable :: face_ids ! face vertex IDs
real(8), dimension(:,:), allocatable :: vert ! unique vertices
real(8), dimension(:,:), allocatable :: norm ! face normals
integer(8), dimension(:), allocatable :: num_face_vert ! number of vertices in each face
integer(8), dimension(:), allocatable :: norm_ids ! face normal ID of each face

! sr makeIncidentBeam
real(8), allocatable, dimension(:,:) :: beamV ! beam vertices
real(8), allocatable, dimension(:,:) :: beamN ! beam normals
real(8), allocatable, dimension(:,:) :: beamMidpoints ! beam  midpoints
integer(8), allocatable, dimension(:,:) :: beamF1 ! beam face vertex indices
integer(8), allocatable, dimension(:) :: beamF2 ! beam face normal indices
complex(8), allocatable, dimension(:,:,:) :: ampl_beam ! amplitude matrix of incident beam

! sr beam_loop
type(outbeamtype), dimension(:), allocatable :: beam_outbeam_tree ! outgoing beams from the beam tracing
type(outbeamtype), dimension(:), allocatable :: ext_diff_outbeam_tree ! outgoing beams from external diffraction
integer(8) beam_outbeam_tree_counter ! counts the current number of beam outbeams

! sr diff_main
complex(8), dimension(:,:), allocatable:: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! total
real(8), dimension(:), allocatable :: theta_vals, phi_vals
complex(8), dimension(:,:), allocatable :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! total

! ############################################################################################################

print*,'========== start main'
start = omp_get_wtime()

! ############# input_mod #############

! read input parameters
call SDATIN(ifn,            & ! <-  input filename
            cfn,            & !  -> crystal filename
            la,             & !  -> wavelength
            rbi,            & !  -> real part of the refractive index
            ibi,            & !  -> imaginary part of the refractive index
            afn,            & !  -> apertures filename
            rec,            & !  -> max number of internal beam recursions
            rot_method,     & !  -> particle rotation method
            is_multithreaded) !  whether ot not code should use multithreading

! read particle file
call PDAL2( cfn,            & ! <-  crystal filename
            num_vert,       & !  -> number of unique vertices
            num_face,       & !  -> number of faces
            num_norm,       & !  -> number of face normals (to do: remove)
            face_ids,       & !  -> face vertex IDs
            vert,           & !  -> unique vertices
            norm,           & !  -> face normals
            num_face_vert,  & !  -> number of vertices in each face
            norm_ids)         !  -> face normal ID of each face (to do: remove)

! rotate particle
call PROT(  ifn,        & ! <-  input filename
            rot_method, & ! <-  particle rotation method
            vert)         ! <-> unique vertices (unrotated in, rotated out) to do: remove inout intent and add a rotated vertices variable
! stop
! write rotated particle to file (optional)            
call PDAS(  vert,       & ! <-  rotated vertices
            face_ids)        ! <-  face vertex IDs
! stop
! fast implementation of the incident beam
call makeIncidentBeam(  beamV,         & ! ->  beam vertices
                        beamF1,        & ! ->  beam face vertex indices
                        beamN,         & ! ->  beam normals
                        beamF2,        & ! ->  beam face normal indices
                        vert,          & ! <-  unique vertices
                        face_ids,      & ! <- face vertex IDs
                        beamMidpoints, & !  -> beam  midpoints
                        ampl_beam)       !  -> amplitude matrix of incident beam       

! beam loop
call beam_loop(   face_ids,                  & ! <-  face vertex IDs
                  vert,                      & ! <-  unique vertices
                  la,                        & ! <-  wavelength
                  rbi,                       & ! <-  real part of the refractive index
                  ibi,                       & ! <-  imaginary part of the refractive index
                  afn,                       & ! <-  apertures filename
                  rec,                       & ! <-  max number of internal beam recursions
                  beamV,                     & ! <-  beam vertices
                  beamF1,                    & ! <-  beam face vertex indices
                  beamN,                     & ! <-  beam normals
                  beamF2,                    & ! <-  beam face normal indices
                  beamMidpoints,             & ! <-  beam  midpoints
                  ampl_beam,                 & ! <-  amplitude matrix of incident beam
                  beam_outbeam_tree,         & !  -> outgoing beams from the beam tracing
                  beam_outbeam_tree_counter, & !  -> counts the current number of beam outbeams
                  ext_diff_outbeam_tree)       !  -> outgoing beams from external diffraction

! diffraction
call diff_main(   beam_outbeam_tree,         & ! <-  outgoing beams from the beam tracing
                  beam_outbeam_tree_counter, & ! <-  counts the current number of beam outbeams
                  la,                        & ! <-  wavelength
                  ampl_far_beam11,           & !  -> amplitude matrix (1,1) due to beam diffraction
                  ampl_far_beam12,           & !  -> amplitude matrix (1,2) due to beam diffraction
                  ampl_far_beam21,           & !  -> amplitude matrix (2,1) due to beam diffraction
                  ampl_far_beam22,           & !  -> amplitude matrix (2,2) due to beam diffraction
                  theta_vals,                & !  -> theta values
                  phi_vals,                  & !  -> phi values
                  ext_diff_outbeam_tree,     & ! <-  outgoing beams from external diffraction
                  ampl_far_ext_diff11,       & !  -> amplitude matrix (1,1) due to external diffraction
                  ampl_far_ext_diff12,       & !  -> amplitude matrix (1,2) due to external diffraction
                  ampl_far_ext_diff21,       & !  -> amplitude matrix (2,1) due to external diffraction
                  ampl_far_ext_diff22,       & !  -> amplitude matrix (2,2) due to external diffraction
                  is_multithreaded)            ! enable or disable multithreading


call outputs(  ampl_far_beam11,     & ! <-  amplitude matrix (1,1) due to beam diffraction
               ampl_far_beam12,     & ! <-  amplitude matrix (1,2) due to beam diffraction
               ampl_far_beam21,     & ! <-  amplitude matrix (2,1) due to beam diffraction
               ampl_far_beam22,     & ! <-  amplitude matrix (2,2) due to beam diffraction
               ampl_far_ext_diff11, & ! <-  amplitude matrix (1,1) due to external diffraction
               ampl_far_ext_diff12, & ! <-  amplitude matrix (1,2) due to external diffraction
               ampl_far_ext_diff21, & ! <-  amplitude matrix (2,1) due to external diffraction
               ampl_far_ext_diff22, & ! <-  amplitude matrix (2,2) due to external diffraction
               theta_vals,          & ! <-  theta values
               phi_vals)              ! <-  phi values



! call test_ref()

finish = omp_get_wtime()
print'(A,f16.8,A)',"total time elapsed: ",finish-start," secs"
print*,'========== end main'

contains



end program main

