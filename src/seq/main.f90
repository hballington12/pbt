! ############################################################################################################

! abt.f90
! main program for physical optics hybrid method beam tracing code developed by Harry Ballington and Evelyn Hesse
! see README file for more details

program main

use input_mod
use misc_submod
use beam_loop_mod
use types_mod
use diff_mod
use omp_lib
use outputs_mod
use ifport

implicit none

! to do:
! add vector re-normalisation checks after rotations (ie. sr rotate_into_aperture_system)
! add quad support
! add support to avoid crash if nan detected
! automatic meshing
! automatic apertures
! get extinction cross section from beam tracing

! ############################################################################################################

! shared
real(8) start, finish ! cpu timing variables
integer(8) i

! input
character(len=*), parameter :: ifn = 'input.txt' ! input filename
character(len=255) :: cwd ! current working directory
character(len=255) :: output_dir ! output directory
integer result ! logical for navigating directories

! sr SDATIN
character(100) cfn ! crystal filename
character(100) cft ! crystal file type
character(100) afn ! apertures filename
real(8) la ! wavelength
real(8) rbi ! real part of the refractive index
real(8) ibi ! imaginary part of the refractive index
integer rec ! max number of internal beam recursions
character(100) rot_method ! rotation method
logical is_multithreaded ! whether or not code should use multithreading
integer num_orients ! number of orientations
logical intellirot ! whether or not to use intelligent euler angle choices for orientation avergaing
character(100) c_method ! method of particle file input
character(100) job_name ! name of job

! sr PDAL2
integer(8) num_vert ! number of unique vertices
integer(8) num_face !  number of faces
integer(8) num_norm ! number of face normals
integer(8), dimension(:,:), allocatable :: face_ids ! face vertex IDs
real(8), dimension(:,:), allocatable :: vert_in ! unique vertices (unrotated)
real(8), dimension(:,:), allocatable :: vert ! unique vertices (rotated)
real(8), dimension(:,:), allocatable :: norm ! face normals
integer(8), dimension(:), allocatable :: num_face_vert ! number of vertices in each face
integer(8), dimension(:), allocatable :: norm_ids ! face normal ID of each face
integer(8), dimension(:), allocatable :: apertures ! apertures asignments for each facet

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
real(8) energy_out_beam
real(8) energy_out_ext_diff

! sr diff_main
complex(8), dimension(:,:), allocatable:: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! total
real(8), dimension(:), allocatable :: theta_vals, phi_vals
complex(8), dimension(:,:), allocatable :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! total

! sr make_mueller
real(8), dimension(:,:,:), allocatable :: mueller, mueller_total, mueller_recv ! mueller matrices
real(8), dimension(:,:), allocatable :: mueller_1d, mueller_1d_total, mueller_1d_recv ! phi-integrated mueller matrices

real(8), dimension(:), allocatable :: alpha_vals, beta_vals, gamma_vals

! ############################################################################################################
! start
print*,'========== start main'
start = omp_get_wtime()
call seed(99)

job_name = read_optional_string(ifn,"jobname") ! get crystal filename
call StripSpaces(job_name) ! remove leading spaces
if(job_name(1:len(trim(job_name))) .eq. "#flagnotfound#") then
    ! print*,'jobname flag was not found, setting as my_job'
    job_name = "my_job"
end if
print*,'job name: "',job_name(1:len(trim(job_name))),'"'
call make_dir(job_name,output_dir)
call StripSpaces(output_dir)
print*,'output directory is "',trim(output_dir),'"'
open(101,file=trim(output_dir)//"/"//"log") ! open global non-standard log file for important records

! ############# input_mod #############

! stop

! read input parameters
call SDATIN(ifn,            & ! <-  input filename
            la,             & !  -> wavelength
            rbi,            & !  -> real part of the refractive index
            ibi,            & !  -> imaginary part of the refractive index
            rec,            & !  -> max number of internal beam recursions
            rot_method,     & !  -> particle rotation method
            is_multithreaded, & !  -> whether ot not code should use multithreading) 
            num_orients, &
            intellirot,  &
            c_method)

! get input particle information
call PDAL2( ifn,            & ! <-  input filename
            c_method,       & ! <-  method of particle file input
            num_vert,       & !  -> number of unique vertices
            num_face,       & !  -> number of faces
            face_ids,       & !  -> face vertex IDs
            vert_in,        & !  -> unique vertices
            num_face_vert,  & !  -> number of vertices in each face
            afn,            & ! <-  apertures filename
            apertures)

call init_loop(num_orients,alpha_vals,beta_vals,gamma_vals,intellirot)

do i = 1, num_orients  

    ! rotate particle
    call PROT_MPI(  ifn,        & ! <-  input filename
                    rot_method, & ! <-  particle rotation method
                    vert_in,    & ! <-> unique vertices (unrotated in, rotated out) to do: remove inout intent and add a rotated vertices variable
                    vert,       &
                    alpha_vals,&
                    beta_vals, &
                    gamma_vals,&
                    i, &
                    num_orients)

    ! write rotated particle to file (optional)            
    call PDAS(  vert,       & ! <-  rotated vertices
                face_ids,   & ! <-  face vertex IDs
                output_dir)   ! <-  output directory
    
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
    call beam_loop( face_ids,                  & ! <-  face vertex IDs
                    vert,                      & ! <-  unique vertices
                    la,                        & ! <-  wavelength
                    rbi,                       & ! <-  real part of the refractive index
                    ibi,                       & ! <-  imaginary part of the refractive index
                    apertures,                 & ! <-  apertures
                    rec,                       & ! <-  max number of internal beam recursions
                    beamV,                     & ! <-  beam vertices
                    beamF1,                    & ! <-  beam face vertex indices
                    beamN,                     & ! <-  beam normals
                    beamF2,                    & ! <-  beam face normal indices
                    beamMidpoints,             & ! <-  beam  midpoints
                    ampl_beam,                 & ! <-  amplitude matrix of incident beam
                    beam_outbeam_tree,         & !  -> outgoing beams from the beam tracing
                    beam_outbeam_tree_counter, & !  -> counts the current number of beam outbeams
                    ext_diff_outbeam_tree,     & !  -> outgoing beams from external diffraction
                    energy_out_beam,           & !  -> total energy out from beams (before diffraction)
                    energy_out_ext_diff)         !  -> total energy out from external diffraction (before diffraction)
    ! stop
    if(num_orients .gt. 1) then
        print'(A15,I8,A3,I8,A20,f8.4,A3)','orientation: ',i,' / ',num_orients,' (total progress: ',dble(i-1)/dble(num_orients)*100,' %)'
        ! print*,'total time elapsed: ',omp_get_wtime()-start
        ! print*,'average time per rotation: ',(omp_get_wtime()-start) / dble(i)
        if (i .gt. 1) then
            print'(A20,F12.4,A5)','est. time remaining: '
            call PROUST(nint(dble(num_orients-i+1)*(omp_get_wtime()-start) / dble(i)))
        end if
    end if

    ! diffraction
    call diff_main( beam_outbeam_tree,         & ! <-  outgoing beams from the beam tracing
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
                    is_multithreaded)            ! <-  enable or disable multithreading

    call finalise(  ampl_far_beam11,     & ! <-  amplitude matrix (1,1) due to beam diffraction
                    ampl_far_beam12,     & ! <-  amplitude matrix (1,2) due to beam diffraction
                    ampl_far_beam21,     & ! <-  amplitude matrix (2,1) due to beam diffraction
                    ampl_far_beam22,     & ! <-  amplitude matrix (2,2) due to beam diffraction
                    ampl_far_ext_diff11, & ! <-  amplitude matrix (1,1) due to external diffraction
                    ampl_far_ext_diff12, & ! <-  amplitude matrix (1,2) due to external diffraction
                    ampl_far_ext_diff21, & ! <-  amplitude matrix (2,1) due to external diffraction
                    ampl_far_ext_diff22, & ! <-  amplitude matrix (2,2) due to external diffraction
                    theta_vals,          & ! <-  theta values
                    phi_vals,            & ! <-  phi values
                    energy_out_beam,     & ! <-  total energy out from beams (before diffraction)
                    energy_out_ext_diff, & ! <-  total energy out from external diffraction (before diffraction)
                    mueller,             & !  -> 2d mueller matrix
                    mueller_1d,          & !  -> 1d mueller matrix
                    la)                    ! <-  wavelength (for optical theorem)

    ! call writeup(mueller, mueller_1d, theta_vals, phi_vals) ! write current mueller to file

    call summation(mueller, mueller_total, mueller_1d, mueller_1d_total)

end do

mueller_total = mueller_total / num_orients
mueller_1d_total = mueller_1d_total / num_orients

! writing to file
call writeup(mueller_total, mueller_1d_total, theta_vals, phi_vals, output_dir) ! write total mueller to file

finish = omp_get_wtime()

print*,'=========='
print'(A,f16.8,A)',"total time elapsed: ",finish-start," secs"
write(101,*),'======================================================'
write(101,'(A,f17.8,A)')," total time elapsed: ",finish-start," secs"
write(101,*),'======================================================'
print*,'========== end main'

close(101) ! close global non-standard output file

contains

end program main
