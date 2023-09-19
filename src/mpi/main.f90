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
! use mpi

implicit none

include 'mpif.h'

! to do:
! add vector re-normalisation checks after rotations (ie. sr rotate_into_aperture_system)
! fix direct forward scattering
! mpi support
! asymmetry parameter, single scat albedo
! add absorption
! add quad support
! add support to avoid crash if nan detected
! automatic meshing
! automatic apertures
! fix phi = 0 numerical error

! ############################################################################################################

! shared
real(8) start, finish ! cpu timing variables
integer(8) i

! input
character(len=*), parameter :: ifn = 'input.txt' ! input filename
character(len=255) :: output_dir ! output directory
character(len=255) :: my_rank_str ! string of my rank
character(len=255) :: my_log_dir ! log file location for each mpi process
integer result ! true if subdirectory was made, false if subdirectory was not made

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
integer(8) offs(1:2) ! off values, used if off rotation method
real(8) eulers(1:3) ! euler angles, used if euler rotation method
type(cc_hex_params_type) cc_hex_params ! parameters for C. Collier Gaussian Random hexagonal columns/plates
type(job_parameters_type) job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details

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
real(8) energy_abs_beam

! sr diff_main
complex(8), dimension(:,:), allocatable:: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! total
real(8), dimension(:), allocatable :: theta_vals, phi_vals
complex(8), dimension(:,:), allocatable :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! total

! sr make_mueller
real(8), dimension(:,:,:), allocatable :: mueller, mueller_total, mueller_recv ! mueller matrices
real(8), dimension(:,:), allocatable :: mueller_1d, mueller_1d_total, mueller_1d_recv ! phi-integrated mueller matrices

! mpi
integer ierr
integer tag
integer p
integer source, dest
integer status(MPI_STATUS_SIZE)
integer my_start, my_end, my_rank
integer n1, n2
real(8), dimension(:), allocatable :: alpha_vals, beta_vals, gamma_vals

! sr finalise
type(output_parameters_type) output_parameters 
type(output_parameters_type) output_parameters_total

real(8) max_area, max_edge_length

! ############################################################################################################

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
tag = 1

print*,'========== start main'
start = omp_get_wtime()
! call seed(99) ! find gfotran replacement for this later!

print*,'my rank: ',my_rank

call parse_command_line(job_params)

! setting up job directory
! rank 0 process broadcasts job directory to other processes
if (my_rank .eq. 0) then
    call make_dir(job_params%job_name,output_dir)
    print*,'output directory is "',trim(output_dir),'"'
    ! result = makedirqq(trim(output_dir)//"/logs") ! make directory for logs
    call system("mkdir "//trim(output_dir)//"/logs")
    print*,'log files directory is "',trim(output_dir)//"/logs/",'"'
    call system("mkdir "//trim(output_dir)//"/tmp") ! make directory for temp files
    print*,'temporary files directory is "',trim(output_dir)//"/tmp/",'"'
    print*,'sending...'
    print*,'output_dir: ',output_dir
    do dest = 1, p-1
        CALL MPI_SEND(output_dir, 255, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
    end do
    ! print*,'sending complete.'
else
    call MPI_RECV(output_dir,255,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,ierr)
end if

write(my_rank_str,*) my_rank
call StripSpaces(my_rank_str)
print*,'my_rank_str: ',my_rank_str
write(my_log_dir,*) trim(output_dir)//"/logs/log",trim(my_rank_str)
my_log_dir = adjustl(my_log_dir)
! print*,'my rank is : ,',my_rank,' , my output directory is "',trim(my_log_dir),'"'

open(101,file=trim(my_log_dir),status="new") ! open global non-standard log file for important records

! get input particle information
call PDAL2( num_vert,       & !  -> number of unique vertices
            num_face,       & !  -> number of faces
            face_ids,       & !  -> face vertex IDs
            vert_in,        & !  -> unique vertices
            num_face_vert,  & !  -> number of vertices in each face
            apertures,      &
            job_params)

num_orients = job_params%num_orients 
n1 = int(num_orients / p)
n2 = mod(num_orients,  p)
my_start = my_rank*(n1+1)+1
my_end = (my_rank+1)*(n1+1)
if (my_rank .ge. n2) then
    my_start = my_start - my_rank + n2
    my_end = my_end - my_rank + n2 - 1
end if
print*,'my rank:',my_rank,'start: ',my_start,'end: ',my_end

allocate(alpha_vals(1:num_orients))
allocate(beta_vals(1:num_orients))
allocate(gamma_vals(1:num_orients))

if (my_rank .eq. 0) then
    call init_loop( alpha_vals, &
                    beta_vals, &
                    gamma_vals, &
                    job_params)
    do dest = 1, p-1
        call MPI_SEND(alpha_vals,size(alpha_vals,1),MPI_REAL8,dest,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(beta_vals,size(beta_vals,1),MPI_REAL8,dest,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(gamma_vals,size(gamma_vals,1),MPI_REAL8,dest,tag,MPI_COMM_WORLD,ierr)
    end do
else
    call MPI_RECV(alpha_vals,size(alpha_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(beta_vals,size(beta_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(gamma_vals,size(gamma_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
end if

max_edge_length = job_params%la*2
max_area = max_edge_length**2*sqrt(3D0)/4D0
print*,'area threshold: ',max_area

if (job_params%tri) then
    print*,'calling triangulate with max edge length: ',job_params%tri_edge_length
    call triangulate(vert_in,face_ids,num_vert,num_face,num_face_vert,job_params%tri_edge_length,'-Q -q',apertures,job_params%tri_roughness, my_rank, output_dir) ! triangulate the particle
    call merge_vertices(vert_in, face_ids, num_vert, num_face, 1D-1) ! merge vertices that are close enough
    call fix_collinear_vertices(vert_in, face_ids, num_vert, num_face, num_face_vert, apertures)
    ! call triangulate(vert_in,face_ids,num_vert,num_face,num_face_vert,max_area,'-Q -q',apertures,0D0) ! retriangulate the particle to have no area greater than 10*lambda
end if

if (my_rank .eq. 0) then
    ! write unrotated particle to file (optional)            
    call PDAS(  vert_in,        & ! <-  rotated vertices
                face_ids,       & ! <-  face vertex IDs
                output_dir,     & ! <-  output directory
                num_face_vert,  & ! <-  number of verices in each face
                "unrotated")    ! <-  filename
end if

do i = my_start, my_end  

    ! rotate particle
    call PROT_MPI(  vert_in,    & ! <-> unique vertices (unrotated in, rotated out) to do: remove inout intent and add a rotated vertices variable
                    vert,       &
                    alpha_vals, &
                    beta_vals,  &
                    gamma_vals, &
                    i,          &
                    job_params)

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
                        apertures,                 & ! <-  apertures
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
                        energy_out_ext_diff,       & !  -> total energy out from external diffraction (before diffraction)
                        energy_abs_beam,           & !  -> total energy absorbed from beams (before diffraction)
                        output_parameters,         & !  -> adds illuminated geometric cross section to output parameters
                        num_face_vert,             & ! <-  number of verices in each face
                        job_params)

    ! if(num_orients .gt. 1) then
    !     print'(A15,I8,A3,I8,A20,f8.4,A3)','orientation: ',i,' / ',num_orients,' (total progress: ',dble(i-1)/dble(num_orients)*100,' %)'
    !     ! print*,'total time elapsed: ',omp_get_wtime()-start
    !     ! print*,'average time per rotation: ',(omp_get_wtime()-start) / dble(i)
    !     if (i .gt. 1) then
    !         print'(A20,F12.4,A5)','est. time remaining: '
    !         call PROUST(nint(dble(num_orients-i+1)*(omp_get_wtime()-start) / dble(i)))
    !     end if
    ! end if
    

    ! diffraction
    call diff_main( beam_outbeam_tree,         & ! <-  outgoing beams from the beam tracing
                    beam_outbeam_tree_counter, & ! <-  counts the current number of beam outbeams
                    ampl_far_beam11,           & !  -> amplitude matrix (1,1) due to beam diffraction
                    ampl_far_beam12,           & !  -> amplitude matrix (1,2) due to beam diffraction
                    ampl_far_beam21,           & !  -> amplitude matrix (2,1) due to beam diffraction
                    ampl_far_beam22,           & !  -> amplitude matrix (2,2) due to beam diffraction
                    ext_diff_outbeam_tree,     & ! <-  outgoing beams from external diffraction
                    ampl_far_ext_diff11,       & !  -> amplitude matrix (1,1) due to external diffraction
                    ampl_far_ext_diff12,       & !  -> amplitude matrix (1,2) due to external diffraction
                    ampl_far_ext_diff21,       & !  -> amplitude matrix (2,1) due to external diffraction
                    ampl_far_ext_diff22,       & !  -> amplitude matrix (2,2) due to external diffraction
                    job_params)

    call finalise(  ampl_far_beam11,     & ! <-  amplitude matrix (1,1) due to beam diffraction
                    ampl_far_beam12,     & ! <-  amplitude matrix (1,2) due to beam diffraction
                    ampl_far_beam21,     & ! <-  amplitude matrix (2,1) due to beam diffraction
                    ampl_far_beam22,     & ! <-  amplitude matrix (2,2) due to beam diffraction
                    ampl_far_ext_diff11, & ! <-  amplitude matrix (1,1) due to external diffraction
                    ampl_far_ext_diff12, & ! <-  amplitude matrix (1,2) due to external diffraction
                    ampl_far_ext_diff21, & ! <-  amplitude matrix (2,1) due to external diffraction
                    ampl_far_ext_diff22, & ! <-  amplitude matrix (2,2) due to external diffraction
                    energy_out_beam,     & ! <-  total energy out from beams (before diffraction)
                    energy_out_ext_diff, & ! <-  total energy out from external diffraction (before diffraction)
                    mueller,             & !  -> 2d mueller matrix
                    mueller_1d,          & !  -> 1d mueller matrix
                    energy_abs_beam,     & ! <-  energy absorbed within the particle
                    output_parameters,   & !  -> some output parameters
                    job_params)

    ! call writeup(mueller, mueller_1d, theta_vals, phi_vals) ! write current mueller to file

    call summation(mueller, mueller_total, mueller_1d, mueller_1d_total,output_parameters,output_parameters_total)

end do

print*,'my rank:',my_rank,'finished'

! sum up mueller 2d and 1d across all mpi processes
print*,'i have finished. my rank = ',my_rank
if (my_rank .ne. 0) then ! if not rank 0 process, send mueller to rank 0
    call MPI_SEND(mueller_1d_total,size(mueller_1d_total,1)*size(mueller_1d_total,2),MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(mueller_total,size(mueller_total,1)*size(mueller_total,2)*size(mueller_total,3),MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
    print*,'sent to rank 0. my rank = ',my_rank
else ! if rank 0 process, receieve from all other ranks
    ! allocate some arrays to hold the received values
    allocate(mueller_1d_recv(1:size(mueller_1d_total,1),1:size(mueller_1d_total,2)))
    allocate(mueller_recv(1:size(mueller_total,1),1:size(mueller_total,2),1:size(mueller_total,3)))
    do source = 1,p-1 ! collect from other processes
        ! print*,'attempting to reveive from ',source,' my rank = ',my_rank
        call MPI_RECV(mueller_1d_recv,size(mueller_1d_recv,1)*size(mueller_1d_recv,2),MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
        mueller_1d_total = mueller_1d_total + mueller_1d_recv ! sum
        call MPI_RECV(mueller_recv,size(mueller_recv,1)*size(mueller_recv,2)*size(mueller_recv,3),MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
        mueller_total = mueller_total + mueller_recv ! sum        
        print*,'received from ',source,' my rank = ',my_rank
    end do
end if


if (my_rank .eq. 0) then

    mueller_total = mueller_total / num_orients
    mueller_1d_total = mueller_1d_total / num_orients
    output_parameters_total%abs = output_parameters_total%abs / num_orients
    output_parameters_total%scatt = output_parameters_total%scatt / num_orients
    output_parameters_total%ext = output_parameters_total%ext / num_orients
    output_parameters_total%albedo = output_parameters_total%albedo / num_orients
    output_parameters_total%asymmetry = output_parameters_total%asymmetry / num_orients
    output_parameters_total%abs_eff = output_parameters_total%abs_eff / num_orients
    output_parameters_total%scatt_eff = output_parameters_total%scatt_eff / num_orients
    output_parameters_total%ext_eff = output_parameters_total%ext_eff / num_orients
    output_parameters_total%geo_cross_sec = output_parameters_total%geo_cross_sec / num_orients
    
    ! writing to file
    call write_outbins(output_dir,job_params%theta_vals,job_params%phi_vals)
    call writeup(mueller_total, mueller_1d_total, output_dir, output_parameters_total, job_params) ! write to file

    ! clean up temporary files
    call system("rm -r "//trim(output_dir)//"/tmp") ! make directory for temp files

    finish = omp_get_wtime()
    print*,'=========='
    print'(A,f16.8,A)',"total time elapsed: ",finish-start," secs"
    write(101,*)'======================================================'
    write(101,'(A,f17.8,A)')" total time elapsed: ",finish-start," secs"
    write(101,*)'======================================================'
    print*,'========== end main'
end if

close(101) ! close global non-standard output file

call MPI_FINALIZE(ierr)

contains

end program main
