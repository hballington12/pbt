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
use mpi_mod

implicit none

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
integer i, j, k, i_loop

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
integer offs(1:2) ! off values, used if off rotation method
real(8) eulers(1:3) ! euler angles, used if euler rotation method
type(cc_hex_params_type) cc_hex_params ! parameters for C. Collier Gaussian Random hexagonal columns/plates
type(job_parameters_type) job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details

! sr PDAL2
integer num_vert ! number of unique vertices
integer num_face !  number of faces
integer num_norm ! number of face normals
integer, dimension(:,:), allocatable :: face_ids ! face vertex IDs
real(8), dimension(:,:), allocatable :: vert_in ! unique vertices (unrotated)
real(8), dimension(:,:), allocatable :: vert ! unique vertices (rotated)
real(8), dimension(:,:), allocatable :: norm ! face normals
integer, dimension(:), allocatable :: num_face_vert ! number of vertices in each face
integer, dimension(:), allocatable :: norm_ids ! face normal ID of each face
integer, dimension(:), allocatable :: apertures ! apertures asignments for each facet

! sr makeIncidentBeam
real(8), allocatable, dimension(:,:) :: beamV ! beam vertices
real(8), allocatable, dimension(:,:) :: beamN ! beam normals
real(8), allocatable, dimension(:,:) :: beamMidpoints ! beam  midpoints
integer, allocatable, dimension(:,:) :: beamF1 ! beam face vertex indices
integer, allocatable, dimension(:) :: beamF2 ! beam face normal indices
complex(8), allocatable, dimension(:,:,:) :: ampl_beam ! amplitude matrix of incident beam

! sr beam_loop
type(outbeamtype), dimension(:), allocatable :: beam_outbeam_tree ! outgoing beams from the beam tracing
type(outbeamtype), dimension(:), allocatable :: ext_diff_outbeam_tree ! outgoing beams from external diffraction
integer beam_outbeam_tree_counter ! counts the current number of beam outbeams
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
type(output_parameters_type) output_parameters_recv
type(output_parameters_type) output_parameters_total

real(8) max_area, max_edge_length
integer seed(1:8)
character(len=255) cache_dir ! cached files directory (if job stops early)
integer, dimension(:), allocatable :: remaining_orients
integer num_remaining_orients

! ############################################################################################################

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
tag = 1

print*,'========== start main'
start = omp_get_wtime()
seed = [0, 0, 0, 0, 0, 0, 0, 0] ! Set the seed values

print*,'my rank: ',my_rank

call parse_command_line(job_params)

if(job_params%resume) then
    print*,'attempting to resume job using cache #',job_params%cache_id
    call resume_job(job_params,num_remaining_orients,remaining_orients,mueller_total,mueller_1d_total,output_parameters_total)
    if(my_rank .ne. 0) then ! only rank 0 should keep summed parameters from the cache, reset vals for all other processes
        deallocate(mueller_total)
        deallocate(mueller_1d_total)
        ! output parameters are set to 0 in sr summation so no need to do that here
    end if
    ! stop
end if

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

max_edge_length = job_params%la*2
max_area = max_edge_length**2*sqrt(3D0)/4D0
print*,'area threshold: ',max_area

if (job_params%tri) then
    print*,'calling triangulate with max edge length: ',job_params%tri_edge_length
    call triangulate(vert_in,face_ids,num_vert,num_face,num_face_vert,job_params%tri_edge_length,'-Q -q',apertures,job_params%tri_roughness, my_rank, output_dir) ! triangulate the particle
    call merge_vertices(vert_in, face_ids, num_vert, 1D-1) ! merge vertices that are close enough
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

if (my_rank .eq. 0) then
    call RANDOM_SEED(put=seed) ! Set the seed for the random number generator
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
    allocate(alpha_vals(1:job_params%num_orients))
    allocate(beta_vals(1:job_params%num_orients))
    allocate(gamma_vals(1:job_params%num_orients))
    call MPI_RECV(alpha_vals,size(alpha_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(beta_vals,size(beta_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(gamma_vals,size(gamma_vals,1),MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierr)
end if

! some stuff
if(job_params%resume) then ! if we are resuming a cached job
    ! do nothing because we've already read this in sr resume_job
else ! if we are not resuming a cached job
    num_remaining_orients = job_params%num_orients ! number of remaining orientations is the total overall
    allocate(remaining_orients(1:num_remaining_orients)) ! allocate array to hold orientation numbers
    do i = 1, num_remaining_orients
        remaining_orients(i) = i
    end do
end if

n1 = int(num_remaining_orients / p)
n2 = mod(num_remaining_orients,  p)
my_start = my_rank*(n1+1)+1
my_end = (my_rank+1)*(n1+1)
if (my_rank .ge. n2) then
    my_start = my_start - my_rank + n2
    my_end = my_end - my_rank + n2 - 1
end if
print*,'my rank:',my_rank,'start: ',my_start,'end: ',my_end

do i = my_start, my_end

    i_loop = remaining_orients(i)

    ! rotate particle
    call PROT_MPI(  vert_in,    & ! <-> unique vertices (unrotated in, rotated out) to do: remove inout intent and add a rotated vertices variable
                    vert,       &
                    alpha_vals, &
                    beta_vals,  &
                    gamma_vals, &
                    i_loop,          &
                    job_params)

    ! fast implementation of the incident beam
    call makeIncidentBeam(  beamV,         & ! ->  beam vertices
                            beamF1,        & ! ->  beam face vertex indices
                            beamN,         & ! ->  beam normals
                            beamF2,        & ! ->  beam face normal indices
                            vert,          & ! <-  unique vertices
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
    !     print'(A15,I8,A3,I8,A20,f8.4,A3)','orientation: ',i_loop,' / ',num_orients,' (total progress: ',dble(i_loop-1)/dble(num_orients)*100,' %)'
    !     ! print*,'total time elapsed: ',omp_get_wtime()-start
    !     ! print*,'average time per rotation: ',(omp_get_wtime()-start) / dble(i_loop)
    !     if (i_loop .gt. 1) then
    !         print'(A20,F12.4,A5)','est. time remaining: '
    !         call PROUST(nint(dble(num_orients-i_loop+1)*(omp_get_wtime()-start) / dble(i_loop)))
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

    if((omp_get_wtime() - start)/3600D0 .gt. job_params%time_limit) then

        print*,'time limit reached. caching job...'

        ! rank 0 makes cache directory and sends to other processes
        if (my_rank .eq. 0) then
            call make_cache_dir("cache/",cache_dir)
            print*,'cache directory is "',trim(cache_dir),'"'
            print*,'sending...'
            do dest = 1, p-1
                CALL MPI_SEND(cache_dir, 255, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
            end do
            print*,'sending complete.'
        else
            call MPI_RECV(cache_dir,255,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,ierr)
        end if

        ! make remaining orientations file with rank 0 process
        if(my_rank .eq. 0) then
            open(unit=10,file=trim(cache_dir)//"/orient_remaining.dat",status="new")
            close(10)
        end if
        ! 1 by 1, write remaining orientations to file...
        do j = 1, p ! for each process
            if(my_rank .eq. j-1) then ! if its my turn to write to the file
                ! print*,'my_rank:',my_rank,'opening file'
                open(unit=10,file=trim(cache_dir)//"/orient_remaining.dat",status="old",access="append") ! open file in append mode
                    ! write my remaining orientations to the file
                    do k = i+1, my_end ! loop from one after my current loop index to my end index
                        write(10,*) remaining_orients(k) ! write my remaining orientations to file
                    end do
                ! print*,'my_rank:',my_rank,'wrote to file'
                close(10) ! close file
                call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! meet other processes at meet point
            else ! if its not my turn to write to file
                call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! do nothing and wait at meet point
            end if
        end do

        ! send to rank 0 process and sum
        call mpi_send_sum(  ierr,                       & ! mpi parameter
                            tag,                        & ! mpi parameter
                            p,                          & ! mpi parameter
                            status,                     & ! mpi parameter
                            my_rank,                    & ! process rank
                            mueller_1d_total,           & ! 1d mueller matrix total for each process
                            mueller_total,              & ! 2d mueller matrix total for each process
                            output_parameters_total)      ! output parameters total for each process

        ! then rank 0 writes to cached files
        if(my_rank .eq. 0) then
            call cache_job( vert_in,                    & ! unrotated vertices
                            face_ids,                   & ! face ids
                            num_face_vert,              & ! num vertices per face
                            apertures,                  & ! apertures
                            job_params,                 & ! job parameters
                            i_loop,                     & ! current loop index (unused)
                            output_parameters_total,    & ! total output parameters
                            mueller_total,              & ! total 2d mueller
                            mueller_1d_total,           & ! total 1d mueller
                            cache_dir)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        call MPI_FINALIZE(ierr)

        stop
    
    end if


end do

print*,'my rank:',my_rank,'finished'

! sum up across all mpi processes
print*,'i have finished. my rank = ',my_rank
call mpi_send_sum(  ierr,                       & ! mpi parameter
                    tag,                        & ! mpi parameter
                    p,                          & ! mpi parameter
                    status,                     & ! mpi parameter
                    my_rank,                    & ! process rank
                    mueller_1d_total,           & ! 1d mueller matrix total for each process
                    mueller_total,              & ! 2d mueller matrix total for each process
                    output_parameters_total)      ! output parameters total for each process


if (my_rank .eq. 0) then

    mueller_total = mueller_total / job_params%num_orients
    mueller_1d_total = mueller_1d_total / job_params%num_orients
    output_parameters_total%abs = output_parameters_total%abs / job_params%num_orients
    output_parameters_total%scatt = output_parameters_total%scatt / job_params%num_orients
    output_parameters_total%ext = output_parameters_total%ext / job_params%num_orients
    output_parameters_total%albedo = output_parameters_total%albedo / job_params%num_orients
    output_parameters_total%asymmetry = output_parameters_total%asymmetry / job_params%num_orients
    output_parameters_total%abs_eff = output_parameters_total%abs_eff / job_params%num_orients
    output_parameters_total%scatt_eff = output_parameters_total%scatt_eff / job_params%num_orients
    output_parameters_total%ext_eff = output_parameters_total%ext_eff / job_params%num_orients
    output_parameters_total%geo_cross_sec = output_parameters_total%geo_cross_sec / job_params%num_orients
    
    ! writing to file
    call write_outbins(output_dir,job_params%theta_vals,job_params%phi_vals)
    call writeup(mueller_total, mueller_1d_total, output_dir, output_parameters_total, job_params) ! write to file

    ! clean up temporary files
    call system("rm -r "//trim(output_dir)//"/tmp") ! remove directory for temp files

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
