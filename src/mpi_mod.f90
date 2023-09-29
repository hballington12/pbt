! abt_mpi_mod.f90
! module containing various subroutines related to MPI
! also loads the mpif.h file

module mpi_mod

use types_mod

implicit none

include 'mpif.h'

contains

subroutine mpi_send_sum(ierr,                       & ! mpi parameter
                        tag,                        & ! mpi parameter
                        p,                          & ! mpi parameter
                        status,                     & ! mpi parameter
                        my_rank,                    & ! process rank
                        mueller_1d_total,           & ! 1d mueller matrix total for each process
                        mueller_total,              & ! 2d mueller matrix total for each process
                        output_parameters_total)      ! output parameters total for each process

    ! sends all stuff to rank 0 process, which then sums
 
    real(8), dimension(:,:), allocatable, intent(inout) :: mueller_1d_total ! phi-integrated mueller matrices
    real(8), dimension(:,:,:), allocatable, intent(inout) :: mueller_total ! mueller matrices
    type(output_parameters_type), intent(inout) :: output_parameters_total
    integer, intent(in) :: ierr
    integer, intent(in) :: tag
    integer, intent(in) :: p
    integer, intent(in) :: status(MPI_STATUS_SIZE)
    integer, intent(in) :: my_rank

    real(8), dimension(:,:), allocatable :: mueller_1d_recv ! phi-integrated mueller matrices
    real(8), dimension(:,:,:), allocatable :: mueller_recv ! mueller matrices
    type(output_parameters_type) output_parameters_recv
    integer source, dest

    if (my_rank .ne. 0) then ! if not rank 0 process, send mueller to rank 0
        call MPI_SEND(mueller_1d_total,size(mueller_1d_total,1)*size(mueller_1d_total,2),MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(mueller_total,size(mueller_total,1)*size(mueller_total,2)*size(mueller_total,3),MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%abs,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%scatt,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%ext,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%albedo,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%asymmetry,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%abs_eff,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%scatt_eff,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%ext_eff,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(output_parameters_total%geo_cross_sec,1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
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
            call MPI_RECV(output_parameters_recv%abs,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%abs = output_parameters_total%abs + output_parameters_recv%abs ! sum  
            call MPI_RECV(output_parameters_recv%scatt,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%scatt = output_parameters_total%scatt + output_parameters_recv%scatt ! sum  
            call MPI_RECV(output_parameters_recv%ext,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%ext = output_parameters_total%ext + output_parameters_recv%ext ! sum  
            call MPI_RECV(output_parameters_recv%albedo,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%albedo = output_parameters_total%albedo + output_parameters_recv%albedo ! sum  
            call MPI_RECV(output_parameters_recv%asymmetry,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%asymmetry = output_parameters_total%asymmetry + output_parameters_recv%asymmetry ! sum  
            call MPI_RECV(output_parameters_recv%abs_eff,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%abs_eff = output_parameters_total%abs_eff + output_parameters_recv%abs_eff ! sum  
            call MPI_RECV(output_parameters_recv%scatt_eff,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%scatt_eff = output_parameters_total%scatt_eff + output_parameters_recv%scatt_eff ! sum  
            call MPI_RECV(output_parameters_recv%ext_eff,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%ext_eff = output_parameters_total%ext_eff + output_parameters_recv%ext_eff ! sum  
            call MPI_RECV(output_parameters_recv%geo_cross_sec,1,MPI_REAL8,source,tag,MPI_COMM_WORLD,status,ierr)
            output_parameters_total%geo_cross_sec = output_parameters_total%geo_cross_sec + output_parameters_recv%geo_cross_sec ! sum  
            print*,'received from ',source,' my rank = ',my_rank
        end do
    end if

end subroutine

end module