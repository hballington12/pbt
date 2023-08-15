
! input_mod.f90
! module containing various subroutines used to initiate aperture beam tracer

module input_mod

use misc_submod
use types_mod
use cc_hex_mod

implicit none

contains

subroutine area_stats(face_areas)

real(8), dimension(:), allocatable, intent(in) :: face_areas ! area of each facet
integer i, num_facets
real(8) min_area, max_area, total_area, avg_area

! init
min_area = face_areas(1)
max_area = face_areas(1)
total_area = 0

num_facets = size(face_areas,1)

do i = 1, num_facets

    total_area = total_area + face_areas(i)
    if(face_areas(i) .gt. max_area) max_area = face_areas(i)
    if(face_areas(i) .lt. min_area) min_area = face_areas(i)

end do

avg_area = total_area / num_facets

print*,'max facet area:      ',max_area
print*,'min facet area:      ',min_area
print*,'total surface area: ',total_area
print*,'avg. facet area:     ',avg_area

write(101,*),'max facet area:      ',max_area
write(101,*),'min facet area:      ',min_area
write(101,*),'total surface area: ',total_area
write(101,*),'avg. facet area:     ',avg_area

end subroutine

subroutine init_loop(num_orients,alpha_vals,beta_vals,gamma_vals,intellirot)

    ! subroutine to pick angles that will be looped over

    integer, intent(in) :: num_orients ! number of orientations
    logical, intent(in) :: intellirot ! whether or not to use intelligent euler angle choices for orientation avergaing
    real(8), dimension(:), allocatable, intent(out) :: alpha_vals, beta_vals, gamma_vals
    integer(8) num_angles, leftover_angles
    real(8), allocatable, dimension(:) :: intelli_vals
    integer(8) i, j, k, counter
    real(8) spacing, rand

    allocate(alpha_vals(1:num_orients))
    allocate(beta_vals(1:num_orients))
    allocate(gamma_vals(1:num_orients))

    if (intellirot) then

        print*,'number of orientations:',num_orients
        num_angles = floor(num_orients**(1d0/3d0))
        print*,'number of intelligent euler angles:',num_angles
        leftover_angles = num_orients - num_angles**3
        print*,'number of (leftover) random euler angles: ',leftover_angles

        allocate(intelli_vals(1:num_angles)) ! allocate array to hold intelligent euler angles
        if(num_angles .eq. 1) then
            intelli_vals(1) = 0.5 ! set to middle if less than than 8 orientations specified
        else
            spacing = 1d0/(num_angles-1)
            do i = 1, num_angles ! for each entry, linear interpolate from 0 to 1
                intelli_vals(i) = spacing * (i-1)
                ! print*,'intelli_vals(i)',intelli_vals(i)
            end do
        end if

        ! loop through and assign intelligent angles
        counter = 0
        do i = 1, num_angles
            do j = 1,num_angles
                do k = 1,num_angles
                    counter = counter + 1 ! count how many orientations we have set up
                    alpha_vals(counter) = intelli_vals(i)
                    beta_vals(counter) = intelli_vals(j)
                    gamma_vals(counter) = intelli_vals(k)

                end do
            end do
        end do

        ! fix numerical errors with a very small amount of extra rotation
        alpha_vals(1:counter) = abs(alpha_vals(1:counter) - 0.0001)
        beta_vals(1:counter) = abs(beta_vals(1:counter) - 0.0001)
        gamma_vals(1:counter) = abs(gamma_vals(1:counter) - 0.0001)

        ! fill in remainining angles with random numbers
        do i = 1, leftover_angles
            counter = counter + 1
            call random_number(rand)
            alpha_vals(counter) = rand
            call random_number(rand)
            beta_vals(counter) = rand
            call random_number(rand)
            gamma_vals(counter) = rand                        
        end do

        ! ! print intelligent euler angles
        ! do i = 1, num_orients
        !     print'(A,f6.4,A,f6.4,A,f6.4)','alpha: ',alpha_vals(i),' beta: ',beta_vals(i),' gamma: ',gamma_vals(i)
        ! end do

    else
        do i = 1, size(alpha_vals,1) ! loop here so that the angles are reproducable regardless of number of orientations
            call random_number(alpha_vals(i))
            call random_number(beta_vals(i))
            call random_number(gamma_vals(i))
        end do
    end if

    end subroutine

subroutine PROT(ifn,rot_method,verts)

    ! rotates particle

    character(len=*), intent(in) :: ifn
    character(100), intent(in) :: rot_method ! rotation method
    real(8), dimension(:,:), allocatable, intent(inout) :: verts ! unique vertices    
    integer(8) offs(1:2)
    real(8) eulers(1:3)
    real(8) vec(1:3) ! off rotation vector
    real(8) hilf0, hilf1
    real(8) rot1(1:3,1:3), rot2(1:3,1:3), rot(1:3,1:3)
    integer i
    real(8) s1, s2, s3, c1, c2, c3
    real(8) rand

    print*,'========== start sr PROT'

    print*,'rotation method: "',rot_method(1:len(trim(rot_method))),'"'

    if(rot_method(1:len(trim(rot_method))) .eq. 'none') then
        ! do nothing
    else if(rot_method(1:len(trim(rot_method))) .eq. 'off') then
        call read_input_vals(ifn,"rot off",offs,2)
        ! print*,'off values: ', offs(1:2)

        if(offs(1) .eq. 30 .and. offs(2) .eq. 0) then
            print*,'off setting: 30x0'
            vec = (/-0.866025,-0.5,0.0/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 10) then
            print*,'off setting: 30x10'
            vec = (/-0.866025,-0.492404,0.0868240/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 20) then
            print*,'off setting: 30x20'
            vec = (/-0.866025,-0.469846,0.171010/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 30) then
            print*,'off setting: 30x30'
            vec = (/-0.866025,-0.433013,0.25/)
        end if

        hilf0 = vec(2)/vec(1)
        hilf1 = cos(atan(hilf0))

        rot1(1,1) = hilf1
        rot1(1,2) = cos(pi/2-atan(hilf0))
        rot1(1,3) = 0
        rot1(2,1) = cos(pi/2+atan(hilf0))
        rot1(2,2) = hilf1
        rot1(2,3) = 0
        rot1(3,1) = 0
        rot1(3,2) = 0
        rot1(3,3) = 1

        rot2(1,1) = cos(pi-acos(vec(3)))
        rot2(1,2) = 0
        rot2(1,3) = cos(pi/2+acos(vec(3)))
        rot2(2,1) = 0
        rot2(2,2) = 1
        rot2(2,3) = 0
        rot2(3,1) = cos(pi/2-acos(vec(3)))
        rot2(3,2) = 0
        rot2(3,3) = cos(pi-acos(vec(3)))

        rot = matmul(rot2,rot1)

        do i = 1, size(verts,1) ! for each vertex
            verts(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do

        ! print*,'verts(2,1:3)',verts(2,1:3)

    else if(rot_method(1:len(trim(rot_method))) .eq. 'euler') then
        call read_input_vals_real(ifn,"rot euler",eulers,3)
        print*,'alpha:',eulers(1)
        print*,'beta:',eulers(2)
        print*,'gamma:',eulers(3)

        eulers = eulers*pi/180 ! convert to rad

        ! mishchenko rotation
        s1 = sin(eulers(1))
        s2 = sin(eulers(2))
        s3 = sin(eulers(3))
        c1 = cos(eulers(1))
        c2 = cos(eulers(2))
        c3 = cos(eulers(3))

        ! make rotation matrix
        rot(1,1) = c1*c2*c3 - s1*s3
        rot(1,2) = -c1*c2*s3 - s1*c3
        rot(1,3) = c1*s2
        rot(2,1) = s1*c2*c3 + c1*s3
        rot(2,2) = -s1*c2*s3 + c1*c3
        rot(2,3) = s1*s2
        rot(3,1) = -s2*c3
        rot(3,2) = s2*s3
        rot(3,3) = c2

        do i = 1, size(verts,1) ! for each vertex
            verts(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do
        
    else if(rot_method(1:len(trim(rot_method))) .eq. 'multi') then
        call random_number(rand)
        eulers(1) = 2*pi*(rand)

        call random_number(rand)
        eulers(2) = acos(1.0 - 2.0*rand)

        call random_number(rand)
        eulers(3) = 2*pi*(rand)

        print*,'alpha:',eulers(1)*180/pi
        print*,'beta:',eulers(2)*180/pi
        print*,'gamma:',eulers(3)*180/pi

        ! mishchenko rotation
        s1 = sin(eulers(1))
        s2 = sin(eulers(2))
        s3 = sin(eulers(3))
        c1 = cos(eulers(1))
        c2 = cos(eulers(2))
        c3 = cos(eulers(3))

        ! make rotation matrix
        rot(1,1) = c1*c2*c3 - s1*s3
        rot(1,2) = -c1*c2*s3 - s1*c3
        rot(1,3) = c1*s2
        rot(2,1) = s1*c2*c3 + c1*s3
        rot(2,2) = -s1*c2*s3 + c1*c3
        rot(2,3) = s1*s2
        rot(3,1) = -s2*c3
        rot(3,2) = s2*s3
        rot(3,3) = c2

        do i = 1, size(verts,1) ! for each vertex
            verts(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do
    end if

    ! stop

    print*,'========== end sr PROT'

end subroutine

subroutine PROT_MPI(ifn,                & ! input filename (so we can read arguments from the "rot" line)
                    rot_method,         & ! rotation method: "euler", "off", "none", "multi"
                    verts,              & ! unrotated vertices
                    verts_rot,          & ! rotated vertices
                    alpha_vals,         & ! list of values for euler alpha angle in range 0 to 1 (for multirot only)
                    beta_vals,          & ! list of values for beta alpha angle in range 0 to 1 (for multirot only)
                    gamma_vals,         & ! list of values for gamma alpha angle in range 0 to 1 (for multirot only)
                    loop_index,         & ! current loop index that each process is on
                    num_orients)          ! total number of orientations for this process

    ! rotates particle
    ! modified to ensure different mpi processes have different 

    character(len=*), intent(in) :: ifn
    character(100), intent(in) :: rot_method ! rotation method
    real(8), dimension(:,:), allocatable, intent(inout) :: verts ! unique vertices    
    real(8), dimension(:,:), allocatable, intent(out) :: verts_rot ! unique vertices    
    real(8), dimension(:), allocatable, intent(in) :: alpha_vals, beta_vals, gamma_vals
    integer(8) offs(1:2)
    real(8) eulers(1:3)
    real(8) vec(1:3) ! off rotation vector
    real(8) hilf0, hilf1
    real(8) rot1(1:3,1:3), rot2(1:3,1:3), rot(1:3,1:3)
    integer i
    real(8) s1, s2, s3, c1, c2, c3
    real(8) rand
    integer(8), intent(in) :: loop_index
    integer, intent(in) ::  num_orients ! number of orientations

    print*,'========== start sr PROT_MPI'
    write(101,*),'======================================================'
    write(101,*),'======================================================'
    write(101,*),'orientation: ',loop_index,' / ',num_orients

    print*,'rotation method: "',rot_method(1:len(trim(rot_method))),'"'

    allocate(verts_rot(1:size(verts,1),1:size(verts,2))) ! allocate array for rotated vertices

    if(rot_method(1:len(trim(rot_method))) .eq. 'none') then
        ! do nothing
    else if(rot_method(1:len(trim(rot_method))) .eq. 'off') then
        call read_input_vals(ifn,"rot off",offs,2)
        ! print*,'off values: ', offs(1:2)

        if(offs(1) .eq. 30 .and. offs(2) .eq. 0) then
            print*,'off setting: 30x0'
            write(101,*),'off setting: "30x0"'
            vec = (/-0.866025,-0.5,0.0/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 10) then
            print*,'off setting: 30x10'
            write(101,*),'off setting: "30x10"'
            vec = (/-0.866025,-0.492404,0.0868240/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 20) then
            print*,'off setting: 30x20'
            write(101,*),'off setting: "30x20"'
            vec = (/-0.866025,-0.469846,0.171010/)
        else if(offs(1) .eq. 30 .and. offs(2) .eq. 30) then
            print*,'off setting: 30x30'
            write(101,*),'off setting: "30x30"'
            vec = (/-0.866025,-0.433013,0.25/)
        end if

        hilf0 = vec(2)/vec(1)
        hilf1 = cos(atan(hilf0))

        rot1(1,1) = hilf1
        rot1(1,2) = cos(pi/2-atan(hilf0))
        rot1(1,3) = 0.0
        rot1(2,1) = cos(pi/2+atan(hilf0))
        rot1(2,2) = hilf1
        rot1(2,3) = 0.0
        rot1(3,1) = 0.0
        rot1(3,2) = 0.0
        rot1(3,3) = 1.0

        rot2(1,1) = cos(pi-acos(vec(3)))
        rot2(1,2) = 0.0
        rot2(1,3) = cos(pi/2+acos(vec(3)))
        rot2(2,1) = 0.0
        rot2(2,2) = 1.0
        rot2(2,3) = 0.0
        rot2(3,1) = cos(pi/2-acos(vec(3)))
        rot2(3,2) = 0.0
        rot2(3,3) = cos(pi-acos(vec(3)))

        rot = matmul(rot2,rot1)

        do i = 1, size(verts,1) ! for each vertex
            verts_rot(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do

        ! print*,'verts(2,1:3)',verts(2,1:3)

    else if(rot_method(1:len(trim(rot_method))) .eq. 'euler') then
        call read_input_vals_real(ifn,"rot euler",eulers,3)
        print*,'alpha:',eulers(1)
        write(101,*),'alpha:',eulers(1)
        print*,'beta: ',eulers(2)
        write(101,*),'beta: ',eulers(2)
        print*,'gamma:',eulers(3)
        write(101,*),'gamma:',eulers(3)

        eulers = eulers*pi/180.0 ! convert to rad

        ! mishchenko rotation
        s1 = sin(eulers(1))
        s2 = sin(eulers(2))
        s3 = sin(eulers(3))
        c1 = cos(eulers(1))
        c2 = cos(eulers(2))
        c3 = cos(eulers(3))

        ! make rotation matrix
        rot(1,1) = c1*c2*c3 - s1*s3
        rot(1,2) = -c1*c2*s3 - s1*c3
        rot(1,3) = c1*s2
        rot(2,1) = s1*c2*c3 + c1*s3
        rot(2,2) = -s1*c2*s3 + c1*c3
        rot(2,3) = s1*s2
        rot(3,1) = -s2*c3
        rot(3,2) = s2*s3
        rot(3,3) = c2

        do i = 1, size(verts,1) ! for each vertex
            verts_rot(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do
        
    else if(rot_method(1:len(trim(rot_method))) .eq. 'multi') then

        ! print*,'rot method was multi'

        rand = alpha_vals(loop_index)    
        eulers(1) = 2*pi*(rand)

        rand = beta_vals(loop_index) 
        eulers(2) = acos(1.0 - 2.0*rand)

        rand = gamma_vals(loop_index) 
        eulers(3) = 2*pi*(rand)

        print*,'alpha:',eulers(1)*180.0/pi
        write(101,*),'alpha:',eulers(1)*180.0/pi
        print*,'beta: ',eulers(2)*180.0/pi
        write(101,*),'beta: ',eulers(2)*180.0/pi
        print*,'gamma:',eulers(3)*180.0/pi
        write(101,*),'gamma:',eulers(3)*180.0/pi

        ! mishchenko rotation
        s1 = sin(eulers(1))
        s2 = sin(eulers(2))
        s3 = sin(eulers(3))
        c1 = cos(eulers(1))
        c2 = cos(eulers(2))
        c3 = cos(eulers(3))

        ! make rotation matrix
        rot(1,1) = c1*c2*c3 - s1*s3
        rot(1,2) = -c1*c2*s3 - s1*c3
        rot(1,3) = c1*s2
        rot(2,1) = s1*c2*c3 + c1*s3
        rot(2,2) = -s1*c2*s3 + c1*c3
        rot(2,3) = s1*s2
        rot(3,1) = -s2*c3
        rot(3,2) = s2*s3
        rot(3,3) = c2

        do i = 1, size(verts,1) ! for each vertex
            verts_rot(i,1:3) = matmul(rot,verts(i,1:3)) ! rotate
        end do
    end if

    ! stop

    print*,'========== end sr PROT_MPI'

end subroutine

subroutine read_input_vals(fn,string,vals,num_vals)

    ! read_input_vals reads values from a line in a text file
    ! a string containing the first part of the line is used to determine which line should be read from
    ! the number of values read from the line is determined by an input integer
    ! to do: add custom delimiter

    character(len=*), intent(in) :: fn ! filename to read from
    character(len=*), intent(in) :: string ! first part of the line to read from
    integer(8), dimension(:) :: vals ! values to be read
    integer(8) num_vals ! number of values to read

    integer(8), parameter :: max_line_length = 150 ! max number of characters in a line of thecrystal file (might need increasing if faces have many vertices)
    character(max_line_length) line ! a line in a file
    integer(8) i, io ! counting variables
    integer(8) num_lines
    logical success
    integer j, k, val_count ! number chars since last delimiter

    open(10,file = fn, status = 'old')

    num_lines = 0  ! initialise line counter
    success = .false.
    val_count = 0

    do
        read(10,*,iostat=io)
        if (io/=0) exit
        num_lines = num_lines + 1
    end do
    
    rewind(10) ! back to top of file
    
    do i = 1,num_lines
        read(10,"(a)",iostat=io)line
        if (line(1:len(string)+1) .eq. string//' ') then
            ! print*,'match found on line:',i
            k = 0
            do j = len(string)+2, len(line) ! read through characters in line
                if(line(j:j) .ne. ' ') then ! if character wasnt space delimiter
                    k = k + 1 ! update characters counted since last space delimiter
                    ! print*,j,line(j:j)
                else ! if delimiter found
                    if(k .gt. 0 .and. val_count .lt. num_vals) then
                        val_count = val_count + 1
                        read(line(j-k:j),*) vals(val_count)
                        ! print*,'value found: ',vals(val_count)
                    end if
                    k = 0
                end if
            end do

        end if
    
    end do

    close(10)

end subroutine

subroutine read_input_vals_real(fn,string,vals,num_vals)

    ! read_input_vals reads values from a line in a text file
    ! a string containing the first part of the line is used to determine which line should be read from
    ! the number of values read from the line is determined by an input integer
    ! to do: add custom delimiter

    character(len=*), intent(in) :: fn ! filename to read from
    character(len=*), intent(in) :: string ! first part of the line to read from
    real(8), dimension(:) :: vals ! values to be read
    integer(8) num_vals ! number of values to read

    integer(8), parameter :: max_line_length = 150 ! max number of characters in a line of thecrystal file (might need increasing if faces have many vertices)
    character(max_line_length) line ! a line in a file
    integer(8) i, io ! counting variables
    integer(8) num_lines
    logical success
    integer j, k, val_count ! number chars since last delimiter

    open(10,file = fn, status = 'old')

    num_lines = 0  ! initialise line counter
    success = .false.
    val_count = 0

    do
        read(10,*,iostat=io)
        if (io/=0) exit
        num_lines = num_lines + 1
    end do
    
    rewind(10) ! back to top of file
    
    do i = 1,num_lines
        read(10,"(a)",iostat=io)line
        if (line(1:len(string)+1) .eq. string//' ') then
            ! print*,'match found on line:',i
            k = 0
            do j = len(string)+2, len(line) ! read through characters in line
                if(line(j:j) .ne. ' ') then ! if character wasnt space delimiter
                    k = k + 1 ! update characters counted since last space delimiter
                    ! print*,j,line(j:j)
                else ! if delimiter found
                    if(k .gt. 0 .and. val_count .lt. num_vals) then
                        val_count = val_count + 1
                        read(line(j-k:j),*) vals(val_count)
                        ! print*,'value found: ',vals(val_count)
                    end if
                    k = 0
                end if
            end do

        end if
    
    end do

    close(10)

end subroutine

subroutine readApertures(afn,apertures, face_ids)

! sr readApertures reads the apertures filename
! will exit if the apertures file does not have the correct number of lines

integer(8), dimension(:), allocatable, intent(out) :: apertures
character(100), intent(in) :: afn ! crystal filename
integer(8), dimension(:,:), allocatable, intent(in) :: face_ids

integer(8), parameter :: max_line_length = 150 ! max number of characters in a line of thecrystal file (might need increasing if faces have many vertices)
character(max_line_length) line ! a line in a file
integer(8) i, io ! counting variables
integer(8) num_lines
    
! print*,'========== start sr readApertures'

! read number of lines in aperture file
num_lines = 0
! print*,'opening apertures file:',afn
open(unit = 10, file = afn, status = 'old')
do  ! scan through the lines in the crystal file...
    read(10,"(a)",iostat=io)line
    if (io/=0) exit
    num_lines = num_lines + 1
end do
   
!print*,'number of lines in apertures file: ',num_lines

if(num_lines .ne. size(face_ids,1)) then
    print*,'number of lines in aperture file did not match number of faces. exiting...'
    print*,'num faces:',size(face_ids,1)
    print*,'num_lines:',num_lines
    stop
end if

allocate(apertures(num_lines)) ! allocate apertures array (with same size as face_ids)

rewind(10)
do  i = 1,num_lines ! read aperture from file
    read(10,*) apertures(i)
    !print*,'i',i,' aperture: ',apertures(i)
end do

close(10)

! print*,'========== end sr readApertures'

end subroutine

subroutine init(face_ids, isVisible, isVisiblePlusShadows, isWithinBeam, distances, beamIDs, &
                isWithinBeam_ps, distances_ps, beamIDs_ps, isShadow, ampl_in, ampl_in_ps, la, waveno, &
                rperp, rpar, tperp, tpar, vk71, vk72, vk73, vk91, vk92, vk93, rot_ampl, new_in_ampl, &
                new_in_ampl_ps, trans_ampl_ps, trans_ampl, refl_ampl, refl_ampl_ps, ampl_diff, &
                beam_outbeam_tree, beam_outbeam_tree_counter, interactionCounter)

! subroutine init initialises many variables for the beam tracing loop

logical, dimension(:), allocatable, intent(out) :: isVisible, isVisiblePlusShadows, isWithinBeam, isWithinBeam_ps, isShadow
integer(8), dimension(:,:), allocatable, intent(in) :: face_ids
real(8), dimension(:), allocatable, intent(out) :: distances
real(8), dimension(:), allocatable, intent(out) :: distances_ps
integer(8), dimension(:), allocatable, intent(out) :: beamIDs, beamIDs_ps
complex(8), dimension(:,:,:), allocatable, intent(out) :: ampl_in
complex(8), dimension(:,:,:), allocatable, intent(out) :: ampl_in_ps
real(8), intent(out) :: waveno
real(8), intent(in) :: la
complex(8), dimension(:), allocatable, intent(out) :: rperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(out) :: rpar ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(out) :: tperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(out) :: tpar ! Fresnel coefficient
real(8), dimension(:), allocatable, intent(out) :: vk71, vk72, vk73 ! reflected e-perp vector
real(8), dimension(:), allocatable, intent(out) :: vk91, vk92, vk93 ! reflected prop vector
real(8), dimension(:,:,:), allocatable, intent(out) :: rot_ampl ! rotation matrix for beams incident on eahc facet
complex(8), dimension(:,:,:), allocatable, intent(out) :: new_in_ampl
complex(8), dimension(:,:,:), allocatable, intent(out) :: new_in_ampl_ps
complex(8), dimension(:,:,:), allocatable, intent(out) :: trans_ampl_ps ! transmitted amplitude matrices, including shadowed facets
complex(8), dimension(:,:,:), allocatable, intent(out) :: trans_ampl ! transmitted amplitude matrices
complex(8), dimension(:,:,:), allocatable, intent(out) :: refl_ampl ! reflected amplitude matrices
complex(8), dimension(:,:,:), allocatable, intent(out) :: refl_ampl_ps ! reflected amplitude matrices
complex(8), dimension(:,:,:), allocatable, intent(out) :: ampl_diff ! external diffraction amplitude matrices
type(outbeamtype), dimension(:), allocatable, intent(out) :: beam_outbeam_tree ! outgoing beams from the beam tracing
integer(8), intent(out) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
integer(8), intent(out) :: interactionCounter ! counts the current number of interactions

! print*,'========== start sr init'

! allocate arrays
allocate(isVisible(size(face_ids,1)))
allocate(isVisiblePlusShadows(size(face_ids,1)))
allocate(isWithinBeam(size(face_ids,1)))
allocate(isWithinBeam_ps(size(face_ids,1)))
allocate(distances(size(face_ids,1)))
allocate(distances_ps(size(face_ids,1)))
allocate(beamIDs(size(face_ids,1)))
allocate(beamIDs_ps(size(face_ids,1)))
allocate(isShadow(size(face_ids,1)))
allocate(ampl_in(1:2,1:2,1:size(face_ids,1)))
allocate(ampl_in_ps(1:2,1:2,1:size(face_ids,1)))
allocate(rperp(size(face_ids,1)))
allocate(rpar(size(face_ids,1)))
allocate(tperp(size(face_ids,1)))
allocate(tpar(size(face_ids,1)))
allocate(vk71(size(face_ids,1)))
allocate(vk72(size(face_ids,1)))
allocate(vk73(size(face_ids,1)))
allocate(vk91(size(face_ids,1)))
allocate(vk92(size(face_ids,1)))
allocate(vk93(size(face_ids,1)))
allocate(rot_ampl(1:2,1:2,1:size(face_ids,1)))
allocate(new_in_ampl(1:2,1:2,1:size(face_ids,1)))
allocate(new_in_ampl_ps(1:2,1:2,1:size(face_ids,1)))
allocate(trans_ampl_ps(1:2,1:2,1:size(face_ids,1)))
allocate(trans_ampl(1:2,1:2,1:size(face_ids,1)))
allocate(refl_ampl(1:2,1:2,1:size(face_ids,1)))
allocate(refl_ampl_ps(1:2,1:2,1:size(face_ids,1)))
allocate(beam_outbeam_tree(1:500000)) ! set to 100000 as guess for max outbeams

waveno = 2*pi/la
beam_outbeam_tree_counter = 0 ! counts the current number of beam outbeams
interactionCounter = 0 ! counts the current number of interactions
ampl_in_ps = 0
ampl_in = 0
isVisible = .false.
isVisiblePlusShadows = .false.
isWithinBeam = .false.
isWithinBeam_ps = .false.
distances = 0
distances_ps = 0
beamIDs = 0
beamIDs_ps = 0
isShadow = .false.
ampl_in = 0
ampl_in_ps = 0
rperp = 0
rpar = 0
tperp = 0
tpar = 0
vk71 = 0
vk72 = 0
vk73 = 0
vk91 = 0
vk92 = 0
vk93 = 0
rot_ampl = 0
new_in_ampl = 0
new_in_ampl_ps = 0
trans_ampl_ps = 0
trans_ampl = 0
refl_ampl = 0
refl_ampl_ps = 0

! print*,'========== end sr init'

end subroutine

subroutine makeIncidentBeam(beamV, beamF1, beamN, beamF2, verts, face_ids, beamMidpoints, ampl_beam)

! subroutine makeIncidentBeam makes a simple square incident beam wavefront at a location above the maximum z value of the particle (currently set to 1000)
! the width and length of the wavefront is set larger than the maximum x and y vertex values of the particle (full illumination)

real(8), allocatable, intent(out), dimension(:,:) :: beamV, beamN
integer(8), allocatable, intent(out), dimension(:,:) :: beamF1
integer(8), allocatable, intent(out), dimension(:) :: beamF2
real(8), dimension(:,:) ,allocatable, intent(in) :: verts
real(8), dimension(:,:) ,allocatable, intent(out) :: beamMidpoints
integer(8), dimension(:,:) ,allocatable, intent(in) :: face_ids
complex(8), allocatable, dimension(:,:,:), intent(out) :: ampl_beam ! amplitude matrix of incident beam

real(8) min_x, min_y, max_x, max_y, min_z, max_z, fac

! print*,'========== start sr makeIncidentBeam'

! allocate arrays
! assume square incident beam, 4 vertices, 1 facet, 1 normal
allocate(beamV(1:4,1:3)) ! 4 vertices, xyz coords
allocate(beamF1(1:1,1:4)) ! 1 facet, 4 vertices
allocate(beamN(1:1,1:3)) ! 1 vertex, xyz normal components
allocate(beamF2(1:1)) ! 1 facet

! get the max xyz and min xy coordinates
max_x = maxval(verts(1:size(verts,1),1))
min_x = minval(verts(1:size(verts,1),1))
max_y = maxval(verts(1:size(verts,1),2))
min_y = minval(verts(1:size(verts,1),2))
max_z = maxval(verts(1:size(verts,1),3))
min_z = minval(verts(1:size(verts,1),3))

! print*,'particle max x: ',max_x
! print*,'particle min x: ',min_x
! print*,'particle max y: ',max_y
! print*,'particle min y: ',min_y
! print*,'particle max z: ',max_z
! print*,'particle min z: ',min_z

! use the max/min coords to define the beam
beamF1(1:1,1) = 1
beamF1(1:1,2) = 2
beamF1(1:1,3) = 3
beamF1(1:1,4) = 4
beamN(1:1,1) = 0
beamN(1:1,2) = 0
beamN(1:1,3) = -1
beamF2(1:1) = 1

! full illumination quick implementation
! looping anti-clockwise around the incident beam
fac = 1.1 ! stretch factor
beamV(1,1) = min_x*fac
beamV(1,2) = min_y*fac
beamV(2,1) = min_x*fac
beamV(2,2) = max_y*fac
beamV(3,1) = max_x*fac
beamV(3,2) = max_y*fac
beamV(4,1) = max_x*fac
beamV(4,2) = min_y*fac
! beamV(1:4,3) = max_z*fac
beamV(1:4,3) = 1000 ! changed for comparing with Matlab

!print*,'using truncated initial beam for testing purposes'
!! partial illumination for test implementation
!beamV(1,1) = -3.5
!beamV(1,2) = -10
!beamV(2,1) = -3.5
!beamV(2,2) = 10
!beamV(3,1) = 3.5
!beamV(3,2) = 10
!beamV(4,1) = 3.5
!beamV(4,2) = -10
!beamV(1:4,3) = 50

! make beam midpoints array
allocate(beamMidpoints(1:1,1:3))
beamMidpoints(1,1:3) = sum(beamV(beamF1(1,1:4),1:3),1)/4

allocate(ampl_beam(1:2,1:2,1:1)) ! 4 (2x2) components, 1 facet
ampl_beam = 0
ampl_beam(1,1,1) = 1 ! set diagonal elements to 1
ampl_beam(2,2,1) = 1

! print'(A16,f6.4,A,f6.4,A,f6.4,A,f6.4,A)',' beam ampl in: (',real(ampl_beam(1,1,1)),' + ',imag(ampl_beam(1,1,1)),'i, ',real(ampl_beam(1,2,1)),' + ',imag(ampl_beam(1,2,1)),'i)'
! print'(A16,f6.4,A,f6.4,A,f6.4,A,f6.4,A)','               (',real(ampl_beam(2,1,1)),' + ',imag(ampl_beam(2,1,1)),'i, ',real(ampl_beam(2,2,1)),' + ',imag(ampl_beam(2,2,1)),'i)'

! print'(A,f10.4,f10.4,f10.4)',' beam verts 1: ', beamV(1,1), beamV(1,2), beamV(1,3)
! print'(A,f10.4,f10.4,f10.4)',' beam verts 2: ', beamV(2,1), beamV(2,2), beamV(2,3)
! print'(A,f10.4,f10.4,f10.4)',' beam verts 3: ', beamV(3,1), beamV(3,2), beamV(3,3)
! print'(A,f10.4,f10.4,f10.4)',' beam verts 4: ', beamV(4,1), beamV(4,2), beamV(4,3)

! print'(A,f10.4,f10.4,f10.4)',' beam midpoint: ', beamMidpoints(1,1), beamMidpoints(1,2), beamMidpoints(1,3)

! print*,'========== end sr makeIncidentBeam'

end subroutine

!subroutine midPointsAndAreas(face_ids, verts, Midpoints, faceAreas)
!
!! subroutine midPointsAndAreas computes midpoints and areas of each facet (assumes triangle facets)
!
!integer(8), dimension(:,:) ,allocatable, intent(in) :: face_ids
!real(8), dimension(:,:), allocatable, intent(in) :: verts
!real(8), dimension(:,:), allocatable, intent(out) :: Midpoints
!real(8), dimension(:), allocatable, intent(out) :: faceAreas
!real(8), dimension(1:3) :: vecA, vecB, AcrossB
!real(8) temp_area
!
!integer i,j
!
!print*,'========== start sr midPointsAndAreas'
!
!allocate(Midpoints(size(face_ids,1),3))
!    
!do i = 1, size(face_ids,1) ! for each face
!    Midpoints(i,1:3) = sum(verts(face_ids(i,1:3),1:3),1)/3
!end do
!
!! allocate faceAreas array
!allocate(faceAreas(1:size(face_ids,1)))
!faceAreas = 0 ! initialise
!
!do i = 1, size(face_ids,1) ! for each face
!    do j = 1,3 ! for each vertex in the face
!        vecA = verts(face_ids(i,j),1:3) - Midpoints(i,1:3) ! vector from midpoint to vertex j
!        if (j .eq. 3) then
!            vecB = verts(face_ids(i,1),1:3) - Midpoints(i,1:3) ! vector from midpoint to vertex j+1
!        else
!            vecB = verts(face_ids(i,j+1),1:3) - Midpoints(i,1:3) ! vector from midpoint to vertex j+1
!        end if
!        call cross(vecA,vecB,AcrossB,.false.) ! cross product, no normalisation, calculates parallelepid area
!        temp_area = sqrt(AcrossB(1)**2 + AcrossB(2)**2 + AcrossB(3)**2)/2 ! triangle area is half the above area
!        faceAreas(i) = faceAreas(i) + temp_area ! add to facet area
!    end do
!    !print*,'i',i
!    !print*,'face area',faceAreas(i)
!end do
!
!print*,'========== end sr midPointsAndAreas'
!
!end subroutine

subroutine make_normals(face_ids, verts, norm_ids, norms)
    
! subroutine make_normals recomputes and returns the normals, as well as the corresponding IDs for each face
    
integer(8), dimension(:,:) ,allocatable, intent(in) :: face_ids ! face vertex IDs
integer(8), dimension(:) ,allocatable, intent(out) :: norm_ids ! face vertex IDs
real(8), dimension(:,:) ,allocatable, intent(in) :: verts ! unique vertices
real(8), dimension(:,:) ,allocatable, intent(out) :: norms ! unique vertices
    
real(8), dimension(1:3) :: vec12, vec23, normal ! temporary vectors used to find the facet normal
real(8) temp_verts(1:3,1:3) ! temporary array to hold the xyz components of the first 3 vertices in each facet

integer i
    
! print*,'========== start sr make_normals'
    
! allocate the arrays for the normals and the IDs
!print*,'number of faces: ',size(face_ids,1)
allocate(norm_ids(size(face_ids,1)))
allocate(norms(size(face_ids,1),3))
    
! for each face
do i = 1,size(face_ids,1)
    temp_verts(1:3,1:3) = verts(face_ids(i,1:3),1:3) ! get first 3 vertices of facet
        
    vec12(1:3) = temp_verts(2,1:3) - temp_verts(1,1:3) ! vector from vertex 1 to vertex 2
    vec23(1:3) = temp_verts(3,1:3) - temp_verts(2,1:3) ! vector from vertex 2 to vertex 3
        
    ! cross product to get facet normal
    call cross(vec12,vec23,normal)
        
    !print*,'i',i
    !print*,'normal', normal
        
    ! save to arrays
    norm_ids(i) = i
    norms(i,1:3) = normal(1:3)
        
end do
    
    
! print*,'========== end sr make_normals'
    
end subroutine

subroutine SDATIN(  ifn,                & ! input filename
                    la,                 & ! wavelength
                    rbi,                & ! real refractive index
                    ibi,                & ! imaginary refractive index
                    rec,                & ! number of beam recursions
                    rot_method,         & ! rotation method
                    is_multithreaded,   & ! multithreading
                    num_orients,        & ! number of orientations
                    intellirot,         & ! intelligent multiple rotations method
                    c_method)             ! method of particle file input (read or cc crystal generation)
    
    ! sr SDATIN reads the input file

    integer, intent(out) :: rec ! max number of internal beam recursions
    real(8), intent(out) :: la, rbi, ibi ! wavelength, real and imaginary parts of the refractive index
    character(len=*), intent(in) :: ifn
    character(100), intent(out) :: rot_method ! rotation method
    character(100), intent(out) :: c_method ! method of particle file input
    logical, intent(out) :: is_multithreaded ! whether or not multithreaded subroutines should be used (where appropriate)
    logical, intent(out) :: intellirot ! whether or not to use intelligent euler angle choices for orientation avergaing
    integer, intent(out) :: num_orients ! number of orientations

    ! subroutine PDAS reads inputs from the input file
    ! cfn - crystal filename
    ! la - wavelength
    ! rbi - real part of refractive index
    ! ibi - imaginary part of refractive index
    ! aps - number of aperturesread_real

    is_multithreaded = .false. ! assume no multithreading, unless read from input file
    intellirot = .false. ! assume no intelligent euler angle choices, unless read from input file
    
    print*,'========== start sr SDATIN'
    write(101,*),'======================================================'
    write(101,*),'=====================JOB SETTINGS====================='
    write(101,*),'======================================================'


    ! cfn = read_string(ifn,"cfn") ! get crystal filename
    ! call StripSpaces(cfn) ! remove leading spaces
    ! print*,'particle filename: "',cfn(1:len(trim(cfn))),'"'
    ! write(101,*),'particle filename:      "',cfn(1:len(trim(cfn))),'"'

    c_method = read_string(ifn,"cmethod") ! get crystal filename
    call StripSpaces(c_method) ! remove leading spaces
    print*,'particle input method: "',c_method(1:len(trim(c_method))),'"'
    write(101,*),'particle input method:      "',c_method(1:len(trim(c_method))),'"'

    ! cft = read_string(ifn,"cft") ! get crystal filename
    ! call StripSpaces(cft) ! remove leading spaces
    ! print*,'particle file type: "',cft(1:len(trim(cft))),'"'  
    ! write(101,*),'particle file type:     "',cft(1:len(trim(cft))),'"' 

    ! afn = read_string(ifn,"afn") ! get crystal filename
    ! call StripSpaces(afn) ! remove leading spaces
    ! print*,'apertures filename: "',afn(1:len(trim(afn))),'"'    
    ! write(101,*),'apertures filename:     "',afn(1:len(trim(afn))),'"'    
    
    la = read_real(ifn,"lambda")
    print*,'lambda: ',la
    write(101,*),'lambda:                ',la
    
    rbi = read_real(ifn,"rbi")
    print*,'refractive index real: ',rbi
    write(101,*),'refractive index real:',rbi
    
    ibi = read_real(ifn,"ibi")
    print*,'refractive index imag: ',ibi
    write(101,*),'refractive index imag: ',ibi

    rec = read_int(ifn,"rec")
    print*,'max beam recursions: ',rec
    write(101,*),'max beam recursions:   ',rec
 
    rot_method = readString2(ifn,"rot") ! get rotation method
    if( rot_method(1:len(trim(rot_method))) .eq. 'off' .or. &
    rot_method(1:len(trim(rot_method))) .eq. 'euler' .or. &
    rot_method(1:len(trim(rot_method))) .eq. 'multi' .or. &
    rot_method(1:len(trim(rot_method))) .eq. 'none') then
        ! print*,'rotation method: "',rot_method(1:len(trim(rot_method))),'"'
    else
        print*,'Error, "',rot_method(1:len(trim(rot_method))),'" is not a valid option for rotation method'
        stop
    end if

    is_multithreaded = read_flag(ifn,"mt") ! check input file for multithread option

    intellirot = read_flag(ifn,"intellirot") ! check input file for intelligent euler angle choices for orientation avergaing

    if(rot_method(1:len(trim(rot_method))) .eq. 'multi') then
        num_orients = read_int(ifn,"rot multi")
    else
        num_orients = 1
    end if

    print*,'num_orients: ',num_orients
    write(101,*),'num_orients:           ',num_orients

    print*,'rotation method: "',rot_method(1:len(trim(rot_method))),'"'
    write(101,*),'rotation method:        "',rot_method(1:len(trim(rot_method))),'"'


    if (is_multithreaded) then
        print*,'multithreading: enabled'
        write(101,*),'multithreading:         ','enabled'
    else
        print*,'multithreading: disabled'
        write(101,*),'multithreading:         ','disabled'
    end if

    if (intellirot .and. num_orients .gt. 1) then
        print*,'multirot method: intelligent'
        write(101,*),'multirot method:        intelligent'
    else if (.not. intellirot .and. num_orients .gt. 1) then
        print*,'multirot method: random'
        write(101,*),'multirot method:        random'
    end if

    print*,'========== end sr SDATIN'
    
end subroutine

subroutine PDAL2(ifn, c_method, num_vert, num_face, face_ids, verts, num_face_vert, afn, apertures)

! subroutine PDAL2 reads a file containing the particle geometry
! accepted file types: "obj" - .obj files, "mrt" - macke ray-tracing style

! the inputs are:
! cfn = particle filename
! cft = particle filetype
! the outputs are:
! num_vert      = number of unique particle vertices
! num_face      = number of crystal faces
! verts         = num_vert x 3 array of real values equal to the vertex coordinates, 
!   where each row corresponds to a unique vertex and each column corresponds to the x, y, and z components
! num_face_vert = num_face array of integers equal to the number of vertices in each face, 
!   where each row corresponds to each face
! face_ids      = num_face x num_face_vert_max array of integers equal to the vertex IDs, 
!   where each row corresponds to each face and each column corresponds to the vertices in the face.
!   A vertex ID is used to point to the row of the vertex in verts
        
    character(len=255) cfn ! particle filename
    character(len=255) cft ! particle filetype
    character(100), intent(out) ::  afn ! apertures filename
    integer(8), intent(out) :: num_vert, num_face ! number of unique vertices, number of faces
    ! integer(8), intent(out) :: num_norm ! number of face normals
    integer(8), dimension(:), allocatable, intent(out) :: num_face_vert ! number of vertices in each face
    ! integer(8), dimension(:), allocatable, intent(out) :: norm_ids ! face normal ID of each face
    real(8), dimension(:,:) ,allocatable, intent(out) :: verts ! unique vertices
    ! real(8), dimension(:,:) ,allocatable, intent(out) :: norms ! unique vertices, face vertex IDs, face normals
    integer(8), dimension(:,:) ,allocatable :: face_ids_temp ! temporary array to hold face vertex IDs
    integer(8), dimension(:,:) ,allocatable, intent(out) :: face_ids ! face vertex IDs (for after excess columns have been truncated)
    character(100), intent(in) :: c_method ! method of particle file input
    character(len=*), intent(in) :: ifn
	integer(8), dimension(:), allocatable, intent(out) :: apertures ! taken as parents parent facets

    integer(8), parameter :: num_face_vert_max_in = 20 ! max number of vertices per face
    integer(8), parameter :: max_line_length = 150 ! max number of characters in a line of thecrystal file (might need increasing if faces have many vertices)
    character(max_line_length) line ! a line in a file
    integer(8) face_string_length
    integer(8) entry_count
    logical is_current_char_slash
    integer(8) vertex_count
    integer(8) num_face_vert_max
    logical has_end_of_line_reached
    integer(8) i, io, j, k, l, m, n, p, o, q ! counting variables
    real(8), dimension(:), allocatable :: faceAreas ! area of each facet
    real(8), dimension(:,:), allocatable :: Midpoints ! face midpoints

    print*,'particle input method: "',c_method(1:len(trim(c_method))),'"'
    print*,'========== start sr PDAL2'
    write(101,*),'======================================================'
    write(101,*),'=================PARTICLE INFORMATION================='
    write(101,*),'======================================================'

    if(c_method(1:len(trim(c_method))) .eq. "read") then ! if particle is to be read from file

        cfn = read_string(ifn,"cfn") ! get crystal filename
        call StripSpaces(cfn) ! remove leading spaces
        print*,'particle filename: "',cfn(1:len(trim(cfn))),'"'
        write(101,*),'particle filename:      "',cfn(1:len(trim(cfn))),'"'

        cft = read_string(ifn,"cft") ! get crystal filename
        call StripSpaces(cft) ! remove leading spaces
        print*,'particle file type: "',cft(1:len(trim(cft))),'"'  
        write(101,*),'particle file type:     "',cft(1:len(trim(cft))),'"'

        afn = read_string(ifn,"afn") ! get crystal filename
        call StripSpaces(afn) ! remove leading spaces
        print*,'apertures filename: "',afn(1:len(trim(afn))),'"'    
        write(101,*),'apertures filename:     "',afn(1:len(trim(afn))),'"'

        print*,'crystal file type: "',trim(cft),'"'
        write(101,*),'particle filename:      "',cfn(1:len(trim(cfn))),'"'    
        write(101,*),'particle file type:     "',trim(cft),'"'

        if (trim(cft) .eq. 'obj') then

            face_string_length = 0
            is_current_char_slash = .false.
            entry_count = 0
            has_end_of_line_reached = .true.
            j = 1 ! reset counting variable
            k = 1 ! reset counting variable
            m = 1 ! reset counting variable
            num_vert = 0 ! rest number of unique vertices
            num_face = 0 ! rest number of faces
            ! num_norm = 0 ! rest number of face normals
            num_face_vert = 0 ! initialise
            
            open(unit = 10, file = cfn, status = 'old')
            
            do  ! scan through the lines in the crystal file...
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                !do i = 1, len(line)
                    !if (line(i:i) == "/") line(i:i) = " "   ! replace "/" with " "
                !enddo
                if (line(1:2) .eq. 'v ') then 
                    num_vert = num_vert + 1 ! count the number of unique vertices
                ! else if (line(1:2) .eq. 'vn') then
                !     num_norm = num_norm + 1 ! count the number of unique normals
                else if (line(1:2) .eq. 'f ') then
                    num_face = num_face + 1 ! count the number of faces
                end if
            end do
            
            rewind(10)
            
            allocate(verts(num_vert,3)) ! allocate an array for the crystal vertex coorindates
            allocate(num_face_vert(num_face)) ! allocate array for the number of vertices in each face
            ! allocate(norm_ids(num_face)) ! allocate array for the face normal ID of each face
            allocate(face_ids_temp(num_face,num_face_vert_max_in)) ! allocate an array for the vertex IDs in each face
            ! allocate(norms(num_norm,3)) ! allocate an array for the face normals
            
            do  ! reading through the lines in the crystal file...
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                !do i = 1, len(line)
                    !if (line(i:i) == "/") line(i:i) = " "   ! replace "/" with " "
                !enddo
                if (line(1:2) .eq. 'v ') then 
                    read(line(2:), *) verts(j,1:3)  ! read all items
                    j = j + 1
                else if (line(1:2) .eq. 'vn') then
                    ! read(line(3:), *) norms(m,1:3)  ! read all items
                    m = m + 1
                else if (line(1:2) .eq. 'f ') then
                    ! print*,'last 2 chars of line: "',line(len(line)-1:len(line)),'"'
                    if (line(len(line)-1:len(line)) .ne. '  ') then
                        print*,'Too many vertices in face: ',k,'. Please increase max_line_length'
                        stop
                    end if
                    ! print*,'line(2:)',line(2:40)
                    o = 3 ! counter to move along the line, starting after 'f' prefix
                    vertex_count = 0 ! number of vertices in this face starts at 0
                    ! print*,'length of line',len(line)
                    do o = 3,len(line)
                        if(line(o:o) .eq. ' ') then
                            ! if (entry_count .gt. 0) then
                            !     entry_count = entry_count + 1
                            !     if (vertex_count .eq. 1) then ! use the normal from vertex #1 (can be changed but for now, these are always the same)
                            !         ! print*,'str length',face_string_length
                            !         ! print*,'number in position',entry_count,' has ',face_string_length,' characters'
                            !         ! print*,'face ',k,', normal has index #',line(o-face_string_length:o-1)
                            !         read(line(o-face_string_length:o-1),*) norm_ids(k)
                            !     end if
                            ! end if

                            entry_count = 0
                            face_string_length = 0
                        else if (line(o:o) .eq. '/') then
                            is_current_char_slash = .true.
                            entry_count = entry_count + 1
                            if (entry_count .gt. 0) then
                            ! print*,'number in position',entry_count,' has ',face_string_length,' characters'
                                if (entry_count .eq. 1) then
                                    vertex_count = vertex_count + 1
                                    if (vertex_count .gt. num_face_vert_max_in) then
                                        print*,'face',k,' has',vertex_count,' vertices which is greater than the max value of',num_face_vert_max_in
                                        print*,'Please increase the max vertices per face "num_face_vert_max_in"'
                                        stop
                                    end if
                                    ! print*,'face ',k,', vertex: ',vertex_count,' has index #',line(o-1-face_string_length:o-1)
                                    read(line(o-1-face_string_length:o-1),*) face_ids_temp(k,vertex_count)
                                end if
                            end if
                            face_string_length = 0
                        else
                            face_string_length = face_string_length + 1

                        end if
                    end do
                    num_face_vert(k) = vertex_count
                    ! print*,'face',k,' has',num_face_vert(k),' vertices'
                    ! print*,face_ids_temp(k,1:vertex_count)
                    k = k + 1
                end if
            end do
            
            close(10)

            num_face_vert_max = maxval(num_face_vert) ! get max vertices in a face

            allocate(face_ids(num_face,num_face_vert_max)) ! allocate an array for the vertex IDs in each face
            face_ids = face_ids_temp(1:num_face,1:num_face_vert_max) ! truncate excess columns
            
            print*,'crystal geometry:'
            print*,'number of unique vertices: ',num_vert
            
            print*,'number of unique faces: ',num_face
            
            !print*,'number of unique normals: ',num_norm
            print*,'max vertices per face: ',num_face_vert_max

            
            !! print the number of vertices in each face
            ! print*,'num_face_vert array:'
            ! do i = 1, num_face 
            !     print*,num_face_vert(i)
            ! end do

            !! print the normal IDs of each face
            ! print*,'num_face_vert array:'
            ! do i = 1, num_face
            !     print*,norm_ids(i)
            ! end do

                !! print the vertex IDs of each face
                !print*,'face_ids array:'
                !do i = 1,num_face 
                !    print*,face_ids(i,1:num_face_vert(i))
                !end do  

            ! print*,'verts(face_ids(8,3),1:3) ',verts(face_ids(8,3),1:3) ! test (get the x,y,z components of the 3rd vertex in the 8th face)
            ! print*,'norms(norm_ids(5),1:3) ',norms(norm_ids(5),1:3) ! test (get the x,y,z components of the 5th face)
            
            ! now compute the face midpoints
            
            !allocate(Midpoints(num_face,3))
            !
            !do i = 1, num_face ! for each face
            !    !print*,'vertices for this face: '
            !    !print*,'vertex IDs for this face: '
            !    !print*,face_ids(i,1:num_face_vert(i))
            !    !print*,'vertex coordinates for this face: '
            !    !print*,verts(face_ids(i,1:num_face_vert(i)),1:3)
            !    !print*,'midpoint: '
            !    !print*,sum(verts(face_ids(i,1:num_face_vert(i)),1:3),1)/num_face_vert(i)
            !    Midpoints(i,1:3) = sum(verts(face_ids(i,1:num_face_vert(i)),1:3),1)/num_face_vert(i)
            !end do
            
            !do i = 1, num_face
            !    print*,'Face ',i,' has midpoint: ',Midpoints(i,1:3)
            !end do

        else if (trim(cft) .eq. 'mrt') then ! if macke ray-tracing style geometry file

            ! print*,'opening crystal file'
            open(unit = 10, file = cfn, status = 'old')
            read(10, *) num_face
            print*,'number of unique faces: ',num_face

            allocate(face_ids_temp(num_face,num_face_vert_max_in)) ! allocate an array for the vertex IDs in each face
            allocate(num_face_vert(num_face)) ! allocate array for the number of vertices in each face

            m = 0 ! counter number of vertices assigned to faces
            num_face_vert_max = 0
            do i = 1, num_face  ! scan through the face lines in crystal file
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                
                read(line,*) j
                num_face_vert(i) = j

                ! print*,'face',i,'had',num_face_vert(i),'vertices'

                if (j .gt. num_face_vert_max) num_face_vert_max = j

                if (j .gt. num_face_vert_max_in) then
                    print*,'face',i,' has',j,' vertices which is greater than the max value of',num_face_vert_max_in
                    print*,'Please increase the max vertices per face "num_face_vert_max_in"'
                    stop
                end if

                ! add to face array
                do k = 1, j
                    m = m + 1 ! update vertex counter
                    face_ids_temp(i,j-k+1) = m
                end do


            end do

            ! print*,'finished reading number of vertices in each face'

            num_vert = m

            print*,'total vertices to be read: ',num_vert

            allocate(face_ids(num_face,num_face_vert_max)) ! allocate an array for the vertex IDs in each face
            face_ids(1:num_face,1:num_face_vert_max) = face_ids_temp(1:num_face,1:num_face_vert_max)

            allocate(verts(num_vert,3)) ! allocate an array for the crystal vertex coorindates

            k = 0 ! counter to count how many vertices have been read
            do i = 1, num_face
                do j = 1, num_face_vert(i)
                    read(10,"(a)",iostat=io)line
                    if (io/=0) exit
                    k = k + 1 ! update vertex counter
                    read(line,*) verts(k,1), verts(k,2), verts(k,3)
                end do
            end do

            if (k .ne. num_vert) then
                print*,'error: expected to read', num_vert,'vertices but found',k,'vertices'
                stop
            else
                ! print*,'read expected number of vertices'
            end if

            print*,'succesfully finished reading mrt particle file'
                

            ! do i = 1, num_face ! for each face
            ! do i = 1, 10 ! for each face
            !         !print*,'vertices for this face: '
            !    print*,'vertex IDs for this face: '
            !    print*,face_ids(i,1:num_face_vert_max)
            !    ! print*,'vertex coordinates for this face: '
            !    ! print*,transpose(verts(face_ids(i,1:num_face_vert(i)),1:3))
            !    !print*,'midpoint: '
            !    !print*,sum(verts(face_ids(i,1:num_face_vert(i)),1:3),1)/num_face_vert(i)
            ! !    Midpoints(i,1:3) = sum(verts(face_ids(i,1:num_face_vert(i)),1:3),1)/num_face_vert(i)
            ! end do
            ! stop

            close(10)
        else
            print*,'error: particle geometry file type "',trim(cft),'" is not supported.'
            stop
        end if

        call readApertures(afn, apertures, face_ids) ! read aperture assignments from file
    
    else if(c_method(1:len(trim(c_method))) .eq. "cc_hex") then ! if particle is to be generated according to Chris Collier hex method
        print*,'attempting to make cc crystal'
        call CC_HEX_MAIN(face_ids,verts,apertures)
        num_vert = size(verts,1)
        num_face = size(face_ids,1)
        num_face_vert_max = 3
        print*,'back from cc_hex_main'
    else
        print*,'error: ',c_method(1:len(trim(c_method))),' is not a valid method of particle file input'
        stop
    end if

    write(101,*),'number of vertices: ',num_vert
    write(101,*),'number of faces:    ',num_face
    write(101,*),'max vertices per face:     ',num_face_vert_max

    call midPointsAndAreas(face_ids, verts, Midpoints, faceAreas) ! calculate particle facet areas (for doing some checks in the following sr)

    call area_stats(faceAreas)

    print*,'========== end sr PDAL2'
    ! stop

end subroutine

character(100) function read_string(ifn,var)

! function read_string(ifn,var) reads the input file defined by variable ifn
! the output is a string contained in the variable var

    character(len=*), intent(in) :: ifn ! input filename
    character(len=*), intent(in) :: var ! variable to read
    character(100) output ! output variable
    character(100) line
    integer num_lines, io, i
    logical success

    open(10,file = ifn, status = 'old')

    num_lines = 0  ! initialise line counter
    success = .false.

    do
        read(10,*,iostat=io)
        if (io/=0) exit
        num_lines = num_lines + 1
    end do

    rewind(10) ! back to top of file

    do i = 1,num_lines
        read(10,"(a)",iostat=io)line
        if (line(1:len(var)+1) .eq. var//' ') then
            output = line(len(var)+1:len(line))
            success = .true.
        end if

    end do

    close(10)

    if (.not. success) then
        print*,'Error: "',var,'" input argument not found'
        stop
    end if

    read_string = output

end function

character(100) function read_optional_string(ifn,var)

! function read_optional_string(ifn,var) reads the input file and looks for an optional flag
! if the flag is found, the flag is read and the output is a string contained in the variable var
! if the flag is not found, the output is a string: "#flagnotfound#" with trailing spaces in the variable var

    character(len=*), intent(in) :: ifn ! input filename
    character(len=*), intent(in) :: var ! variable to read
    character(100) output ! output variable
    character(100) line
    integer num_lines, io, i
    logical success

    open(10,file = ifn, status = 'old')

    num_lines = 0  ! initialise line counter
    success = .false.

    do
        read(10,*,iostat=io)
        if (io/=0) exit
        num_lines = num_lines + 1
    end do

    rewind(10) ! back to top of file

    do i = 1,num_lines
        read(10,"(a)",iostat=io)line
        if (line(1:len(var)+1) .eq. var//' ') then
            output = line(len(var)+1:len(line))
            success = .true.
        end if

    end do

    close(10)

    if (.not. success) then
        output = "#flagnotfound#"
    end if

    read_optional_string = output

end function

character(100) function readString2(ifn,var)

! function read_string(ifn,var) reads the input file defined by variable ifn
! the output is a string contained in the variable var
! only reads until the first space (useful for more complicated inputs)

    character(len=*), intent(in) :: ifn ! input filename
    character(len=*), intent(in) :: var ! variable to read
    character(100) output ! output variable
    character(100) line
    integer num_lines, io, i, j
    logical success, success2

    open(10,file = ifn, status = 'old')

    num_lines = 0  ! initialise line counter
    success = .false.
    success2 = .false.

    do
        read(10,*,iostat=io)
        if (io/=0) exit
        num_lines = num_lines + 1
    end do

    rewind(10) ! back to top of file

    do i = 1,num_lines
        read(10,"(a)",iostat=io)line
        if (line(1:len(var)+1) .eq. var//' ') then ! if match found
            j = 0
            do while (success2 .eq. .false.)
                j = j + 1
                ! print*,line(len(var)+1+j:len(var)+1+j)
                if(line(len(var)+1+j:len(var)+1+j) .eq. ' ') success2 = .true.
            end do
            ! print*,'j',j
            output = line(len(var)+2:len(var)+j)
            success = .true.
        end if

    end do

    close(10)

    if (.not. success) then
        print*,'no input for rotation detected. Assuming no rotation required...'
        output = "none"
    end if

    readString2 = output

end function

real(8) function read_real(ifn,var)

! function read_string(ifn,var) reads the input file defined by variable ifn
! the output is a string contained in the variable var

character(len=*), intent(in) :: ifn ! input filename
character(len=*), intent(in) :: var ! variable to read
character(100) string_to_convert
real(8) output ! output variable
character(100) line
integer num_lines, io, i
logical success

open(10,file = ifn, status = 'old')

num_lines = 0  ! initialise line counter
success = .false.

do
    read(10,*,iostat=io)
    if (io/=0) exit
    num_lines = num_lines + 1
end do

rewind(10) ! back to top of file

do i = 1,num_lines
    read(10,"(a)",iostat=io)line
    if (line(1:len(var)+1) .eq. var//' ') then
        string_to_convert = line(len(var)+1:len(line))
        call StripSpaces(string_to_convert)
        read(string_to_convert(1:len(trim(string_to_convert))),*) output
        success = .true.
    end if

end do

close(10)

if (.not. success) then
    print*,'Error getting crystal filename: "',var,'" argument not found'
    stop
end if

read_real = output

end function

integer function read_int(ifn,var)

! function read_string(ifn,var) reads the input file defined by variable ifn
! the output is a string contained in the variable var

character(len=*), intent(in) :: ifn ! input filename
character(len=*), intent(in) :: var ! variable to read
character(100) string_to_convert
integer output ! output variable
character(100) line
integer num_lines, io, i
logical success

open(10,file = ifn, status = 'old')

num_lines = 0  ! initialise line counter
success = .false.

do
    read(10,*,iostat=io)
    if (io/=0) exit
    num_lines = num_lines + 1
end do

rewind(10) ! back to top of file

do i = 1,num_lines
    read(10,"(a)",iostat=io)line
    if (line(1:len(var)+1) .eq. var//' ') then
        string_to_convert = line(len(var)+1:len(line))
        call StripSpaces(string_to_convert)
        read(string_to_convert(1:len(trim(string_to_convert))),*) output
        success = .true.
    end if

end do

close(10)

if (.not. success) then
    print*,'Error getting crystal filename: "',var,'" argument not found'
    stop
end if

read_int = output

end function

logical function read_flag(ifn,var)

! function read_string(ifn,var) reads the input file defined by variable ifn
! the output is a string contained in the variable var

character(len=*), intent(in) :: ifn ! input filename
character(len=*), intent(in) :: var ! variable to read
! character(100) string_to_convert
! integer(8) output ! output variable
character(100) line
integer num_lines, io, i
logical success

open(10,file = ifn, status = 'old')

read_flag = .false. ! assume false
num_lines = 0  ! initialise line counter
success = .false.

do
    read(10,*,iostat=io)
    if (io/=0) exit
    num_lines = num_lines + 1
end do

rewind(10) ! back to top of file

do i = 1,num_lines
    read(10,"(a)",iostat=io)line
    if (line(1:len(var)+1) .eq. var//' ') then
        ! string_to_convert = line(len(var)+1:len(line))
        ! call StripSpaces(string_to_convert)
        ! read(string_to_convert(1:len(trim(string_to_convert))),*) output
        print*,'enabled optional flag: "'//var//'"'
        read_flag = .true.
        success = .true.
    end if

end do

close(10)
! stop
! if (.not. success) then
!     print*,'Error getting crystal filename: "',var,'" argument not found'
!     stop
! end if

! read_int = output

end function

end module input_mod