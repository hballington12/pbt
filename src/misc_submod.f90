! misc_submod.f90
! submodule for misc subroutines and functions used across all modules and the main program
! also contains fundamental constants
    
    module misc_submod
        
        ! use ifport
        use types_mod
        
        implicit none
        
        real(8), parameter, public :: pi = 3.14159265358979 
        
        contains
        
        subroutine output_eulers(alpha_vals,beta_vals,gamma_vals,output_dir,job_params)

            ! writes the euler angles to a file

        real(8), dimension(:), allocatable, intent(in) :: alpha_vals, beta_vals, gamma_vals    
        character(len=255), intent(in) :: output_dir ! output directory
        type(job_parameters_type), intent(in) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details

        integer i

        print*,'outputting eulers to file...'
        open(unit=10,file=trim(output_dir)//"/"//"eulers.dat")
        do i = 1, job_params%num_orients
            write(10,"(f14.8,f14.8,f14.8)") 360D0*alpha_vals(i), 180D0/pi*acos(1-2*beta_vals(i)), 360D0*gamma_vals(i)
        end do
        close(10)  

        end subroutine

        subroutine resume_job(job_params,num_remaining_orients,remaining_orients,mueller_total,mueller_1d_total,output_parameters_total)
            
            ! attempts to pull cached files from cache directory
            
            type(job_parameters_type), intent(inout) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details
            integer, intent(out) :: num_remaining_orients
            real(8), dimension(:,:,:), allocatable, intent(out) :: mueller_total ! mueller matrices
            real(8), dimension(:,:), allocatable, intent(out) :: mueller_1d_total ! mueller matrices
            type(output_parameters_type), intent(out) :: output_parameters_total
            integer, dimension(:), allocatable, intent(out) :: remaining_orients
            
            integer cache_id, i, nlines, io, j
            character(len=32) cache_id_string
            character(len=255) line
            real(8) junk
            
            cache_id = job_params%cache_id
            
            if(cache_id .eq. -1) then
                print*,'error: invalid cache id. did you forget to specify flag "-cache_id <value>"?'
                stop
            else
                print*,'cache id: ',cache_id
            end if
            
            print*,'opening cache files...'
            
            write(cache_id_string,*) cache_id
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/job_params.txt"
            
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/job_params.txt") ! open job_params file
            read(10,*) junk
            ! read(10,*) loop_start
            ! loop_start = loop_start + 1 ! add 1 because loop_start was where the last job part finished
            read(10,*) job_params%cfn
            read(10,*) job_params%cft
            read(10,*) job_params%afn
            read(10,*) job_params%la
            read(10,*) job_params%rbi
            read(10,*) job_params%ibi
            read(10,*) job_params%rec
            read(10,*) job_params%rot_method
            read(10,*) job_params%is_multithreaded
            read(10,*) job_params%num_orients
            read(10,*) job_params%intellirot
            read(10,*) job_params%c_method
            read(10,*) job_params%job_name
            read(10,*) job_params%suppress_2d
            read(10,*) job_params%tri
            job_params%tri = .false. ! dont re-triangulate
            read(10,*) job_params%tri_edge_length
            read(10,*) job_params%tri_roughness
            close(10)
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/output_params.txt"
            
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/output_params.txt") ! open job_params file
            read(10,*) output_parameters_total%abs
            read(10,*) output_parameters_total%scatt
            read(10,*) output_parameters_total%ext
            read(10,*) output_parameters_total%albedo
            read(10,*) output_parameters_total%asymmetry
            read(10,*) output_parameters_total%abs_eff
            read(10,*) output_parameters_total%scatt_eff
            read(10,*) output_parameters_total%ext_eff
            read(10,*) output_parameters_total%geo_cross_sec
            close(10)
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/theta_vals.dat"
            
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/theta_vals.dat") ! open job_params file
            nlines = 0
            do  ! scan through the lines in the crystal file...
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                nlines = nlines + 1
            end do
            ! deallocate if necessary, then reallocate theta vals to match cached file
            if(allocated(job_params%theta_vals)) deallocate(job_params%theta_vals)
            allocate(job_params%theta_vals(1:nlines))
            print*,'num lines in theta file: ',nlines
            rewind(10)
            do i = 1, nlines
                read(10,*) job_params%theta_vals(i)
                ! print*,'job_params%theta_vals(i)',job_params%theta_vals(i)
            end do
            close(10)
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/phi_vals.dat"
            
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/phi_vals.dat") ! open job_params file
            nlines = 0
            do  ! scan through the lines in the crystal file...
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                nlines = nlines + 1
            end do
            ! deallocate if necessary, then reallocate theta vals to match cached file
            if(allocated(job_params%phi_vals)) deallocate(job_params%phi_vals)
            allocate(job_params%phi_vals(1:nlines))
            print*,'num lines in phi file: ',nlines
            rewind(10)
            do i = 1, nlines
                read(10,*) job_params%phi_vals(i)
                ! print*,'job_params%phi_vals(i)',job_params%phi_vals(i)
            end do
            close(10)
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/orient_remaining.dat"
            
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/orient_remaining.dat") ! open job_params file
            nlines = 0
            do  ! scan through the lines in the crystal file...
                read(10,"(a)",iostat=io)line
                if (io/=0) exit
                nlines = nlines + 1
            end do
            ! deallocate if necessary, then reallocate theta vals to match cached file
            num_remaining_orients = nlines
            allocate(remaining_orients(1:num_remaining_orients))
            print*,'num orients remaining: ',num_remaining_orients
            rewind(10)
            do i = 1, num_remaining_orients
                read(10,*) remaining_orients(i)
                print*,'remaining_orients(i)',remaining_orients(i)
            end do
            close(10)
            
            ! read mueller
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/mueller_scatgrid_1d"
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/mueller_scatgrid_1d") ! open job_params file
            
            allocate(mueller_1d_total(1:size(job_params%theta_vals,1),1:16)) ! 1:1 is for each element
            
            do i = 1, size(job_params%theta_vals,1)
                read(10,'(f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
                junk, &
                mueller_1d_total(i,1), mueller_1d_total(i,2), mueller_1d_total(i,3), mueller_1d_total(i,4), &
                mueller_1d_total(i,5), mueller_1d_total(i,6), mueller_1d_total(i,7), mueller_1d_total(i,8), &
                mueller_1d_total(i,9), mueller_1d_total(i,10), mueller_1d_total(i,11), mueller_1d_total(i,12), &
                mueller_1d_total(i,13), mueller_1d_total(i,14), mueller_1d_total(i,15), mueller_1d_total(i,16)
            end do
            
            close(10)
            
            print*,'trying to open: "',"cache/"//trim(adjustl(cache_id_string))//"/mueller_scatgrid"
            open(10,file="cache/"//trim(adjustl(cache_id_string))//"/mueller_scatgrid") ! open job_params file
            
            allocate(mueller_total(1:size(job_params%phi_vals,1),1:size(job_params%theta_vals,1),1:16)) ! 1:1 is for each element
            
            do i = 1, size(job_params%theta_vals,1)
                do j = 1, size(job_params%phi_vals,1)
                    read(10,'(f12.4,f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
                    junk, junk, &
                    mueller_total(j,i,1), mueller_total(j,i,2), mueller_total(j,i,3), mueller_total(j,i,4), &
                    mueller_total(j,i,5), mueller_total(j,i,6), mueller_total(j,i,7), mueller_total(j,i,8), &
                    mueller_total(j,i,9), mueller_total(j,i,10), mueller_total(j,i,11), mueller_total(j,i,12), &
                    mueller_total(j,i,13), mueller_total(j,i,14), mueller_total(j,i,15), mueller_total(j,i,16)                                                                             
                end do
            end do
            
            close(10)
            
            print*,'read cached files. numer of orients remaining: ',num_remaining_orients,"/",job_params%num_orients
            
        end subroutine
        
        subroutine save_params(job_params,loop_index,output_dir,output_parameters_total)
            
            ! saves the job parameters to a file
            
            character(len=255), intent(in) :: output_dir ! cached files directory (if job stops early)
            integer, intent(in) :: loop_index ! current loop index of the main loop
            type(job_parameters_type), intent(in) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details
            type(output_parameters_type), intent(in) :: output_parameters_total
            
            integer i
            
            open(10,file=trim(output_dir)//"/job_params.txt") ! open job_params file
            write(10,*) loop_index
            write(10,*) job_params%cfn
            write(10,*) job_params%cft
            write(10,*) job_params%afn
            write(10,*) job_params%la
            write(10,*) job_params%rbi
            write(10,*) job_params%ibi
            write(10,*) job_params%rec
            write(10,*) job_params%rot_method
            write(10,*) job_params%is_multithreaded
            write(10,*) job_params%num_orients
            write(10,*) job_params%intellirot
            write(10,*) job_params%c_method
            write(10,*) job_params%job_name
            write(10,*) job_params%suppress_2d
            write(10,*) job_params%tri
            write(10,*) job_params%tri_edge_length
            write(10,*) job_params%tri_roughness
            close(10)
            
            open(10,file=trim(output_dir)//"/theta_vals.dat") ! open theta_vals file
            do i = 1, size(job_params%theta_vals,1)
                write(10,*) job_params%theta_vals(i)
            end do
            close(10)
            
            open(10,file=trim(output_dir)//"/phi_vals.dat") ! open phi_vals file
            do i = 1, size(job_params%phi_vals,1)
                write(10,*) job_params%phi_vals(i)
            end do
            close(10)
            
            open(10,file=trim(output_dir)//"/output_params.txt") ! open output_params file
            write(10,*) output_parameters_total%abs
            write(10,*) output_parameters_total%scatt
            write(10,*) output_parameters_total%ext
            write(10,*) output_parameters_total%albedo
            write(10,*) output_parameters_total%asymmetry
            write(10,*) output_parameters_total%abs_eff
            write(10,*) output_parameters_total%scatt_eff
            write(10,*) output_parameters_total%ext_eff
            write(10,*) output_parameters_total%geo_cross_sec
            close(10)
            
        end subroutine
        
        subroutine save_apertures(  apertures,  &
                                    output_dir)
            
            ! saves the apertures to a file
            
            character(len=255), intent(in) :: output_dir ! cached files directory (if job stops early)
            integer, dimension(:), allocatable, intent(in) :: apertures ! apertures asignments for each facet
            
            integer i
            
            open(10,file=trim(output_dir)//"/apertures.dat") ! open pertures file
            
            do i = 1,size(apertures,1)
                write(10,*) apertures(i)
            end do
            
            close(10)
            
        end subroutine
        
        subroutine write_job_params(job_params)
            
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
            character(255) job_name ! name of job / output directory
            type(cc_hex_params_type) cc_hex_params ! parameters for C. Collier Gaussian Random hexagonal columns/plates
            type(job_parameters_type), intent(in) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details
            
            integer i
            
            cfn = job_params%cfn
            cft = job_params%cft
            afn = job_params%afn
            la = job_params%la
            rbi = job_params%rbi
            ibi = job_params%ibi
            rec = job_params%rec
            rot_method = job_params%rot_method
            is_multithreaded = job_params%is_multithreaded
            num_orients = job_params%num_orients
            intellirot = job_params%intellirot
            c_method = job_params%c_method
            job_name = job_params%job_name
            cc_hex_params = job_params%cc_hex_params
            
            write(101,*)'======================================================'
            write(101,*)'=====================JOB SETTINGS====================='
            write(101,*)'======================================================'
            write(101,*)'job_name:                      "',trim(job_name),'"'
            write(101,*)'lambda:                      ',la
            write(101,*)'refractive index real:      ',rbi
            write(101,*)'refractive index imag:       ',ibi
            write(101,*)'max beam recursions:',rec
            write(101,*)'rotation method:               "',rot_method(1:len(trim(rot_method))),'"'
            write(101,*)'no. of orientations:',num_orients
            if (is_multithreaded) then
                write(101,*)'multithreading:                ','enabled'
            else
                write(101,*)'multithreading:                ','disabled'
            end if
            write(101,*)'particle input method:         "',c_method(1:len(trim(c_method))),'"'
            if (intellirot .and. num_orients .gt. 1) then
                write(101,*)'multirot method:        intelligent'
            else if (.not. intellirot .and. num_orients .gt. 1) then
                write(101,*)'multirot method:        random'
            end if
            if(c_method(1:len(trim(c_method))) .eq. "read") then ! if particle is to be read from file
                write(101,*)'particle filename:             "',cfn(1:len(trim(cfn))),'"'
                write(101,*)'particle file type:            "',cft(1:len(trim(cft))),'"'
                write(101,*)'apertures filename:            "',afn(1:len(trim(afn))),'"'
            else if(c_method(1:len(trim(c_method))) .eq. "cc_hex") then ! if particle is to be generated according to Chris Collier hex method
                write(101,*)'gaussian random particle parameters...'
                write(101,*)'L:                                ',cc_hex_params%l
                write(101,*)'hexagon edge length:              ',cc_hex_params%hr
                write(101,*)'no. facets per hex. edge: ',cc_hex_params%nfhr
                write(101,*)'prism edge length:                ',cc_hex_params%pfl
                write(101,*)'no. facets per prism edge:',cc_hex_params%nfpl
                write(101,*)'no. rotations at p-b:     ',cc_hex_params%pher
                write(101,*)'no. rotations at p-p:     ',cc_hex_params%pper
                write(101,*)'no. roughness scales:     ',cc_hex_params%nscales
                do i = 1, cc_hex_params%nscales
                    write(101,'(A,I1,A,f16.8)')' corr. length ',i,':                ',cc_hex_params%cls(i)
                    write(101,'(A,I1,A,f16.8)')' std. deviation ',i,':              ',cc_hex_params%sds(i)
                end do
            end if
            
        end subroutine
        
        subroutine write_outbins(output_dir,theta_vals,phi_vals)
            
            ! writes the far-field angular bins to a file
            
            real(8), dimension(:), allocatable, intent(in) :: theta_vals, phi_vals
            character(len=*), intent(in) :: output_dir
            
            integer i, j
            
            open(10,file=trim(output_dir)//"/"//"outbins")
            
            do i = 1, size(theta_vals,1)
                do j = 1, size(phi_vals,1)
                    write(10,'(f12.4,f12.4)') theta_vals(i)*180.0/pi, phi_vals(j)*180.0/pi
                end do
            end do
            
            close(10)
            
        end subroutine
        
        subroutine make_dir(dir_path_in,cwd_out)
            
            ! makes job directory - quick and dirty method
            
            character(len=*), intent(in) :: dir_path_in
            character(len=255) dir_path
            character(len=255), intent(out) :: cwd_out
            logical exists ! for checking if directories or files exist
            logical result ! true if subdirectory was made, false if subdirectory was not made
            integer job_num ! job number
            character(len=255) job_num_string ! job number
            
            job_num = 1
            result = .false.
            
            ! attempt to make job with specified name
            write(dir_path,*) trim(dir_path_in)
            call StripSpaces(dir_path)
            ! print*,'enquiring at: "',trim(dir_path),'"'
            inquire(file = trim(dir_path), exist = exists)
            ! print*,'exists?',exists
            ! stop
            if(exists .eqv. .false.) then ! if job name is available
                ! print*,'Creating job directory at "',trim(dir_path),'"'
                write(cwd_out,*) trim(dir_path)
                ! result = makedirqq(trim(dir_path))
                call system("mkdir "//trim(dir_path))
            else ! if job name is not available, append numbers until an available name is found
                do while (result .eqv. .false.) ! while a new directory has not been made
                    
                    ! append job_num to directory name
                    write(job_num_string,"(I3)") job_num
                    call StripSpaces(job_num_string)
                    ! print*,'job_num_string:',trim(job_num_string)
                    ! print*,'dir_path: ',trim(dir_path_in)//trim(job_num_string)
                    ! stop
                    write(dir_path,*) trim(dir_path_in)//"_"//trim(job_num_string)
                    call StripSpaces(dir_path)
                    ! print*,'attempting to make directory with name "',trim(dir_path),'"'
                    
                    inquire(file = trim(dir_path), exist = exists)
                    if(exists .eqv. .false.) then 
                        ! print*,'Creating job directory at "',trim(dir_path),'"'
                        write(cwd_out,*) trim(dir_path)
                        ! result = makedirqq(trim(dir_path))
                        call system("mkdir "//trim(dir_path))
                        result = .true.
                    else
                        ! if already exists, add 1 to job number and try again
                        job_num = job_num + 1
                        ! print*,'error: job directory already exists'
                        ! stop
                    end if
                    
                end do
                
            end if
            
            call StripSpaces(cwd_out)
            
            
            
        end subroutine
        
        subroutine make_cache_dir(dir_path_in,cwd_out)
            
            ! makes the a cache directory
            
            character(len=*), intent(in) :: dir_path_in
            character(len=255) dir_path
            character(len=255), intent(out) :: cwd_out
            logical exists ! for checking if directories or files exist
            logical result ! true if subdirectory was made, false if subdirectory was not made
            integer job_num ! job number
            character(len=255) job_num_string ! job number
            
            job_num = 1
            result = .false.
            
            do while (result .eqv. .false.) ! while a new directory has not been made
                
                ! append job_num to directory name
                write(job_num_string,"(I3)") job_num
                call StripSpaces(job_num_string)
                ! print*,'job_num_string:',trim(job_num_string)
                ! print*,'dir_path: ',trim(dir_path_in)//trim(job_num_string)
                ! stop
                write(dir_path,*) trim(dir_path_in)//""//trim(job_num_string)
                call StripSpaces(dir_path)
                ! print*,'attempting to make directory with name "',trim(dir_path),'"'
                
                inquire(file = trim(dir_path), exist = exists)
                if(exists .eqv. .false.) then 
                    ! print*,'Creating job directory at "',trim(dir_path),'"'
                    write(cwd_out,*) trim(dir_path)
                    ! result = makedirqq(trim(dir_path))
                    call system("mkdir "//trim(dir_path))
                    result = .true.
                else
                    ! if already exists, add 1 to job number and try again
                    job_num = job_num + 1
                    ! print*,'error: job directory already exists'
                    ! stop
                end if
                
            end do
            
            call StripSpaces(cwd_out)
            
            
            
        end subroutine
        
        subroutine fix_collinear_vertices(verts, face_ids, num_vert, num_face, num_face_vert, apertures)
            
            ! removes edges which have vertices that lie along them
            ! needs fixing to accommodate multiple collinear vertices lying along 1 edge
            
            real(8), dimension(:,:) ,allocatable, intent(inout) :: verts ! unique vertices
            integer, dimension(:,:) ,allocatable, intent(inout) :: face_ids ! face vertex IDs
            integer, intent(inout) :: num_vert, num_face ! number of unique vertices, number of faces
            integer, dimension(:), allocatable, intent(inout) :: num_face_vert ! number of vertices in each face
            integer, dimension(:), allocatable, intent(inout) :: apertures ! apertures asignments for each facet
            
            integer i, j ,k
            real(8) edge_vector(1:3), a_vector(1:3) ! some vectors
            integer next_index, prev_index
            real(8) edge_vector_length, a_vector_length, a_dot_product, closest_vector_length
            integer num_collinear_vertices ! number of vertices found along an edge
            integer collinear_vertex_id ! id of a vertex found to be collinear
            integer temp_face_ids(1000000,1:10) ! face vertex IDs
            integer temp_num_face_vert(1000000) ! number of vertices in each face
            integer temp_apertures(1000000) ! apertures asignments for each facet
            integer max_verts ! max vertices per face
            logical success
            
            success = .true.
            
            ! save the input face_id and num_face_vert arrays
            do i = 1, num_face
                temp_num_face_vert(i) = num_face_vert(i)
                temp_apertures(i) = apertures(i)
                do j =1, num_face_vert(i)
                    temp_face_ids(i,j) = face_ids(i,j)
                end do
            end do
            
            do i = 1, size(face_ids,1) ! for each face
                ! print*,'i:',i,'num faces:',num_face_vert(i)
                do j = 1, num_face_vert(i) ! for each vertex in the face
                    num_collinear_vertices = 0 ! init counter
                    closest_vector_length = HUGE(1D0) ! get a large number
                    ! get the edge vector
                    if(j .eq. num_face_vert(i)) then
                        next_index = 1
                        prev_index = j - 1
                    else if(j .eq. 1) then
                        next_index = j + 1
                        prev_index = num_face_vert(i)
                    else
                        prev_index = j - 1
                        next_index = j + 1
                    end if
                    edge_vector(1:3) = verts(face_ids(i,next_index),1:3) - verts(face_ids(i,j),1:3) ! get edge vector
                    edge_vector_length = sqrt(edge_vector(1)**2 + edge_vector(2)**2 + edge_vector(3)**2) ! get the length of the edge vector
                    
                    do k = 1, num_vert ! loop over all other vertices
                        if(k .ne. face_ids(i,next_index) .and. k .ne. face_ids(i,j)) then ! ignore vertices part of this edge
                            a_vector(1:3) = verts(k,1:3) - verts(face_ids(i,j),1:3) ! get vector from start of edge to this vertex
                            a_vector_length = sqrt(a_vector(1)**2 + a_vector(2)**2 + a_vector(3)**2) ! get the length of the edge vector
                            a_dot_product = dot_product(edge_vector,a_vector)/(a_vector_length*edge_vector_length) ! get normalised dot product of edge vector and vector to vertex
                            if(a_dot_product .gt. 0.999) then
                                if(a_vector_length .lt. edge_vector_length) then
                                    ! print*,'detected that vertex #',k,'lies on the edge between vertex #',face_ids(i,next_index),'and vertex #',face_ids(i,j)
                                    ! if we pass this point, we have have found a vertex from which we must divide facet
                                    ! we should take care to ensure that this edge does not contain any other vertices (this would require re-calling this subroutine)
                                    num_collinear_vertices = num_collinear_vertices + 1
                                    if(a_vector_length .lt. closest_vector_length) then ! if it was the closest collinear vertex so far
                                        collinear_vertex_id = k ! update a tracker to hold the vertex id
                                        closest_vector_length = a_vector_length ! update the distance
                                    end if
                                end if
                            end if
                        end if
                    end do
                    if(num_collinear_vertices .eq. 0) then
                        ! do nothing
                    else if(num_collinear_vertices .gt. 1) then
                        ! print*,'warning: need to recall collinear fix sr'
                        success = .false.
                    else ! else, if only 1 collinear vertex along this edge
                        ! print*,'first vertex of this edge is #',j,'which is vertex ID:',face_ids(i,j)
                        ! print*,'prev vertex of this edge is #',prev_index,'which is vertex ID:',face_ids(i,prev_index)
                        ! print*,'the collinear vertex along this edge was vertex ID:',collinear_vertex_id
                        ! in the old face, replace the first edge vertex with the collinear one we found
                        temp_face_ids(i,j) = collinear_vertex_id
                        ! also make a new triangular face to fill the gap
                        num_face = num_face + 1
                        temp_apertures(num_face) = apertures(i) ! set as the same aperture assignment
                        temp_num_face_vert(num_face) = 3 ! triangular face, so has 3 vertices
                        temp_face_ids(num_face,3) = face_ids(i,j) ! third vertex is the first vertex of the edge
                        temp_face_ids(num_face,2) = collinear_vertex_id ! second vertex is the collinear vertex we found
                        temp_face_ids(num_face,1) = face_ids(i,prev_index) ! first vertex is the previous vertex from the original facet
                        ! stop
                    end if
                end do
            end do
            
            ! print*,'faces in: ',size(face_ids,1)
            ! print*,'faces after fix: ',num_face
            
            if(success) then
                print*,'sr fix_collinear_vertices: added',num_face - size(face_ids,1),' extra faces to remove some collinear vertices'
                
                max_verts = maxval(num_face_vert)
                
                ! print*,'max number of vertices per face = ',max_verts
                
                ! now reassign
                deallocate(face_ids)
                deallocate(num_face_vert)
                deallocate(apertures)
                ! print*,'deallocated stuff.'
                allocate(face_ids(1:num_face,1:max_verts))
                ! print*,'allocated stuff.'
                allocate(num_face_vert(1:num_face))
                ! print*,'allocated stuff.'
                allocate(apertures(1:num_face))
                ! print*,'allocated stuff.'
                
                do i = 1, num_face
                    num_face_vert(i) = temp_num_face_vert(i)
                    apertures(i) = temp_apertures(i)
                    do j =1, temp_num_face_vert(i)
                        face_ids(i,j) = temp_face_ids(i,j)
                    end do
                end do
            else
                print*,'sr fix_collinear_vertices failed because multiple vertices were found to lie along a single edge'
                print*,'to fix this, please refine triangle settings'
                ! print*,'results should be ok, expect some small energy losses.'
            end if
            
            ! stop
            
            
            
        end subroutine
        
        subroutine merge_vertices(verts, face_ids, num_vert, distance_threshold)
            
            ! merges vertices that are close enough together
            ! general outline:
            ! for each vertex, calculate the distance to all other vertices
            ! if the distance is less than the specified distance threshold, merge them at the centre
            ! else, keep the vertex as it is
            ! add the vertex to a new vertex array and appropriately assign the face_ids
            
            real(8), dimension(:,:) ,allocatable, intent(inout) :: verts ! unique vertices
            integer, dimension(:,:) ,allocatable, intent(inout) :: face_ids ! face vertex IDs (for after excess columns have been truncated)
            integer, intent(inout) :: num_vert ! number of unique vertices, number of faces
            real(8), intent(in) :: distance_threshold ! max distance between two vertices for them to be merged
            
            logical, dimension(:), allocatable :: has_vertex_been_checked ! whether this vertex has been checked for merge
            real(8) this_vertex(1:3), another_vertex(1:3), midpoint(1:3)
            real(8) distance
            integer i, j, k, num_vertices_to_merge
            integer vertices_to_merge_with(1:100)
            real(8) merged_verts(1000000,1:3)
            integer merged_vert_counter
            integer, dimension(:,:) ,allocatable :: merged_face_ids ! face vertex IDs (for after excess columns have been truncated)
            integer, dimension(:) ,allocatable :: vertex_mapping ! maps the old vertices to the merged vertices
            
            print*,'start vertex merge...'
            
            allocate(vertex_mapping(1:size(verts,1)))
            allocate(merged_face_ids(1:size(face_ids,1),1:size(face_ids,2)))
            allocate(has_vertex_been_checked(1:size(verts,1))) ! allocate as same size as number of vertices
            has_vertex_been_checked = .false.
            merged_vert_counter = 0
            
            do i = 1, num_vert ! for each vertex
                if(.not. has_vertex_been_checked(i)) then ! if we havent already checked this vertex
                    num_vertices_to_merge = 0 ! set the counter
                    this_vertex(1:3) = verts(i,1:3) ! get the components of this vertex
                    do j = 1, num_vert ! check against all other vertices
                        if(j .eq. i) then ! ignore self
                            ! do nothing
                        else
                            another_vertex(1:3) = verts(j,1:3) ! get the components of another vertex
                            distance = sqrt((this_vertex(1)-another_vertex(1))**2 + &
                            (this_vertex(2)-another_vertex(2))**2 + &
                            (this_vertex(3)-another_vertex(3))**2) ! distance between the 2 vertices
                            if(distance .lt. distance_threshold) then ! if the 2 vertices are close enough...
                                num_vertices_to_merge = num_vertices_to_merge + 1 ! update counter
                                vertices_to_merge_with(num_vertices_to_merge) = j ! record the vertex id
                            end if
                        end if
                    end do
                    ! print*,'vertex:',i,' - number of vertices to merge with: ',num_vertices_to_merge
                    merged_vert_counter = merged_vert_counter + 1 ! update counter
                    
                    ! now get the midpoint of the vertices we wish to merge
                    midpoint(1:3) = this_vertex(1:3) ! start by setting the midpoint to this vertex
                    vertex_mapping(i) = merged_vert_counter ! record the mapping from old to merged vertices
                    ! print*,'i',i,'vertex_mapping(i)',vertex_mapping(i)
                    do k = 1, num_vertices_to_merge
                        ! print*,vertices_to_merge_with(k)
                        midpoint(1:3) = midpoint(1:3) + verts(vertices_to_merge_with(k),1:3) ! add the other vertices
                        vertex_mapping(vertices_to_merge_with(k)) = merged_vert_counter ! record the mapping from old to merged vertices
                        has_vertex_been_checked(vertices_to_merge_with(k)) = .true.
                    end do
                    midpoint(1:3) = midpoint(1:3) / (num_vertices_to_merge + 1) ! divide by number of vertices to get the midpoint
                    
                    merged_verts(merged_vert_counter,1:3) = midpoint(1:3) ! save the new vertex into a new array
                    
                    has_vertex_been_checked(i) = .true. ! record that this vertex has been checked
                end if
            end do
            
            ! do i = 1, num_vert
            !     print*,'vertex',i,' was mapped to merged vertex',vertex_mapping(i)
            ! end do
            
            ! make the new face_ids array
            do i = 1, size(face_ids,1)
                do j = 1, size(face_ids,2)
                    if((face_ids(i,j)) .ne. 0) then ! if its not a null
                        merged_face_ids(i,j) = vertex_mapping(face_ids(i,j))
                    end if
                end do
            end do
            
            ! repopulate the input arrays with the merged arrays
            deallocate(verts)
            allocate(verts(1:merged_vert_counter,1:3))
            
            do i = 1, merged_vert_counter
                verts(i,1:3) = merged_verts(i,1:3)
            end do
            do i = 1, size(face_ids,1)
                do j = 1, size(face_ids,2)
                    face_ids(i,j) = merged_face_ids(i,j)
                end do
            end do
            
            print*,'merged',num_vert,' vertices into ',merged_vert_counter,' vertices.'
            
            num_vert = merged_vert_counter
            
            ! stop
            
        end subroutine
        
        subroutine triangulate(verts,face_ids,num_vert,num_face,num_face_vert,max_edge_length,flags,apertures,roughness,rank,output_dir)
            
            ! calls triangle to triangulate and subdivide a surface
            ! returns the subdivided surface
            
            real(8), dimension(:,:) ,allocatable, intent(inout) :: verts ! unique vertices
            integer, dimension(:,:) ,allocatable, intent(inout) :: face_ids ! face vertex IDs (for after excess columns have been truncated)
            integer, intent(inout) :: num_vert, num_face ! number of unique vertices, number of faces
            integer, dimension(:), allocatable, intent(inout) :: num_face_vert
            character(len=*) flags
            real(8), intent(in) :: max_edge_length
            integer, dimension(:), allocatable, intent(inout) :: apertures ! apertures asignments for each facet
            real(8), intent(in) :: roughness
            character(len=255), intent(in) :: output_dir ! output directory
            
            integer i, j, junk, vert_counter, vert_counter_total, face_counter, face_counter_total
            real(8), dimension(:,:), allocatable :: v, v0, v1 ! vertices of facet (after lr flip)
            ! integer, dimension(:), allocatable :: my_array ! points some numbers to some other numbers
            real(8) com(1:3) ! facet centre of mass
            real(8) rot(1:3,1:3) ! rotation matrix
            ! integer, dimension(:), allocatable :: num_face_vert_out
            integer, dimension(:), allocatable :: apertures_temp ! apertures asignments for each facet
            
            integer num_verts
            integer num_nodes, num_faces
            real(8) verts_out(1:1000000,1:3) ! vertices after triangulation
            integer face_ids_out(1:1000000,1:3) ! face_ids after triangulation
            integer apertures_out(1:1000000) ! apertures after triangulation
            
            ! real(8), dimension(:,:), allocatable :: verts_out_final
            ! integer, dimension(:,:), allocatable :: face_ids_out_final
            character(255) triangle_cmd ! command for calling triangle
            character(100) a_string
            ! logical, dimension(:), allocatable :: is_boundary ! whether or not a vertex lies on the boundary of an aperture
            integer is_vertex_on_boundary
            logical exists ! for checking if directories or files exist
            integer rank ! for avoiding clashes between processes of different rank
            character(len=16) rank_string
            
            allocate(apertures_temp(1:size(apertures,1))) ! make an array to hold the input apertures
            apertures_temp = apertures ! save the input apertures
            
            vert_counter_total = 0
            face_counter_total = 0
            
            do i = 1, num_face ! for each face
                ! do i = 1, 1 ! for each face
                
                vert_counter = 0
                face_counter = 0
                
                print*,'face:',i
                num_verts = num_face_vert(i)
                print*,'num verts in this face: ',num_verts
                
                ! rotate into x-y plane
                if (allocated(v)) deallocate(v)
                if (allocated(v0)) deallocate(v0)
                if (allocated(v1)) deallocate(v1)
                allocate(v(1:3,1:num_verts))
                allocate(v0(1:3,1:num_verts))
                allocate(v1(1:3,1:num_verts))
                
                do j = 1,num_verts
                    v(1,j) = verts(face_ids(i,j),1)
                    v(2,j) = verts(face_ids(i,j),2)
                    v(3,j) = verts(face_ids(i,j),3)
                    ! print*,v(1:3,j)
                end do
                
                com = sum(v,2)/num_verts ! get centre of mass
                
                ! print*,'com:',com
                
                v0(1,1:num_verts) = v(1,1:num_verts) - com(1) ! translate aperture to centre of mass system
                v0(2,1:num_verts) = v(2,1:num_verts) - com(2)
                v0(3,1:num_verts) = v(3,1:num_verts) - com(3)
                
                call get_rotation_matrix3(v0,rot)
                
                v1 = matmul(rot,v0) ! rotate vertices
                
                ! print*,'rot:',rot
                
                print*,'rank: ',rank
                
                write(rank_string,*) rank
                
                
                call write_face_poly(num_verts,v1,face_ids,i,rank_string,output_dir)
                
                ! stop
                
                ! print*,'start triangulate...'
                ! print*,'===================================='
                ! print*,'===================================='
                
                ! max_edge_length = 0.25
                
                write(a_string,'(f12.6)') max_edge_length**2*sqrt(3D0)/4D0 ! based on area -> length of equilateral triangle
                
                ! print*,'extra flags: ',trim(adjustl(flags))
                
                ! print*,trim(a_string)
                
                ! write(triangle_cmd,*) './src/tri/triangle' // ' -p face.poly -a'//trim(adjustl(a_string))//" "//trim(adjustl(flags))
                write(triangle_cmd,*) './src/tri/triangle' // ' -p ',trim(output_dir)//'/tmp/face'//trim(adjustl(rank_string))//'.poly'//' -a'//trim(adjustl(a_string))//" "//trim(adjustl(flags))
                ! write(triangle_cmd,*) './triangle' // ' -p face.poly -q -Q -B -a'//'0.005'//" "//trim(adjustl(flags))
                
                print*,'triangle command: "',trim(adjustl(triangle_cmd)),'"'
                
                ! check that triangle has been compiled by the user...
                inquire(file="./src/tri/triangle", exist=exists)
                
                if (exists) then
                    print*,'triangle executable found'
                else
                    print*,'error: triangle exectuable not found. please compile triangle at ./src/tri or disable triangulation'
                    print*,'reminder: triangulation is forced automatically if the input particle has facets with more than 3 vertices.'
                    stop
                end if
                
                ! call execute_command_line('./triangle' // ' -p face.poly -q -Q -B -a0.25',wait=.true.) ! remove -Q to reenable triangle print commands
                call execute_command_line(trim(adjustl(triangle_cmd)),wait=.true.) ! remove -Q to reenable triangle print commands
                
                ! print*,'===================================='
                ! print*,'===================================='
                ! print*,'end triangulate'
                
                ! print*,'reading triangulated face back in...'
                
                open(unit=10,file=trim(output_dir)//'/tmp/face'//trim(adjustl(rank_string))//'.1.node',status="old") ! open the triangulated node file
                read(10,*) num_nodes
                print*,'number of vertices after triangulation: ',num_nodes
                do j = 1, num_nodes ! read in each node
                    vert_counter = vert_counter + 1 ! update new vertex counter
                    read(10,*) junk, verts_out(vert_counter + vert_counter_total,1), verts_out(vert_counter + vert_counter_total,2), junk, is_vertex_on_boundary
                    ! print*,'is_vertex_on_boundary:',is_vertex_on_boundary
                    verts_out(vert_counter + vert_counter_total,3) = 0 ! always 0 because triangle works in 2D
                    if(is_vertex_on_boundary .eq. 0) then
                        verts_out(vert_counter + vert_counter_total,3) = verts_out(vert_counter + vert_counter_total,3) + roughness*(rand()-0.5) ! apply roughness to z-component
                    end if
                    ! rotate the vertices back
                    verts_out(vert_counter + vert_counter_total,1:3) = matmul(transpose(rot),verts_out(vert_counter + vert_counter_total,1:3))
                    ! add com back
                    verts_out(vert_counter + vert_counter_total,1) = verts_out(vert_counter + vert_counter_total,1) + com(1)
                    verts_out(vert_counter + vert_counter_total,2) = verts_out(vert_counter + vert_counter_total,2) + com(2)
                    verts_out(vert_counter + vert_counter_total,3) = verts_out(vert_counter + vert_counter_total,3) + com(3)
                    ! print*,'applying roughness'
                    
                    ! print*,verts_out(vert_counter + vert_counter_total,1), verts_out(vert_counter + vert_counter_total,2), verts_out(vert_counter + vert_counter_total,3)
                end do
                close(10)
                
                open(unit=10,file=trim(output_dir)//'/tmp/face'//trim(adjustl(rank_string))//'.1.ele',status="old") ! open the triangulated ele file
                read(10,*) num_faces
                print*,'number of faces after triangulation: ',num_faces
                do j = 1, num_faces ! read in each face
                    face_counter = face_counter + 1 ! update new vertex counter
                    read(10,*) junk, face_ids_out(face_counter + face_counter_total,1), face_ids_out(face_counter + face_counter_total,2), face_ids_out(face_counter + face_counter_total,3)
                    face_ids_out(face_counter + face_counter_total,1) = face_ids_out(face_counter + face_counter_total,1) + vert_counter_total
                    face_ids_out(face_counter + face_counter_total,2) = face_ids_out(face_counter + face_counter_total,2) + vert_counter_total
                    face_ids_out(face_counter + face_counter_total,3) = face_ids_out(face_counter + face_counter_total,3) + vert_counter_total
                    
                    apertures_out(face_counter + face_counter_total) = apertures(i) ! this face aperture set as the same as the input face that it was part of
                    ! print*,face_ids_out(face_counter + face_counter_total,1), face_ids_out(face_counter + face_counter_total,2), face_ids_out(face_counter + face_counter_total,3)
                end do
                close(10)
                
                ! save the total vertices and faces read after this loop iteration
                vert_counter_total = vert_counter_total + vert_counter
                face_counter_total = face_counter_total + face_counter
                
                ! print*,'total vertices read so far: ',vert_counter_total
                ! print*,'total faces read so far: ',face_counter_total
                ! print*,'read triangulated face back in'
                print*,'******************************'
                
            end do
            
            print*,'finished triangulation'
            print*,'total vertices: ',vert_counter_total
            print*,'total faces: ',face_counter_total
            
            num_vert = vert_counter_total
            num_face = face_counter_total
            
            deallocate(verts)
            deallocate(face_ids)
            deallocate(num_face_vert)
            deallocate(apertures)
            
            allocate(num_face_vert(1:face_counter_total))
            allocate(verts(1:vert_counter_total,1:3))
            allocate(face_ids(1:face_counter_total,1:3))
            allocate(apertures(1:face_counter_total))
            
            num_face_vert = 3 ! triangulated
            
            do i = 1, vert_counter_total
                verts(i,1:3) = verts_out(i,1:3)
            end do
            
            do i = 1, face_counter_total
                face_ids(i,1:3) = face_ids_out(i,1:3)
                apertures(i) = apertures_out(i)
            end do
            
            
            
            
            
            
        end subroutine
        
        subroutine get_rotation_matrix3(v0,rot)
            
            ! computes a rotation matrix for rotating into an aperture located in the xy plane
            
            real(8), intent(in) :: v0(1:3,1:3) ! vertices of facet
            real(8), intent(out) :: rot(1:3,1:3) ! rotation matrix
            
            real(8) a1(1:3) ! vector from com to vertex 1
            real(8) b1(1:3) ! vector from com to vertex 2
            real(8) theta1
            real(8) rot1(1:3,1:3) ! rotation matrix (about z-axis)
            real(8) a2(1:3) ! a1 after first rotation
            real(8) b2(1:3) ! b1 after first rotation
            real(8) theta2
            real(8) rot2(1:3,1:3) ! rotation matrix (about x-axis)
            real(8) a3(1:3) ! a1 after second rotation
            real(8) b3(1:3) ! b1 after second rotation
            real(8) theta3
            real(8) rot3(1:3,1:3) ! rotation matrix (about y-axis)
            real(8) a4(1:3) ! a1 after third rotation
            real(8) b4(1:3) ! b1 after third rotation
            real(8) rot4(1:3,1:3) ! flip rotation matrix
            
            a1(1:3) = v0(1:3,1)
            b1(1:3) = v0(1:3,2)
            
            ! create the first rotation matrix (about z-axis)
            theta1 = atan(a1(1)/a1(2))
            if(isnan(theta1)) theta1 = 0 ! no rotation if NaN
            rot1(1,1) = cos(theta1)
            rot1(1,2) = -sin(theta1)
            rot1(1,3) = 0
            rot1(2,1) = sin(theta1)
            rot1(2,2) = cos(theta1)
            rot1(2,3) = 0
            rot1(3,1) = 0
            rot1(3,2) = 0
            rot1(3,3) = 1
            
            ! do the first rotation (about z-axis)
            a2 = matmul(rot1,a1)
            b2 = matmul(rot1,b1)
            
            ! create the second rotation matrix (about x-axis)
            theta2 = -atan(a2(3)/a2(2))
            rot2(1,1) = 1
            rot2(1,2) = 0
            rot2(1,3) = 0
            rot2(2,1) = 0
            rot2(2,2) = cos(theta2)
            rot2(2,3) = -sin(theta2)
            rot2(3,1) = 0
            rot2(3,2) = sin(theta2)
            rot2(3,3) = cos(theta2)
            
            ! do the second rotation (about x-axis)
            a3 = matmul(rot2,a2)
            b3 = matmul(rot2,b2)
            
            ! create the third rotation matrix (about y-axis)
            theta3 = atan(b3(3)/b3(1))
            rot3(1,1) = cos(theta3)
            rot3(1,2) = 0
            rot3(1,3) = sin(theta3)
            rot3(2,1) = 0
            rot3(2,2) = 1
            rot3(2,3) = 0
            rot3(3,1) = -sin(theta3)
            rot3(3,2) = 0
            rot3(3,3) = cos(theta3)
            
            ! print*,'isnan(rot1)',isnan(rot1)
            ! print*,'isnan(rot2)',isnan(rot2)
            ! print*,'isnan(rot3)',isnan(rot3)
            
            ! do the third rotation (about y-axis)
            a4 = matmul(rot3,a3)
            b4 = matmul(rot3,b3)
            
            if (a4(1)*b4(2)-a4(2)*b4(1) .lt. 0) then ! rotate about y-axis (any axis should work) by 180 degrees if facing up
                rot4(1,1) = -1
                rot4(1,2) = 0
                rot4(1,3) = 0
                rot4(2,1) = 0
                rot4(2,2) = 1
                rot4(2,3) = 0
                rot4(3,1) = 0
                rot4(3,2) = 0
                rot4(3,3) = -1
                rot = matmul(rot4,matmul(rot3,matmul(rot2,rot1)))
            else
                rot = matmul(rot3,matmul(rot2,rot1))
            end if
            
        end subroutine
        
        subroutine write_face_poly(num_verts,v1,face_ids,i,rank_string,output_dir)
            
            ! writes a face to a temporary file, ready for triangulation
            
            integer, dimension(:), allocatable :: my_array ! points some numbers to some other numbers
            integer, intent(in) :: num_verts
            real(8), dimension(:,:), allocatable, intent(in) :: v1
            integer, dimension(:,:) ,allocatable, intent(in) :: face_ids ! face vertex IDs (for after excess columns have been truncated)
            character(len=16), intent(in) :: rank_string
            character(len=255), intent(in) :: output_dir ! output directory
            
            integer j, i
            
            ! print*,'my rank string: ',trim(adjustl(rank_string))
            ! stop
            ! print*,'writing face to file'
            allocate(my_array(1:num_verts))
            ! open(unit=10,file="face.poly"//trim(adjustl(rank_string)))
            open(unit=10,file=trim(output_dir)//"/tmp/face"//trim(adjustl(rank_string))//".poly")
            write(10,'(I3,I3,I3,I3)') num_verts, 2, 1, 0 ! vertex header
            do j = 1, num_verts ! write vertices of this face
                write(10,'(I3,F16.8,F16.8)') j, v1(1,j), v1(2,j)
                my_array(j) = face_ids(i,j)
                ! print*,'my_array(j)',my_array(j)
                ! print*,'v1:',v1(1,j), v1(2,j)
            end do
            write(10,'(I3,I3)') num_verts, 0 ! segments header (enforced edges)
            do j = 1, num_verts
                if(j .eq. num_verts) then
                    ! print*,'face_ids(i,j)',face_ids(i,j)
                    ! print*,'findloc(my_array,face_ids(i,j))',findloc(my_array,face_ids(i,j))
                    ! write(10,'(I3,I3,I3)') j, face_ids(i,j), face_ids(i,1)
                    write(10,'(I3,I3,I3)') j, findloc(my_array,face_ids(i,j)), findloc(my_array,face_ids(i,1))
                else
                    ! write(10,'(I3,I3,I3)') j, face_ids(i,j), face_ids(i,j+1)
                    write(10,'(I3,I3,I3)') j, findloc(my_array,face_ids(i,j)), findloc(my_array,face_ids(i,j+1))
                end if
            end do
            write(10,'(I3)') 0 ! number of holes
            close(10)
            
            ! print*,'finished writing face to file'
            
        end subroutine
        
        subroutine PDAS(verts, &
            face_ids, &
            output_dir, &
            num_face_vert,  & ! <-  number of verices in each face
            filename)
            
            ! writes rotated particle to file
            ! to do: add Macke-style output
            
            real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
            integer, dimension(:,:), allocatable, intent(in) :: face_ids ! face vertex IDs
            character(len=*), intent(in) :: output_dir
            character(len=*), intent(in) :: filename
            integer, dimension(:), allocatable, intent(in) :: num_face_vert ! number of vertices in each face
            
            integer num_verts, num_faces, i, j, k
            character(100) my_string, my_string2
            
            print*,'========== start sr PDAS'
            
            num_verts = size(verts,1)
            num_faces = size(face_ids,1)
            
            print*,'writing rotated particle to file...'
            
            write(my_string,*) 'v '
            write(my_string2,*) verts(1,1)
            
            ! print*,'my_string: ',my_string
            ! print*,'my_string: ',my_string2
            
            call StripSpaces(my_string2)
            
            ! print*,'my_string222: ',trim(my_string2)
            
            my_string = "v "//trim(my_string2)//" "//trim(my_string2)//" "//trim(my_string2)
            
            ! print*,'my_string: ',my_string
            ! print*,'output dir: "',trim(output_dir),'"'
            
            open(10,file=trim(output_dir)//"/"//trim(filename)//".obj") ! wavefront format
            ! write(10,*) '# rotated particle'
            do i = 1, num_verts
                ! ! make string for line in file
                ! my_string = "v "
                ! write(my_string2,*) verts(i,1)
                ! write(my_string3,*) verts(i,2)
                ! write(my_string4,*) verts(i,3)
                ! call StripSpaces(my_string2)
                ! call StripSpaces(my_string3)
                ! call StripSpaces(my_string4)
                ! my_string = "v "//trim(my_string2)//" "//trim(my_string3)//" "//trim(my_string4)
                ! write(10,*) trim(my_string)
                
                write(10,'(A1,f16.8,f16.8,f16.8)') "v ", verts(i,1), verts(i,2), verts(i,3)
            end do
            do i = 1, num_faces
                my_string = "f "
                call StripSpaces(my_string)
                do j = 1, num_face_vert(i)
                    write(my_string2,*) face_ids(i,j)
                    call StripSpaces(my_string2)
                    ! print*,'my_string2: ',trim(my_string2)
                    my_string = trim(my_string)//" "//trim(my_string2)
                    ! print*,'my_string: ',trim(my_string) 
                end do
                ! print*,'my_string: ',trim(my_string)
                ! stop
                ! write(10,'(A1,I10,I10,I10)') 'f', face_ids(i,1), face_ids(i,2), face_ids(i,3)
                write(10,'(A100)') adjustl(my_string)
            end do
            close(10)
            
            open(10,file=trim(output_dir)//"/"//trim(filename)//".cry") ! macke format (only for triangulated at the moment)
            ! write(10,*) '# rotated particle'
            write(10,'(I8)') num_faces
            do i = 1, num_faces
                ! print*,'face #',i,' had ',num_face_vert(i),' vertices'
                write(10,*) num_face_vert(i)
            end do
            ! print*,'maxval(num_face_vert)',maxval(num_face_vert)
            ! stop
            do i = 1, num_faces
                do j = 1, num_face_vert(i)
                    k = num_face_vert(i) - j + 1
                    ! k = j
                    write(10,'(f16.8,f16.8,f16.8)') verts(face_ids(i,k),1), verts(face_ids(i,k),2), verts(face_ids(i,k),3)
                end do
            end do
            close(10)
            
            print*,'finished writing rotated particle to file'
            
            print*,'========== end sr PDAS'
            
        end subroutine
        
        subroutine flip_vert_lr(verts)
            
            ! flips a real(8) array called verts with dimension(1:3,1:3) along the 2nd dimension
            
            real(8), intent(inout) :: verts(1:3,1:3)
            real(8) vert_flip(1:3,1:3)
            
            vert_flip(1:3,1) = verts(1:3,3) ! flip
            vert_flip(1:3,2) = verts(1:3,2)
            vert_flip(1:3,3) = verts(1:3,1)
            
            verts(1:3,1:3) = vert_flip(1:3,1:3) ! reassign
            
        end subroutine
        
        subroutine meshgrid_real(array_1, array_2, mesh_1, mesh_2)
            
            ! makes a meshgrid from two 1-dimensional input arrays of type real(8)
            ! attempts to copy the Matlab meshgrid() function
            ! each 1-d array contains the values of each axis of a 2-d grid
            ! ie. x = 0, 1, 2, 3, 4, 5, ...
            !     y = 0, 0.1, 0.2, ...
            ! the outputs are 2-d arrays which contain the values at each grid point
            ! ie. mesh_x = 0, 1, 2, 3, 4, 5, ...
            !              0, 1, 2, 3, 4, 5, ...
            !              0, 1, 2, 3, 4, 5, ...
            !     mesh_y = 0, 0, 0, 0, 0, 0, ...
            !              0.1, 0.1, 0.1, 0.1, 0.1, ...
            !              0.2, 0.2, 0.2, 0.2, 0.2, ...
            ! version for real arrays
            
            real(8), dimension(:), allocatable, intent(in) :: array_1, array_2
            real(8), dimension(:,:), allocatable, intent(out) :: mesh_1, mesh_2
            
            integer array_1_dim, array_2_dim
            integer i, j
            
            array_1_dim = size(array_1,1)
            array_2_dim = size(array_2,1)
            
            ! allocate arrays
            allocate(mesh_1(1:array_2_dim,1:array_1_dim))
            allocate(mesh_2(1:array_2_dim,1:array_1_dim))
            
            do i = 1, array_1_dim
                do j = 1, array_2_dim
                    mesh_1(j,i) = array_1(i)
                    mesh_2(j,i) = array_2(j)
                end do
            end do
            
        end subroutine
        
        subroutine test_ref()
            
            ! compares the test output file with the reference output file to check for coding errors
            
            real(8), dimension(:,:), allocatable :: ref_data, test_data
            integer i, j, num_lines, num_rows1, num_cols1, num_rows2, num_cols2, io, error_count
            
            print*,'start test-ref comparison...'
            
            ! read ref data
            open(10,file = "ref.txt", status = 'old')
            num_lines = 0
            do
                read(10,*,iostat=io)
                if (io/=0) exit
                num_lines = num_lines + 1
            end do
            num_rows1 = num_lines
            print*,'num_lines',num_lines
            num_cols1 = 3
            allocate(ref_data(1:num_rows1,1:num_cols1))
            rewind(10)
            do i = 1, num_rows1
                read(10,*)  ref_data(i,1), ref_data(i,2), ref_data(i,3)
                ! ref_data(i,4), &   
                ! ref_data(i,5), ref_data(i,6), ref_data(i,7), ref_data(i,8)
            end do
            close(10)
            
            ! read test data
            open(10,file = "test.txt", status = 'old')
            num_lines = 0
            do
                read(10,*,iostat=io)
                if (io/=0) exit
                num_lines = num_lines + 1
            end do
            num_rows2 = num_lines
            num_cols2 = 3
            allocate(test_data(1:num_rows2,1:num_cols2))
            rewind(10)
            do i = 1, num_rows2
                read(10,*)  test_data(i,1), test_data(i,2), test_data(i,3)
                ! , test_data(i,4), &   
                !             test_data(i,5), test_data(i,6), test_data(i,7), test_data(i,8)
            end do
            close(10)
            
            ! compare data
            
            if(num_rows1 .ne. num_rows2) then
                print*,'error, inconsistent number of rows. stopping...'
                stop
            end if
            
            error_count = 0
            do i = 1, num_rows1
                do j = 1, num_cols1
                    if(abs(ref_data(i,j) - test_data(i,j)) .gt. 0.00001) then
                        ! print*,'fail: j, i',j,i
                        error_count = error_count + 1
                    end if
                    if(error_count .gt. 1000) then
                        print*,'catastrophic number of errors (> 1000). stopping...'
                        stop
                    end if
                end do
            end do
            
            if(error_count .gt. 0) then
                print*,'fail: ',error_count,'/',num_rows1*num_cols1
            else
                print*,'pass'
            end if
            
        end subroutine
        
        subroutine cat_prop(propagation_vectors2, propagation_vectors3)
            
            ! concatenates propagation_vectors2 onto the 3rd dimension of propagation_vectors3
            ! wrote this in a giffy so its a bit messy
            
            real(8), dimension(:,:,:), allocatable, intent(in) :: propagation_vectors2
            real(8), dimension(:,:,:), allocatable, intent(inout) :: propagation_vectors3
            
            real(8), dimension(:,:,:), allocatable :: temp_var
            
            integer dim1, dim2, dim3, i, j, k, m
            
            dim1 = size(propagation_vectors3,1)
            dim2 = size(propagation_vectors3,2)
            dim3 = size(propagation_vectors3,3)
            
            allocate(temp_var(1:dim1,1:dim2,1:dim3))
            temp_var = 0 ! init
            
            do i = 1, dim1
                do j = 1, dim2
                    do k = 1, dim3
                        temp_var(i,j,k) = propagation_vectors3(i,j,k) ! read in propagation_vectors3 and save it
                    end do
                end do
            end do
            
            deallocate(propagation_vectors3)
            allocate(propagation_vectors3(1:dim1,1:dim2,1:dim3+size(propagation_vectors2,3))) ! make new array size
            
            ! write propagation_vectors3 back in to newly-sized array
            do i = 1, dim1
                do j = 1, dim2
                    do k = 1, dim3
                        propagation_vectors3(i,j,k) = temp_var(i,j,k) ! read in propagation_vectors3 and save it
                    end do
                end do
            end do
            
            ! write propagation_vectors2 on the end
            do i = 1, dim1
                do j = 1, dim2
                    m = 0
                    do k = dim3 + 1, dim3+size(propagation_vectors2,3)
                        m = m + 1
                        propagation_vectors3(i,j,k) = propagation_vectors2(i,j,m) ! read in propagation_vectors3 and save it
                    end do
                end do
            end do
            
        end subroutine
        
        subroutine cat_int_var(var1, var2)
            
            ! takes 2 2-d input arrays and stitches them together
            ! extra space is padded with zeros
            ! arrays are concatenated along the 2nd dimension (ie. column to column)
            ! var1 is what is to be concatenated
            ! var2 is what is to be concatenated on to
            ! version for integer arrays
            
            integer, dimension(:,:), allocatable, intent(in) :: var1
            integer, dimension(:,:), allocatable, intent(inout) :: var2
            
            integer, dimension(:,:), allocatable :: var1_temp
            integer, dimension(:,:), allocatable :: var2_temp
            integer var1_size1, var1_size2 ! array dimension sizes
            integer var2_size1, var2_size2
            integer i, j, max_rows, max_cols, icol, irow
            
            var1_size1 = size(var1,1)
            var1_size2 = size(var1,2)
            var2_size1 = size(var2,1)
            var2_size2 = size(var2,2)
            
            ! print*,'var1_size1',var1_size1
            ! print*,'var1_size2',var1_size2
            ! print*,'var2_size1',var2_size1
            ! print*,'var2_size2',var2_size2
            
            ! allocate and save inputs
            allocate(var1_temp(1:var1_size1,1:var1_size2))
            allocate(var2_temp(1:var2_size1,1:var2_size2))
            var1_temp(1:var1_size1,1:var1_size2) = var1(1:var1_size1,1:var1_size2)
            var2_temp(1:var2_size1,1:var2_size2) = var2(1:var2_size1,1:var2_size2)
            
            deallocate(var2) ! dealllocate so we can concatenate the two arrays into this one
            
            
            max_rows = max(var1_size1,var2_size1) ! find out max rows
            max_cols = var1_size2+var2_size2 ! get totall columns
            ! print*,'allocating var2 with ',max_rows,' rows and ',max_cols,' columns'
            allocate(var2(1:max_rows,1:max_cols)) ! allocate var2 to have this many rows and combined number of columns
            var2 = 0 ! init
            
            ! loop over var2_temp, add to var2
            icol = 0 ! init column counter
            do i = 1, var2_size2 ! loop over cols of var2 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var2_size1 ! loop over rows of var2 in
                    irow = irow + 1
                    var2(irow,icol) = var2_temp(j,i)
                end do
            end do
            do i = 1, var1_size2 ! loop over cols of var1 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var1_size1 ! loop over rows of var1 in
                    irow = irow + 1
                    var2(irow,icol) = var1_temp(j,i)
                end do
            end do
            
        end subroutine
        
        subroutine cat_real_var(var1, var2)
            
            ! takes 2 2-d input arrays and stitches them together
            ! extra space is padded with zeros
            ! arrays are concatenated along the 2nd dimension (ie. column to column)
            ! var1 is what is to be concatenated
            ! var2 is what is to be concatenated on to
            ! version for real arrays
            
            real(8), dimension(:,:), allocatable, intent(in) :: var1
            real(8), dimension(:,:), allocatable, intent(inout) :: var2
            
            real(8), dimension(:,:), allocatable :: var1_temp
            real(8), dimension(:,:), allocatable :: var2_temp
            integer var1_size1, var1_size2 ! array dimension sizes
            integer var2_size1, var2_size2
            integer i, j, max_rows, max_cols, icol, irow
            
            var1_size1 = size(var1,1)
            var1_size2 = size(var1,2)
            var2_size1 = size(var2,1)
            var2_size2 = size(var2,2)
            
            ! print*,'var1_size1',var1_size1
            ! print*,'var1_size2',var1_size2
            ! print*,'var2_size1',var2_size1
            ! print*,'var2_size2',var2_size2
            
            ! allocate and save inputs
            allocate(var1_temp(1:var1_size1,1:var1_size2))
            allocate(var2_temp(1:var2_size1,1:var2_size2))
            var1_temp(1:var1_size1,1:var1_size2) = var1(1:var1_size1,1:var1_size2)
            var2_temp(1:var2_size1,1:var2_size2) = var2(1:var2_size1,1:var2_size2)
            
            deallocate(var2) ! dealllocate so we can concatenate the two arrays into this one
            
            
            max_rows = max(var1_size1,var2_size1) ! find out max rows
            max_cols = var1_size2+var2_size2 ! get totall columns
            ! print*,'allocating var2 with ',max_rows,' rows and ',max_cols,' columns'
            allocate(var2(1:max_rows,1:max_cols)) ! allocate var2 to have this many rows and combined number of columns
            var2 = 0 ! init
            
            ! loop over var2_temp, add to var2
            icol = 0 ! init column counter
            do i = 1, var2_size2 ! loop over cols of var2 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var2_size1 ! loop over rows of var2 in
                    irow = irow + 1
                    var2(irow,icol) = var2_temp(j,i)
                end do
            end do
            do i = 1, var1_size2 ! loop over cols of var1 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var1_size1 ! loop over rows of var1 in
                    irow = irow + 1
                    var2(irow,icol) = var1_temp(j,i)
                end do
            end do
            
        end subroutine
        
        subroutine cat_complex_var(var1, var2)
            
            ! takes 2 2-d input arrays and stitches them together
            ! extra space is padded with zeros
            ! arrays are concatenated along the 2nd dimension (ie. column to column)
            ! var1 is what is to be concatenated
            ! var2 is what is to be concatenated on to
            ! version for complex arrays
            
            complex(8), dimension(:,:), allocatable, intent(in) :: var1
            complex(8), dimension(:,:), allocatable, intent(inout) :: var2
            
            complex(8), dimension(:,:), allocatable :: var1_temp
            complex(8), dimension(:,:), allocatable :: var2_temp
            integer var1_size1, var1_size2 ! array dimension sizes
            integer var2_size1, var2_size2
            integer i, j, max_rows, max_cols, icol, irow
            
            var1_size1 = size(var1,1)
            var1_size2 = size(var1,2)
            var2_size1 = size(var2,1)
            var2_size2 = size(var2,2)
            
            ! print*,'var1_size1',var1_size1
            ! print*,'var1_size2',var1_size2
            ! print*,'var2_size1',var2_size1
            ! print*,'var2_size2',var2_size2
            
            ! allocate and save inputs
            allocate(var1_temp(1:var1_size1,1:var1_size2))
            allocate(var2_temp(1:var2_size1,1:var2_size2))
            var1_temp(1:var1_size1,1:var1_size2) = var1(1:var1_size1,1:var1_size2)
            var2_temp(1:var2_size1,1:var2_size2) = var2(1:var2_size1,1:var2_size2)
            
            deallocate(var2) ! dealllocate so we can concatenate the two arrays into this one
            
            
            max_rows = max(var1_size1,var2_size1) ! find out max rows
            max_cols = var1_size2+var2_size2 ! get totall columns
            ! print*,'allocating var2 with ',max_rows,' rows and ',max_cols,' columns'
            allocate(var2(1:max_rows,1:max_cols)) ! allocate var2 to have this many rows and combined number of columns
            var2 = 0 ! init
            
            ! loop over var2_temp, add to var2
            icol = 0 ! init column counter
            do i = 1, var2_size2 ! loop over cols of var2 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var2_size1 ! loop over rows of var2 in
                    irow = irow + 1
                    var2(irow,icol) = var2_temp(j,i)
                end do
            end do
            do i = 1, var1_size2 ! loop over cols of var1 in
                icol = icol + 1 ! update column counter
                irow = 0 ! init row counter
                do j = 1, var1_size1 ! loop over rows of var1 in
                    irow = irow + 1
                    var2(irow,icol) = var1_temp(j,i)
                end do
            end do
            
        end subroutine
        
        subroutine trim_complex_1d_sr1(array, mask)
            
            ! trims an complex(8) 1d array according to a logical mask
            ! mask and array should have the same length
            ! pseudocode:
            ! do i = 1, length of array
            !   if mask(i) is true, keep array(i)
            !   else if mask(i) is false, remove array(i)
            ! end do
            
            logical, dimension(:), allocatable, intent(in) :: mask
            complex(8), dimension(:), allocatable, intent(inout) :: array
            
            complex(8), dimension(:), allocatable :: array_temp
            integer i, counter
            
            allocate(array_temp(1:size(array,1))) ! allocate
            array_temp = array ! save input
            
            counter = 0 ! init
            do i = 1, size(array,1) ! looping through mask
                if(mask(i)) counter = counter + 1
            end do
            
            ! deallocate input array and allocate trimmed array
            deallocate(array)
            allocate(array(1:counter))
            array = 0 ! init
            counter = 0 ! init
            do i = 1, size(array_temp,1) ! looping through mask
                if(mask(i)) then
                    counter = counter + 1
                    array(counter) = array_temp(i)
                end if
            end do
            
        end subroutine
        
        subroutine trim_real_1d_sr1(array, mask)
            
            ! trims an real(8) 1d array according to a logical mask
            ! mask and array should have the same length
            ! pseudocode:
            ! do i = 1, length of array
            !   if mask(i) is true, keep array(i)
            !   else if mask(i) is false, remove array(i)
            ! end do
            
            logical, dimension(:), allocatable, intent(in) :: mask
            real(8), dimension(:), allocatable, intent(inout) :: array
            
            real(8), dimension(:), allocatable :: array_temp
            integer i, counter
            
            allocate(array_temp(1:size(array,1))) ! allocate
            array_temp = array ! save input
            
            counter = 0 ! init
            do i = 1, size(array,1) ! looping through mask
                if(mask(i)) counter = counter + 1
            end do
            
            ! deallocate input array and allocate trimmed array
            deallocate(array)
            allocate(array(1:counter))
            array = 0 ! init
            counter = 0 ! init
            do i = 1, size(array_temp,1) ! looping through mask
                if(mask(i)) then
                    counter = counter + 1
                    array(counter) = array_temp(i)
                end if
            end do
            
        end subroutine
        
        subroutine trim_int_1d_sr1(array, mask)
            
            ! trims an integer 1d array according to a logical mask
            ! mask and array should have the same length
            ! pseudocode:
            ! do i = 1, length of array
            !   if mask(i) is true, keep array(i)
            !   else if mask(i) is false, remove array(i)
            ! end do
            
            logical, dimension(:), allocatable, intent(in) :: mask
            integer, dimension(:), allocatable, intent(inout) :: array
            
            integer, dimension(:), allocatable :: array_temp
            integer i, counter
            
            allocate(array_temp(1:size(array,1))) ! allocate
            array_temp = array ! save input
            
            counter = 0 ! init
            do i = 1, size(array,1) ! looping through mask
                if(mask(i)) counter = counter + 1
            end do
            
            ! deallocate input array and allocate trimmed array
            deallocate(array)
            allocate(array(1:counter))
            array = 0 ! init
            counter = 0 ! init
            do i = 1, size(array_temp,1) ! looping through mask
                if(mask(i)) then
                    counter = counter + 1
                    array(counter) = array_temp(i)
                end if
            end do
            
        end subroutine
        
        subroutine normalise_vec(vec)
            
            ! normalises a real(8) 3-vector
            
            real(8), intent(inout) :: vec(1:3)
            real(8) nf
            
            nf = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
            
            vec = vec/nf
            
        end subroutine
        
        subroutine unique_int(array_in, array_unique)
            
            ! sorts and outputs unique integers from a 1D array
            ! input and output should be of type "integer, dimension(:), allocatable"
            ! inspired by: https://stackoverflow.com/questions/44198212/a-fortran-equivalent-to-unique
            
            integer, dimension(:), allocatable, intent(in) :: array_in
            integer, dimension(:), allocatable :: unique
            integer, dimension(:), allocatable, intent(out) :: array_unique
            integer min_val, max_val, i
            
            ! allocate unique array
            allocate(unique(1:size(array_in,1)))
            ! get array min and max values
            min_val = minval(array_in)-1
            max_val = maxval(array_in)
            i = 0 ! set counter
            do while (min_val < max_val) ! count up through values in array
                i = i + 1
                min_val = minval(array_in, mask = array_in > min_val) ! pick out the next highest unique value
                unique(i) = min_val
            end do
            ! allocate final unique array
            allocate(array_unique(i), source = unique(1:i))
            ! print "(10i5:)", array_unique
            
        end subroutine
        
        complex(8) function exp2cmplx(arg)
        
        ! exp2complex returns the complex number for a given complex exponential
        ! ie. e^{1i*arg} gives an output of cos(arg) + 1i*sin(arg)
        
        real(8), intent(in) :: arg
        
        exp2cmplx = cmplx(cos(arg),sin(arg),8)
        
    end function
    
    subroutine StripSpaces(string) 
        
        ! taken from: https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
        ! removes leading spaces from a character array
        
        character(len=*) :: string
        integer :: stringLen 
        integer :: last, actual
        
        stringLen = len (string)
        last = 1
        actual = 1
        
        do while (actual < stringLen)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ' '
            else
                last = last + 1
                if (actual < last) &
                actual = last
            endif
        end do
        
    end subroutine
    
    subroutine midPointsAndAreas(face_ids, verts, midpoints, face_areas, num_face_vert)
        
        ! subroutine midPointsAndAreas computes midpoints and areas of each facet
        ! assumes triangle facets
        ! to do: extend to quads
        
        integer, dimension(:,:) ,allocatable, intent(in) :: face_ids
        real(8), dimension(:,:), allocatable, intent(in) :: verts
        real(8), dimension(:,:), allocatable, intent(out) :: midpoints
        real(8), dimension(:), allocatable, intent(out) :: face_areas
        real(8), dimension(1:3) :: vec_a, vec_b, a_cross_b
        real(8) temp_area
        integer, dimension(:), allocatable, intent(in) :: num_face_vert ! number of vertices in each face
        
        integer i,j, num_verts
        
        !print*,'========== start sr midPointsAndAreas'
        
        allocate(midpoints(size(face_ids,1),3))
        
        num_verts = size(face_ids,2) ! assumes all faces have the same number of vertices
        
        do i = 1, size(face_ids,1) ! for each face
            num_verts = num_face_vert(i)
            midpoints(i,1:3) = sum(verts(face_ids(i,1:num_verts),1:3),1)/num_verts
        end do
        
        ! allocate face_areas array
        allocate(face_areas(1:size(face_ids,1)))
        face_areas = 0 ! initialise
        
        do i = 1, size(face_ids,1) ! for each face
            num_verts = num_face_vert(i)
            do j = 1,num_verts ! for each vertex in the face
                vec_a = verts(face_ids(i,j),1:3) - midpoints(i,1:3) ! vector from midpoint to vertex j
                if (j .eq. num_verts) then
                    vec_b = verts(face_ids(i,1),1:3) - midpoints(i,1:3) ! vector from midpoint to vertex j+1
                else
                    vec_b = verts(face_ids(i,j+1),1:3) - midpoints(i,1:3) ! vector from midpoint to vertex j+1
                end if
                call cross(vec_a,vec_b,a_cross_b,.false.) ! cross product, no normalisation, calculates parallelepid area
                temp_area = sqrt(a_cross_b(1)**2 + a_cross_b(2)**2 + a_cross_b(3)**2)/2 ! triangle area is half the above area
                face_areas(i) = face_areas(i) + temp_area ! add to facet area
            end do
            !print*,'i',i
            !print*,'face area',face_areas(i)
        end do
        
        !print*,'========== end sr midPointsAndAreas'
        
    end subroutine
    
    subroutine cross(a,b,c,normalise)
        
        ! calculates a cross b and returns c
        ! assumes normalisation, unless optional logical is given stating otherwise
        
        real(8), dimension(1:3), intent(in) :: a, b
        real(8), dimension(1:3), intent(out) :: c
        real(8) nf
        logical, intent(in), optional :: normalise
        
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
        
        nf = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
        
        if(present(normalise)) then
            if(normalise .eqv. .false.) nf = 1
        end if
        
        c = c/nf
        
        
    end subroutine
    
    subroutine simpne ( ntab, x, y, result )
        
        ! a method for integration
        
        !*****************************************************************************80
        !
        !! SIMPNE approximates the integral of unevenly spaced data.
        !
        !  Discussion:
        !
        !    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
        !    to the data and integrates that exactly.
        !
        !  Modified:
        !
        !    10 February 2006
        !
        !  Reference:
        !
        !    Philip Davis, Philip Rabinowitz,
        !    Methods of Numerical Integration,
        !    Second Edition,
        !    Dover, 2007,
        !    ISBN: 0486453391,
        !    LC: QA299.3.D28.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) NTAB, number of data points.  
        !    NTAB must be at least 3.
        !
        !    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
        !    in order.
        !
        !    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
        !
        !    Output, real ( kind = 8 ) RESULT.
        !    RESULT is the approximate value of the integral.
        !
        implicit none
        
        integer ( kind = 4 ) ntab
        
        real ( kind = 8 ) del(3)
        real ( kind = 8 ) e
        real ( kind = 8 ) f
        real ( kind = 8 ) feints
        real ( kind = 8 ) g(3)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) n
        real ( kind = 8 ) pi(3)
        real ( kind = 8 ) result
        real ( kind = 8 ) sum1
        real ( kind = 8 ) x(ntab)
        real ( kind = 8 ) x1
        real ( kind = 8 ) x2
        real ( kind = 8 ) x3
        real ( kind = 8 ) y(ntab)
        
        result = 0.0D+00
        
        if ( ntab <= 2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SIMPNE - Fatal error!'
            write ( *, '(a)' ) '  NTAB <= 2.'
            stop 1
        end if
        
        n = 1
        
        do
            
            x1 = x(n)
            x2 = x(n+1)
            x3 = x(n+2)
            e = x3 * x3- x1 * x1
            f = x3 * x3 * x3 - x1 * x1 * x1
            feints = x3 - x1
            
            del(1) = x3 - x2
            del(2) = x1 - x3
            del(3) = x2 - x1
            
            g(1) = x2 + x3
            g(2) = x1 + x3
            g(3) = x1 + x2
            
            pi(1) = x2 * x3
            pi(2) = x1 * x3
            pi(3) = x1 * x2
            
            sum1 = 0.0D+00
            do i = 1, 3
                sum1 = sum1 + y(n-1+i) * del(i) &
                * ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
            end do
            result = result - sum1 / ( del(1) * del(2) * del(3) )
            
            n = n + 2
            
            if ( ntab <= n + 1 ) then
                exit
            end if
            
        end do
        
        if ( mod ( ntab, 2 ) /= 0 ) then
            return
        end if
        
        n = ntab - 2
        x3 = x(ntab)
        x2 = x(ntab-1)
        x1 = x(ntab-2)
        e = x3 * x3 - x2 * x2
        f = x3 * x3 * x3 - x2 * x2 * x2
        feints = x3 - x2
        
        del(1) = x3 - x2
        del(2) = x1 - x3
        del(3) = x2 - x1
        
        g(1) = x2 + x3
        g(2) = x1 + x3
        g(3) = x1 + x2
        
        pi(1) = x2 * x3
        pi(2) = x1 * x3
        pi(3) = x1 * x2
        
        sum1 = 0.0D+00
        do i = 1, 3
            sum1 = sum1 + y(n-1+i) * del(i) * &
            ( f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints )
        end do
        
        result = result - sum1 / ( del(1) * del(2) * del(3) )
        
        return
        end
        
        SUBROUTINE PROUST(T) !Remembrance of time passed.
            
            ! taken from: https://rosettacode.org/wiki/Convert_seconds_to_compound_duration#Fortran
            ! used to format the predicted remaining time
            
            INTEGER T !The time, in seconds. Positive only, please.
            INTEGER NTYPES !How many types of time?
            PARAMETER (NTYPES = 5) !This should do.
            INTEGER USIZE(NTYPES) !Size of the time unit.
            CHARACTER*3 UNAME(NTYPES)!Name of the time unit.
            PARAMETER (USIZE = (/7*24*60*60, 24*60*60, 60*60,   60,    1/)) !The compiler does some arithmetic.
            PARAMETER (UNAME = (/      "wk ",      "d  ",  "hr ","min","sec"/)) !Approved names, with trailing spaces.
            CHARACTER*28 TEXT !A scratchpad.
            INTEGER I,L,N,S !Assistants.
            S = T !A copy I can mess with.
            L = 0 !No text has been generated.
            DO I = 1,NTYPES !Step through the types to do so.
                N = S/USIZE(I) !Largest first.
                IF (N.GT.0) THEN !Above the waterline?
                    S = S - N*USIZE(I) !Yes! Remove its contribution.
                    IF (L.GT.0) THEN !Is this the first text to be rolled?
                        L = L + 2  !No.
                        TEXT(L - 1:L) = ", " !Cough forth some punctuation.
                    END IF !Now ready for this count.
                    WRITE (TEXT(L + 1:),1) N,UNAME(I) !Place, with the unit name.
                    1       FORMAT (I0,1X,A) !I0 means I only: variable-length, no leading spaces.
                    L = LEN_TRIM(TEXT) !Find the last non-blank resulting.
                END IF !Since I'm not keeping track.
            END DO !On to the next unit.
            !  Cast forth the result.
            print*, ">> ",TEXT(1:L)," <<" !With annotation.
            write(101,*) ">> ",TEXT(1:L)," <<" !With annotation.
            END !Simple enough with integers.
            
            
            subroutine progress_bar(iteration, maximum)

                ! https://github.com/macie/fortran-libs/blob/master/progressbar.f95
                !
                ! Prints progress bar.
                !
                ! Args: 
                !     iteration - iteration number
                !     maximum - total iterations
                !
                    implicit none
                    integer :: iteration, maximum
                    integer :: counter
                    integer :: step, done
                
                    step = nint(iteration * 100 / (1.0 * maximum))
                    done = floor(step / 10.0)  ! mark every 10%
                
                    do counter = 1, 36                    ! clear whole line - 36 chars
                        write(6,'(a)',advance='no') '\b'  ! (\b - backslash)
                    end do
                
                    write(6,'(a)',advance='no') ' -> In progress... ['
                    if (done .LE. 0) then
                        do counter = 1, 10
                            write(6,'(a)',advance='no') '='
                        end do
                    else if ((done .GT. 0) .and. (done .LT. 10)) then
                        do counter = 1, done
                            write(6,'(a)',advance='no') '#'
                        end do
                        do counter = done+1, 10
                            write(6,'(a)',advance='no') '='
                        end do 
                    else
                        do counter = 1, 10
                            write(6,'(a)',advance='no') '#'
                        end do
                    end if
                    write(6,'(a)',advance='no') '] '
                    write(6,'(I3.1)',advance='no') step
                    write(6,'(a)',advance='no') '%'
                end

        end module misc_submod