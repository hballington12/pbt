! diff_mod.f90
! module containing various subroutines used for diffraction :)

module diff_mod

    use misc_submod
    use types_mod
    use omp_lib
    
    implicit none
        
    contains

    subroutine getRotationMatrix(rot, vk71, vk72, vk73, ev11, ev12, ev13, ev31, ev32, ev33)

        ! get rotation matrix for rotation about vector ev3 from plane perp to ev1 to plane perp to vk7
        
        real(8), intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
        real(8), intent(in) :: ev11, ev12, ev13 ! incident e-perp vector
        real(8), intent(in) :: ev31, ev32, ev33 ! incident propagation vector
        real(8), intent(inout) :: rot(1:2,1:2) ! rotation matrix
        
        real(8) myfunc_arg
        real(8) evo21, evo22, evo23
        
        myfunc_arg = vk71*ev11 + vk72*ev12 + vk73*ev13
        
        evo21 = ev32*ev13 - ev33*ev12
        evo22 = ev33*ev11 - ev31*ev13
        evo23 = ev31*ev12 - ev32*ev11
        
        rot(1,1) = myfunc_arg
        rot(1,2) = -(vk71*evo21 + vk72*evo22 + vk73*evo23)
        rot(2,1) = +(vk71*evo21 + vk72*evo22 + vk73*evo23)
        rot(2,2) = myfunc_arg
        
        end subroutine

    subroutine random_rotation(cos_rot,sin_rot)

        ! unused

        real(8), intent(out) :: cos_rot
        real(8), intent(out) :: sin_rot
        real(8) rand, test1, help1, w1, w2

        test1 = 2
        do while (test1 .gt. 1) ! effective search for sin/cos(phi)...
          call random_number(rand)
          w1 = 1. - 2.*rand
          call random_number(rand)
          w2 = 1. - 2.*rand ! ... if phi is random
          test1 = w1*w1 + w2*w2 
        end do
        help1 = sqrt(test1)
        cos_rot = w1/help1
        sin_rot = w2/help1
      end subroutine random_rotation

    subroutine trim_outbeam_tree(beam_tree,beam_tree_counter)

    ! trims the outbeam tree
    ! outgoing beams are ignored if the amplitude matrix has less than 1e-6 energy

    type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_tree ! beam_tree to be trimmed
    integer(8), intent(inout) :: beam_tree_counter ! counts the current number of beam outbeams

    integer i, j
    type(outbeamtype), dimension(:), allocatable :: beam_tree_temp ! beam_tree temporary copy
    real(8) energy
    
    allocate(beam_tree_temp(1:size(beam_tree,1)))

    j = 0
    ! do i = 1, size(beam_tree,1) ! loop over beam_tree entries
    do i = 1, beam_tree_counter ! loop over beam_tree entries
            energy = 0.5*(  beam_tree(i)%ampl(1,1)*conjg(beam_tree(i)%ampl(1,1)) + &
                            beam_tree(i)%ampl(1,2)*conjg(beam_tree(i)%ampl(1,2)) + &
                            beam_tree(i)%ampl(2,1)*conjg(beam_tree(i)%ampl(2,1)) + &
                            beam_tree(i)%ampl(2,2)*conjg(beam_tree(i)%ampl(2,2))) ! calc energy out
        if(energy .gt. 1e-6) then ! if significant energy
            j = j + 1 ! update total number of outbeams
            beam_tree_temp(j) = beam_tree(i) ! add entry to trimmed outbeam tree
        end if
    end do

    print*,'trimmed outbeam tree from',beam_tree_counter,' to ',j,' outbeams'

    beam_tree_counter = j ! replace old beam counter with adjusted one
    beam_tree = beam_tree_temp ! replace in beam_tree with trimmed one

    end subroutine

    subroutine transform_bins2(x,y,z,com,rot,x3,y3,z3)
    
        ! translates and rotates bins to aperture system (loop method)
    
    real(8), dimension(:,:), allocatable, intent(in) :: x, y, z ! far-field bin positions
    real(8), intent(in) :: com(1:3) ! aperture centre of mass
    real(8), intent(in) :: rot(1:3,1:3) ! update to new rotation matrix which aligns with incidence
    real(8), dimension(:,:), allocatable, intent(out) :: x3, y3, z3 ! far-field bin positions
    
    real(8), dimension(:,:), allocatable :: x1, y1, z1 ! far-field bin positions after translation to new com
    
    integer(8) i, j
    
    ! ===== start transform bins =====
    
    ! allocate translated far-field bins
    allocate(x1(1:size(x,1),1:size(x,2)))
    allocate(y1(1:size(x,1),1:size(x,2)))
    allocate(z1(1:size(x,1),1:size(x,2)))
    ! allocate arrays to hold reshaped, rotated arrays
    allocate(x3(1:size(x,1),1:size(x,2)))
    allocate(y3(1:size(x,1),1:size(x,2)))
    allocate(z3(1:size(x,1),1:size(x,2)))
    
    ! translate far-field bins to new com
    x1 = x - com(1)
    y1 = y - com(2)
    z1 = z - com(3)
    
    ! rotate
    do i = 1, size(x,1)
        do j = 1, size(x,2)
            x3(i,j) = rot(1,1)*x1(i,j) + rot(1,2)*y1(i,j) + rot(1,3)*z1(i,j)
            y3(i,j) = rot(2,1)*x1(i,j) + rot(2,2)*y1(i,j) + rot(2,3)*z1(i,j)
            z3(i,j) = rot(3,1)*x1(i,j) + rot(3,2)*y1(i,j) + rot(3,3)*z1(i,j)
        end do
    end do
    
    end subroutine
    
    subroutine transform_bins(x,y,z,com,rot,x3,y3,z3)
    
        ! translates and rotates bins to aperture system (matmul method)
        ! unused
    
    real(8), dimension(:,:), allocatable, intent(in) :: x, y, z ! far-field bin positions
    real(8), intent(in) :: com(1:3) ! aperture centre of mass
    real(8), intent(in) :: rot(1:3,1:3) ! update to new rotation matrix which aligns with incidence
    real(8), dimension(:,:), allocatable, intent(out) :: x3, y3, z3 ! far-field bin positions
    
    real(8), dimension(:,:), allocatable :: x1, y1, z1 ! far-field bin positions after translation to new com
    real(8), dimension(:), allocatable :: x2, y2, z2
    real(8), dimension(:,:), allocatable :: binvecs
    real(8), dimension(:,:), allocatable :: rotbinvecs
    
    integer(8) i
    
    ! ===== start transform bins =====
    
    ! allocate translated far-field bins
    allocate(x1(1:size(x,1),1:size(x,2)))
    allocate(y1(1:size(x,1),1:size(x,2)))
    allocate(z1(1:size(x,1),1:size(x,2)))
    
    ! translate far-field bins to new com
    x1 = x - com(1)
    y1 = y - com(2)
    z1 = z - com(3)
    
    ! allocate arrays to hold reshaped arrays, ready for rotating
    allocate(x2(1:size(x1,1)*size(x1,2)))
    allocate(y2(1:size(x1,1)*size(x1,2)))
    allocate(z2(1:size(x1,1)*size(x1,2)))
    
    ! reshape for rotating
    x2 = reshape(x1,shape(x2))
    y2 = reshape(y1,shape(y2))
    z2 = reshape(z1,shape(z2))
    
    ! allocate bin vectors array
    allocate(binvecs(1:3,1:size(x2)))
    
    ! stitch together for matrix multiplication
    do i = 1, size(x2)
        binvecs(1,i) = x2(i)
        binvecs(2,i) = y2(i)
        binvecs(3,i) = z2(i)
    end do
    
    ! allocate rotated bin vectors array
    allocate(rotbinvecs(1:3,1:size(x2)))
    
    ! rotate
    rotbinvecs = matmul(rot,binvecs)
    
    ! allocate arrays to hold reshaped, rotated arrays
    allocate(x3(1:size(x,1),1:size(x,2)))
    allocate(y3(1:size(x,1),1:size(x,2)))
    allocate(z3(1:size(x,1),1:size(x,2)))
    
    x3 = reshape(rotbinvecs(1,1:size(rotbinvecs,2)),shape(x3))
    y3 = reshape(rotbinvecs(2,1:size(rotbinvecs,2)),shape(x3))
    z3 = reshape(rotbinvecs(3,1:size(rotbinvecs,2)),shape(x3))
    
    end subroutine
    
    
    subroutine karczewski(  diff_ampl, & ! polarisation matrix
                            m, & ! KW perp field vector
                            k, & ! KW par field vector
                            prop2, & ! outgoing propagation direction in aperture system
                            x3, y3, z3) ! x, y, z coordinate of far-field bin
        
        ! sr karczewski computes the diff ampl matrix for polarisation of far-field diffraction at a far-field bin
        ! to do: add a debug mode with the numerical checks found in matlab code
    
    real(8), intent(out) :: diff_ampl(1:2,1:2)
    real(8), intent(out) :: m(1:3)
    real(8), intent(out) :: k(1:3)
    real(8), intent(in) :: prop2(1:3)
    real(8), intent(in) :: x3, y3, z3 ! far-field bin positions
    real(8) bin_vec_size ! far-field bin distances
    real(8) big_kx, big_ky, big_kz ! outward propagation vector in aperture system
    real(8) frac
    real(8) a1m, b2m, a1e, b2e, b1m, a2e, a1em, a2em, b1em, b2em
    
    ! print*,'start karcewski...'

    bin_vec_size = sqrt(x3**2 + y3**2 + z3**2) ! distance to each far-field bin

    k(1) = x3/bin_vec_size ! propagation vector components for each bin vector in aperture system
    k(2) = y3/bin_vec_size
    k(3) = z3/bin_vec_size
    
    big_kx = prop2(1) ! propagation direction in aperture system
    big_ky = prop2(2)
    big_kz = prop2(3)
    
    m(1) = -k(1)*k(2)/sqrt(1-k(2)**2) ! perp field direction
    m(2) = sqrt(1-k(2)**2)
    m(3) = -k(2)*k(3)/sqrt(1-k(2)**2)
    
    frac = sqrt((1-k(2)**2)/(1-big_ky**2)) ! precalculate a factor
    
    ! KW coefficients
    a1m = -big_kz*frac
    b2m = -k(3)/frac
    a1e = b2m
    b2e = a1m
    b1m = -k(1)*k(2)/frac + big_kx*big_ky*frac
    a2e = -b1m
    
    ! combined (e-m theory) KW coefficients
    a1em = 0.5*(a1m+a1e)
    a2em = 0.5*a2e
    b1em = 0.5*b1m 
    b2em = 0.5*(b2m+b2e)

    ! ##### e-m theory #####
    diff_ampl(1,1) = a1em
    diff_ampl(1,2) = b1em
    diff_ampl(2,1) = a2em
    diff_ampl(2,2) = b2em
    ! ######################

    ! ! ###### e theory ######
    ! diff_ampl11 = a1e
    ! diff_ampl12 = 0
    ! diff_ampl21 = a2e
    ! diff_ampl22 = b2e
    ! ! ######################    

    ! ! ###### m theory ######
    ! diff_ampl11 = a1m
    ! diff_ampl12 = b1m
    ! diff_ampl21 = 0
    ! diff_ampl22 = b2m
    ! ! ###################### 

    ! print*,'end karcewski...'
    
    end subroutine
    
    subroutine get_rotation_matrix2(v0,rot)
    
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
    
    ! do the third rotation (about y-axis)
    a4 = matmul(rot3,a3)
    b4 = matmul(rot3,b3)
    
    if (a4(1)*b4(2)-a4(2)*b4(1) .gt. 0) then ! rotate about y-axis (any axis should work) by 180 degrees if facing down
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
    
    subroutine contour_integral(lambda,area_facs2,rot1,rot2,v0,prop2,x3,y3,z3)
    
        ! computes the area factor aka scalar fraunhofer diffraction pattern using contour integral method
    
    real(8), intent(in) :: lambda ! wavelength
    complex(8), dimension(:,:), allocatable, intent(inout) :: area_facs2
    real(8) bin_vec_size_k ! distance to far-field bins - important for accurate phase
    real(8) kxx, kyy
    real(8) waveno ! wave number
    real(8), intent(in) :: v0(1:3,1:3) ! vertices of facet after translating to com system
    real(8), intent(in) :: rot1(1:3,1:3) ! rotation matrix #1
    real(8), intent(in) :: rot2(1:3,1:3) ! rotation matrix (about x-axis)
    real(8) rot(1:3,1:3) ! update to new rotation matrix which aligns with incidence
    real(8) v1(1:3,1:3) ! vertices of facet after rotating
    integer(8) j, i1, i2
    real(8), dimension(:,:), allocatable, intent(in) :: x3, y3, z3
    real(8), dimension(:,:), allocatable :: kxxs, kyys, bin_vec_size_ks
    real(8), intent(in) :: prop2(1:3)
    real(8) kInc(1:3) ! incoming propagation vector * wavenumber
    real(8) x(1:3), y(1:3) ! x and y coordinates of aperture after rotating into x-y plane
    real(8) m(1:3), n(1:3) ! gradient and inverse gradient between vertex j and j + 1
    real(8) mj, nj, xj, yj, xj_plus1, yj_plus1
    real(8) dx, dy
    real(8) alpha, beta, delta, omega1, omega2
    real(8) sumre, sumim
    real(8) delta1, delta2
    real(8) nf

    ! allocate
    allocate(kxxs(1:size(x3,2),1:size(x3,1)))
    allocate(kyys(1:size(x3,2),1:size(x3,1)))
    allocate(bin_vec_size_ks(1:size(x3,2),1:size(x3,1)))

    waveno = 2.0*pi/lambda
    
    rot = matmul(rot2,rot1) ! rotate vertices
    v1 = matmul(rot,v0) ! translate vertices to new com
    
    kInc = prop2*waveno ! kincrot in matlab
    
    ! ===========================
    ! aperture moved into xy plane
    ! now do contour integral
    ! ===========================
    
    x(1:3) = v1(1,1:3)
    y(1:3) = v1(2,1:3)
    
    do j = 1, 3 ! loop over vertices, compute gradient
        if(j .eq. 3) then
            m(j) = (y(1) - y(j)) / (x(1) - x(j))
        else
            m(j) = (y(j+1) - y(j)) / (x(j+1) - x(j))
        end if
    end do
    
    n = 1/m ! get  inverse gradient
    area_facs2 = 0
    
    do i2 = 1, size(x3,2)
        do i1 = 1, size(x3,1)
            bin_vec_size_ks(i2,i1) = waveno*sqrt(x3(i1,i2)**2 + y3(i1,i2)**2 + z3(i1,i2)**2) ! distance in k-space to far-field bin
            kxxs(i2,i1) = kinc(1) - waveno*waveno*x3(i1,i2)/bin_vec_size_ks(i2,i1) ! kx' in derivation
            kyys(i2,i1) = kinc(2) - waveno*waveno*y3(i1,i2)/bin_vec_size_ks(i2,i1) ! ky' in derivation

            if (abs(kxxs(i2,i1)) .lt. 1e-6) kxxs(i2,i1) = 1e-6
            if (abs(kyys(i2,i1)) .lt. 1e-6) kyys(i2,i1) = 1e-6

        end do
    end do

    ! omp stuff below works but loop isnt big enough to give a speedup
    !!$OMP PARALLEL DEFAULT(private) num_threads(1) SHARED(xfar,x3,y3,z3,kinc,waveno,m,n,area_facs2,x,y)
    !!$OMP DO
    ! contour summation
    do j = 1, 3
        mj = m(j)
        nj = n(j)
        xj = x(j)
        yj = y(j)

        if(abs(mj) .gt. HUGE(1D0)) mj = 1e6 ! fix infinite gradient (these terms have little contribution becaese x_(j+1) = x(j) and they cancel)
        if(abs(nj) .gt. HUGE(1D0)) nj = 1e6 ! fix infinite 1/gradient

        if(j .eq. 3) then
            xj_plus1 = x(1)
            yj_plus1 = y(1)
        else
            xj_plus1 = x(j+1)
            yj_plus1 = y(j+1)
        end if

        dx = xj_plus1 - xj
        dy = yj_plus1 - yj

        do i2 = 1, size(x3,2)
            do i1 = 1, size(x3,1)

                ! 11/7/23 reduced number of flops here with some simplification and saving of values

                ! bin_vec_size_k = waveno*sqrt(x3(i1,i2)**2 + y3(i1,i2)**2 + z3(i1,i2)**2) ! distance in k-space to far-field bin
                bin_vec_size_k = bin_vec_size_ks(i2,i1) ! distance in k-space to far-field bin
                ! kxx = kinc(1) - waveno*waveno*x3(i1,i2)/bin_vec_size_k ! kx' in derivation
                ! kyy = kinc(2) - waveno*waveno*y3(i1,i2)/bin_vec_size_k ! ky' in derivation

                kxx = kxxs(i2,i1) ! kx' in derivation
                kyy = kyys(i2,i1) ! ky' in derivation

                delta = kxx*xj+kyy*yj
                delta1 = kyy*mj+kxx
                delta2 = kxx*nj+kyy
                omega1 = dx*delta1
                omega2 = dy*delta2

                alpha = 1/(2*kyy*delta1) ! denominator in p summand
                beta =  1/(2*kxx*delta2) ! denominator in q summand

                sumim = +alpha*(cos(delta)-cos(delta+omega1)) - beta*(cos(delta)-cos(delta+omega2))
                sumre = -alpha*(sin(delta)-sin(delta+omega1)) + beta*(sin(delta)-sin(delta+omega2))

                area_facs2(i1,i2) = area_facs2(i1,i2) + cmplx(cos(bin_vec_size_k),sin(bin_vec_size_k)) * cmplx(sumre,sumim) / lambda

            end do
        end do

        



    end do
    !!$OMP END PARALLEL

    ! loop to ignore bins which have numerical errors ie. when kyy -> 0 or kxx -> 0
    ! this appears to just be a flaw in the theory -> merits further investigation
    ! abs(area_facs2) appears to take value in the range 0 to 1 (check this!)
    ! if it goes over this, rescale it
    ! do i2 = 1, size(x3,2)
    !     do i1 = 1, size(x3,1)
    !         if (abs(area_facs2(i1,i2)) .gt. 10) then
    !             nf = abs(area_facs2(i1,i2))
    !             area_facs2(i1,i2) = area_facs2(i1,i2) / nf
    !         end if
    !     end do
    ! end do

    end subroutine
    
    subroutine diffraction( ampl,       & ! the outgoing amplitude matrix on the particle surface
                            v_in,       & ! the vertices of the aperture
                            prop0,      & ! the outgoing beam propagation direction (in lab system)
                            perp0,      & ! the outgoing beam perpendicular field direction (in lab system)
                            xfar,       & ! the x coordinate of each far-field bin (meshgrid style)
                            yfar,       & ! the y coordinate of each far-field bin (meshgrid style)
                            zfar,       & ! the z coordinate of each far-field bin (meshgrid style)
                            lambda,     & ! wavelength
                            amplC11s,   & ! the far-field diffraction amplitude matrix to be calculated (1,1)
                            amplC12s,   & ! the far-field diffraction amplitude matrix to be calculated (1,2)
                            amplC21s,   & ! the far-field diffraction amplitude matrix to be calculated (2,1)
                            amplC22s,   & ! the far-field diffraction amplitude matrix to be calculated (2,2)
                            phi_vals,   & ! phi values
                            theta_vals)   ! theta values
    
        ! main diffraction subroutine
        ! takes in a few arguments and outputs the far-field scattering pattern
    
    complex(8), intent(inout) :: ampl(1:2,1:2) ! amplitude matrix to pass to diffraction sr
    real(8), intent(in) :: perp0(1:3) ! perp field direction to pass to diffraction sr
    real(8), intent(in) :: prop0(1:3) ! outgoing propagation direction to pass to diffraction sr
    real(8), intent(in) :: v_in(1:3,1:3) ! vertices of facet to pass to diffraction sr
    real(8), dimension(:,:), allocatable, intent(in) :: xfar, yfar, zfar ! far-field bin positions
    real(8), intent(in) :: lambda ! wavelength
    complex(8), dimension(:,:), allocatable, intent(inout) :: amplC11s, amplC12s, amplC21s, amplC22s
    real(8), dimension(:), allocatable, intent(in) :: phi_vals
    real(8), dimension(:), allocatable, intent(in) :: theta_vals

    complex(8), dimension(:,:), allocatable :: area_facs2
    real(8) diff_ampl(1:2,1:2)
    real(8) rot(1:3,1:3), rot2(1:3,1:3)
    real(8) m(1:3)
    real(8) k(1:3)
    logical anti_parallel
    real(8) incidence(1:3), incidence2(1:3)
    real(8) v(1:3,1:3) ! vertices of facet (after lr flip)
    real(8) hc(1:3)
    real(8) evo2(1:3)
    real(8) rot4(1:2,1:2)
    integer i, j
    real(8) temp_rot1(1:2,1:2) 
    real(8) temp_vec3(1:3)
    complex(8) ampl_temp2(1:2,1:2)
    real(8) v0(1:3,1:3) ! vertices of facet after translating to com system
    real(8) com(1:3) ! aperture centre of mass
    real(8) prop1(1:3), perp1(1:3)
    real(8) angle
    real(8) prop2(1:3), perp2(1:3), e_par2(1:3)
    real(8), dimension(:,:), allocatable:: x3, y3, z3 ! far-field bin positions
    real(8) cos_rot, sin_rot
    real(8) bin_vec_size ! distance to far-field bins - important for accurate phase
    real(8) phi
    real(8) my_k(1:3)
    real(8), dimension(:), allocatable :: my_phi_vals ! phi values in aperture system for this aperture only
    logical success
    integer(8) theta_i

    allocate(my_phi_vals(1:size(phi_vals,1)))
    allocate(area_facs2(1:size(xfar,1),1:size(xfar,2)))
    
    ! flip vertex order to make normals face the other way (bit of a bodge) to do: find a way to remove this
    call flip_vert_lr(v)

    v(1:3,1:3) = v_in(1:3,1:3) ! store vertices
    
    com = sum(v,2)/3 ! get centre of mass
    v0(1,1:3) = v(1,1:3) - com(1) ! translate aperture to centre of mass system
    v0(2,1:3) = v(2,1:3) - com(2)
    v0(3,1:3) = v(3,1:3) - com(3)
    
    call get_rotation_matrix2(v0,rot)

    prop1 = matmul(rot,prop0)
    perp1 = matmul(rot,perp0)
    
    angle = -atan2(prop1(2),prop1(1))
    
    rot2(1,1) = cos(angle)
    rot2(1,2) = -sin(angle)
    rot2(1,3) = 0
    rot2(2,1) = sin(angle)
    rot2(2,2) = cos(angle)
    rot2(2,3) = 0
    rot2(3,1) = 0
    rot2(3,2) = 0
    rot2(3,3) = 1
    prop1 = matmul(rot,prop0)
    perp1 = matmul(rot,perp0)
    prop2 = matmul(rot2,prop1)
    perp2 = matmul(rot2,perp1)

    call cross(perp2,prop2,e_par2) ! to do: turn on normalisation here
    
    if(e_par2(3) .gt. 0.0001) then
        anti_parallel = .true.
    else
        anti_parallel = .false.
    end if

    ! required because karczewski theory requires perpendicular field direction to be along +ive y in aperture system
    ! this bodge reverses the amplitude matrix components if this is the case
    if(anti_parallel) ampl = -ampl

    incidence = (/0, 0, -1/) ! incident illumination)
    incidence2 = matmul(rot2,matmul(rot,incidence))

    call transform_bins2(xfar,yfar,zfar,com,matmul(rot2,rot),x3,y3,z3)
    
    ! bodge to fix numerical error for rot4 rotation matrix
    ! if incidence direction is almost parallel to bin vector, phi is not well defined
    ! therefore, get phi for the direct forwards direction from the phi values of the closest bin with theta value greater than 1
    success = .false.
    i = 0
    do while (i .lt. size(theta_vals,1) .and. .not. success)
        i = i + 1
        if(theta_vals(i) .gt. 1*pi/180) then
            theta_i = i
            success = .true.
        end if
    end do

    do i = 1, size(xfar,1)
        do j = 1, size(xfar,2)

            phi = phi_vals(i)
            ! phi = 0d0

            ! karczewski theory
            call karczewski(diff_ampl,m,k,prop2,x3(i,j), y3(i,j),z3(i,j))

            ! get the vector perpendicular to the scattering plane as viewed in the aperture system
            if(abs(dot_product(incidence2,k)) .lt. 0.999) then ! if bin is not in direct forwards
                call cross(incidence2,k,hc)
            else ! rotate vector perp to scattering plane (y-z) into aperture system
                ! call random_rotation(cos_rot,sin_rot)
                ! hc = matmul(rot2,matmul(rot,(/1d0,0d0,0d0/)))
                ! hc = matmul(rot2,matmul(rot,(/cos(my_phi_vals(i)),sin(my_phi_vals(i)),0d0/)))
                ! hc = matmul((/cos(phi),sin(phi),0d0/),matmul(rot2,rot))
                my_k = (/x3(i,theta_i),y3(i,theta_i),z3(i,theta_i)/)/sqrt(x3(i,theta_i)**2 + y3(i,theta_i)**2 + z3(i,theta_i)**2)
                call cross(incidence2,my_k,hc)
                ! hc(1) = -hc(1)
                ! hc(2) = -hc(2)
                ! hc = (/cos(phi),sin(phi),0d0/)
                ! hc = matmul(rot2,matmul(rot,(/cos_rot,sin_rot,0d0/)))
            end if

            call cross(k,m,evo2)

            ! make rotation matrix to rotate from KW plane to scattering plane
            rot4(1,1) = dot_product(hc,m)
            rot4(2,2) = dot_product(hc,m)
            rot4(1,2) = -dot_product(hc,evo2)
            rot4(2,1) = +dot_product(hc,evo2)          

            ! get the premultiplication matrix to rotate the incident amplitude matrix into the scattering plane
            ! get x-y normalised e-perp direction of scattering plane (used lab system here)
            ! temp_vec3 = (/xfar(i,j),yfar(i,j),zfar(i,j)/) / sqrt(xfar(i,j)**2 + yfar(i,j)**2)
            
            ! bin_vec_size = sqrt(xfar(i,j)**2 + yfar(i,j)**2 + zfar(i,j)**2)
            ! ! if scattering bin is close to direct forwards, get perp vector from phi value
            ! if(abs(zfar(i,j)/bin_vec_size) .gt. 0.99999) then
            !     ! call random_rotation(cos_rot,sin_rot)
            !     ! temp_vec3 = (/cos_rot,-sin_rot,0d0/)
            !     ! temp_vec3 = (/1d0,0d0,0d0/)
            !     phi = phi_vals(i)
            !     temp_vec3 = (/cos(phi),sin(phi),0d0/)

            ! end if

            temp_vec3 = (/cos(phi),sin(phi),0d0/)
            ! temp_vec3 = (/0d0,1d0,0d0/)

            ! get rotation matrix to rotate from incidence plane (x-y) about incidence direction to scattering plane
            call getRotationMatrix(temp_rot1,-temp_vec3(2),temp_vec3(1),0d0,1d0,0d0,0d0,0d0,0d0,-1d0)

            ! bodge (rotation was wrong way)
            temp_rot1 = transpose(temp_rot1)

            ! apply rotation matrices
            ampl_temp2 = matmul(rot4,matmul(diff_ampl,matmul(ampl,temp_rot1)))
            ! ampl_temp2 = matmul(diff_ampl,matmul(ampl,temp_rot1))
            ! ampl_temp2 = ampl ! disable this!!

            amplC11s(i,j) = ampl_temp2(1,1)
            amplC12s(i,j) = ampl_temp2(1,2)
            amplC21s(i,j) = ampl_temp2(2,1)
            amplC22s(i,j) = ampl_temp2(2,2)

        end do
    end do

    ! scalar fraunhofer integral, output contained in area_facs2
    call contour_integral(lambda,area_facs2,rot,rot2,v0,prop2,x3,y3,z3)

    amplC11s = amplC11s * area_facs2
    amplC12s = amplC12s * area_facs2
    amplC21s = amplC21s * area_facs2
    amplC22s = amplC22s * area_facs2

    end subroutine
    
    subroutine diff_main(   beam_outbeam_tree,          & ! tree of outgoing beams
                            beam_outbeam_tree_counter,  & ! total number of outoing beams in beam tree
                            lambda,                     & ! wavelength
                            ampl_far_beam11,            & ! far-field diffracted beam amplitude matrix (1,1)
                            ampl_far_beam12,            & ! far-field diffracted beam amplitude matrix (1,2)
                            ampl_far_beam21,            & ! far-field diffracted beam amplitude matrix (2,1)
                            ampl_far_beam22,            & ! far-field diffracted beam amplitude matrix (2,2)
                            theta_vals,                 & ! theta values (in rad)
                            phi_vals,                   & ! phi values (in rad)
                            ext_diff_outbeam_tree,      & ! tree of outgoing external diffraction beams
                            ampl_far_ext_diff11,        & ! far-field external diffraction amplitude matrix (1,1)
                            ampl_far_ext_diff12,        & ! far-field external diffraction amplitude matrix (1,2)
                            ampl_far_ext_diff21,        & ! far-field external diffraction amplitude matrix (2,1)
                            ampl_far_ext_diff22,        & ! far-field external diffraction amplitude matrix (2,2)
                            is_multithreaded)             ! whether multithreaded diffraction should be enabled
    
        ! sr diff_main is the main shell for diffraction of all beams + external diffraction at a fixed orientation

    type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
    integer(8), intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
    real(8), intent(in) :: lambda ! wavelength
    real(8), dimension(:), allocatable, intent(out) :: theta_vals, phi_vals
    real(8), dimension(:,:), allocatable :: xfar, yfar, zfar ! far-field bin positions
    logical, intent(in) :: is_multithreaded
    integer j
    complex(8), dimension(:,:), allocatable, intent(out) :: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! total
    complex(8), dimension(:,:), allocatable, intent(out) :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! total
    complex(8), dimension(:,:), allocatable :: amplC11s, amplC12s, amplC21s, amplC22s ! summand
    complex(8) ampl(1:2,1:2) ! amplitude matrix to pass to diffraction sr
    real(8) perp0(1:3) ! perp field direction to pass to diffraction sr
    real(8) prop0(1:3) ! outgoing propagation direction to pass to diffraction sr
    real(8) v(1:3,1:3) ! vertices of facet to pass to diffraction sr
    complex(8), dimension(:,:), allocatable :: area_facs2
    type(outbeamtype), dimension(:), allocatable, intent(in) :: ext_diff_outbeam_tree
    real(8) start, finish ! cpu timing variables
    real(8) progressReal
    integer progressInt 

    start = omp_get_wtime()

    call make_far_field_bins(xfar,yfar,zfar,theta_vals,phi_vals) ! get meshgrid-style far-field bins x, y, z

    ! allocate some far-field bins and init
    allocate(ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_ext_diff11(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_ext_diff12(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_ext_diff21(1:size(xfar,1),1:size(xfar,2)))
    allocate(ampl_far_ext_diff22(1:size(xfar,1),1:size(xfar,2)))
    allocate(area_facs2(1:size(xfar,1),1:size(xfar,2)))
    allocate(amplC11s(1:size(xfar,1),1:size(xfar,2)))
    allocate(amplC12s(1:size(xfar,1),1:size(xfar,2)))
    allocate(amplC21s(1:size(xfar,1),1:size(xfar,2)))
    allocate(amplC22s(1:size(xfar,1),1:size(xfar,2)))
    
    area_facs2 = 0
    ampl_far_beam11 = 0
    ampl_far_beam12 = 0
    ampl_far_beam21 = 0
    ampl_far_beam22 = 0
    ampl_far_ext_diff11 = 0
    ampl_far_ext_diff12 = 0
    ampl_far_ext_diff21 = 0
    ampl_far_ext_diff22 = 0
    progressInt = 0

    ! print*,'========== start sr diff_main'
    
    call trim_outbeam_tree(beam_outbeam_tree,beam_outbeam_tree_counter) ! removes very low energy outbeams from the beam tree

    if (is_multithreaded) then ! multi-threaded diffraction

        ! print*,'multithreading: enabled'
        ! print*,'max threads:',omp_get_max_threads()
        !$OMP PARALLEL num_threads(omp_get_max_threads()) PRIVATE(amplC11s,amplC12s,amplC21s,amplC22s,ampl,perp0,prop0,v)
        !$OMP DO
        do j = 1, beam_outbeam_tree_counter

            ! progressReal = j*100/beam_outbeam_tree_counter          ! percent completion
            ! if(int(progressReal) .gt. progressInt .and. mod(int(progressReal),10) .eq. 0) then  ! if at least 1% progress has been made
            !     progressInt = int(progressReal)           ! update progress counter
            ! print*,'Completion: ',progressInt,'%'     ! print progress
            ! end if
    
            ampl(1:2,1:2) = beam_outbeam_tree(j)%ampl(1:2,1:2)
            perp0(1:3) = beam_outbeam_tree(j)%vk7(1:3)
            prop0(1:3) = beam_outbeam_tree(j)%prop_out(1:3)
            v(1:3,1:3) = beam_outbeam_tree(j)%verts(1:3,1:3)
        
            ! get far-field contribution from this outbeam
            call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC11s,amplC12s,amplC21s,amplC22s,phi_vals,theta_vals)
        
            !$OMP CRITICAL
            ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) + amplC11s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) + amplC12s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) + amplC21s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) + amplC22s(1:size(xfar,1),1:size(xfar,2))
            !$OMP END CRITICAL
        
        end do
        
        ! if(omp_get_thread_num() .eq. 0) then
        !     print*,'end diff beam loop...'
        !     print*,'start ext diff loop...'
        ! end if
        
        ! progressInt = 0
        !$OMP DO
        do j = 1, size(ext_diff_outbeam_tree,1)

            ! progressReal = j*100/size(ext_diff_outbeam_tree,1)        ! percent completion
            ! if(int(progressReal) .gt. progressInt .and. mod(int(progressReal),10) .eq. 0) then  ! if at least 1% progress has been made
            !     progressInt = int(progressReal)           ! update progress counter
            ! print*,'Completion: ',progressInt,'%'     ! print progress
            ! end if

            ampl(1:2,1:2) = ext_diff_outbeam_tree(j)%ampl(1:2,1:2)
            perp0(1:3) = ext_diff_outbeam_tree(j)%vk7(1:3)
            prop0(1:3) = ext_diff_outbeam_tree(j)%prop_out(1:3)
            v(1:3,1:3) = ext_diff_outbeam_tree(j)%verts(1:3,1:3)
        
            call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC11s,amplC12s,amplC21s,amplC22s,phi_vals,theta_vals)

            !$OMP CRITICAL
            ampl_far_ext_diff11(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff11(1:size(xfar,1),1:size(xfar,2)) + amplC11s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff12(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff12(1:size(xfar,1),1:size(xfar,2)) + amplC12s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff21(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff21(1:size(xfar,1),1:size(xfar,2)) + amplC21s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff22(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff22(1:size(xfar,1),1:size(xfar,2)) + amplC22s(1:size(xfar,1),1:size(xfar,2))
            !$OMP END CRITICAL
                
        end do
        !$OMP END PARALLEL


    else ! single-threaded diffraction

        ! print*,'multithreading: disabled'

        do j = 1, beam_outbeam_tree_counter
        ! do j = 1, 0 ! disable beam diffraction

            ! progressReal = j*100/beam_outbeam_tree_counter          ! percent completion
            ! if(int(progressReal) .gt. progressInt .and. mod(int(progressReal),10) .eq. 0) then  ! if at least 1% progress has been made
            !     progressInt = int(progressReal)           ! update progress counter
            ! print*,'Completion: ',progressInt,'%'     ! print progress
            ! end if
    
            ampl(1:2,1:2) = beam_outbeam_tree(j)%ampl(1:2,1:2)
            perp0(1:3) = beam_outbeam_tree(j)%vk7(1:3)
            prop0(1:3) = beam_outbeam_tree(j)%prop_out(1:3)
            v(1:3,1:3) = beam_outbeam_tree(j)%verts(1:3,1:3)
        
            ! ! get far-field contribution from this outbeam
            call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC11s,amplC12s,amplC21s,amplC22s,phi_vals,theta_vals)
        
            ! !$OMP CRITICAL
            ! ! sum
            ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) + amplC11s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) + amplC12s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) + amplC21s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) + amplC22s(1:size(xfar,1),1:size(xfar,2))
            ! !$OMP END CRITICAL
        
        end do

        ! if (is_multithreaded) then
            !!$OMP END PARALLEL
        ! end if
        ! stop
        !!$OMP DO
        ! do j = 1, beam_outbeam_tree_counter
        ! ! do j = 1, 100
        !     ! print*,'j = ',j

        !     progressReal = j*100/beam_outbeam_tree_counter          ! percent completion
        !     if(int(progressReal) .gt. progressInt .and. mod(int(progressReal),10) .eq. 0) then  ! if at least 1% progress has been made
        !         progressInt = int(progressReal)           ! update progress counter
        !     print*,'Completion: ',progressInt,'%'     ! print progress
        !     end if

        !     ampl(1:2,1:2) = beam_outbeam_tree(j)%ampl(1:2,1:2)
        !     perp0(1:3) = beam_outbeam_tree(j)%vk7(1:3)
        !     prop0(1:3) = beam_outbeam_tree(j)%prop_out(1:3)
        !     v(1:3,1:3) = beam_outbeam_tree(j)%verts(1:3,1:3)
        
        !     ! get far-field contribution from this outbeam
        !     call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC11s,amplC12s,amplC21s,amplC22s)
        
        !     !!$OMP CRITICAL
        !     ! sum
        !     ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam11(1:size(xfar,1),1:size(xfar,2)) + amplC11s(1:size(xfar,1),1:size(xfar,2))
        !     ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam12(1:size(xfar,1),1:size(xfar,2)) + amplC12s(1:size(xfar,1),1:size(xfar,2))
        !     ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam21(1:size(xfar,1),1:size(xfar,2)) + amplC21s(1:size(xfar,1),1:size(xfar,2))
        !     ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) = ampl_far_beam22(1:size(xfar,1),1:size(xfar,2)) + amplC22s(1:size(xfar,1),1:size(xfar,2))
        !     !!$OMP END CRITICAL
        
        ! end do
        
        !!$OMP END PARALLEL
    
        ! print*,'end diff beam loop...'
        
        ! print*,'start ext diff loop...'
        
        ! print*,'size(ext_diff_outbeam_tree,1)',size(ext_diff_outbeam_tree,1)

        ! progressInt = 0
        do j = 1, size(ext_diff_outbeam_tree,1)
        ! do j = 1, 0 ! disable external diff
        ! do j = 90, 90 ! test

            ! progressReal = j*100/size(ext_diff_outbeam_tree,1)        ! percent completion
            ! if(int(progressReal) .gt. progressInt .and. mod(int(progressReal),10) .eq. 0) then  ! if at least 1% progress has been made
            !     progressInt = int(progressReal)           ! update progress counter
            ! print*,'Completion: ',progressInt,'%'     ! print progress
            ! end if

            ampl(1:2,1:2) = ext_diff_outbeam_tree(j)%ampl(1:2,1:2)
            perp0(1:3) = ext_diff_outbeam_tree(j)%vk7(1:3)
            prop0(1:3) = ext_diff_outbeam_tree(j)%prop_out(1:3)
            v(1:3,1:3) = ext_diff_outbeam_tree(j)%verts(1:3,1:3)

            ! if(abs(ampl(1,1)) .gt. 1e-5 .or. abs(ampl(1,2)) .gt. 1e-5 .or. abs(ampl(2,1)) .gt. 1e-5 .or. abs(ampl(2,2)) .gt. 1e-5) then
        
            call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC11s,amplC12s,amplC21s,amplC22s,phi_vals,theta_vals)
        
            ampl_far_ext_diff11(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff11(1:size(xfar,1),1:size(xfar,2)) + amplC11s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff12(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff12(1:size(xfar,1),1:size(xfar,2)) + amplC12s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff21(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff21(1:size(xfar,1),1:size(xfar,2)) + amplC21s(1:size(xfar,1),1:size(xfar,2))
            ampl_far_ext_diff22(1:size(xfar,1),1:size(xfar,2)) = ampl_far_ext_diff22(1:size(xfar,1),1:size(xfar,2)) + amplC22s(1:size(xfar,1),1:size(xfar,2))
        
            ! end if
        
        end do
    end if


















    ! call kmp_set_stacksize_s(109715200)
    ! print*,'stack size: (bytes)',kmp_get_stacksize_s()
    
    
    ! CALL get_environment_variable("OMP_STACKSIZE", stack_size)
    ! print*,'stack size: ',trim(stack_size)

    ! print*,'max omp threads: ',OMP_GET_MAX_THREADS()

    !if (is_multithreaded) then
        !!OMP_GET_MAX_THREADS
        !!$OMP PARALLEL num_threads(omp_get_max_threads()) PRIVATE(amplC11s,amplC12s,amplC21s,amplC22s,ampl,perp0,prop0,v)
    !end if
    !!$OMP PARALLEL num_threads(4)
    ! stop
    !print*,'hello from thread #',OMP_GET_THREAD_NUM()

    !!$OMP DO
    
    
    ! print*,'end ext diff loop...'
    
    finish = omp_get_wtime()
    ! print*,'=========='
    print'(A,f16.8,A)',"end diffraction - total time taken: ",finish-start," secs"
    ! print*,'=========='

    ! print*,'========== end sr diff_main'

    end subroutine
    
    subroutine make_far_field_bins(x, y, z,theta_vals, phi_vals, & ! non-optional args
        theta_start_in, theta_end_in, phi_start_in, phi_end_in) ! optional args
    
        ! makes the far field x, y, z bins at which the diffracted field will be computed
    
    integer(8) theta_dim, phi_dim ! theta and phi discretisation
    real(8), dimension(:,:), allocatable, intent(out) :: x, y, z ! far-field bin positions
    real(8), optional, intent(in) :: theta_start_in, theta_end_in, phi_start_in, phi_end_in ! 
    
    real(8) r ! distance to far-field
    real(8), dimension(:), allocatable, intent(out) :: theta_vals, phi_vals
    real(8), dimension(:,:), allocatable :: theta_vals_mesh, phi_vals_mesh
    real(8) theta_start, theta_end, phi_start, phi_end ! 
    integer i

    call read_theta_vals(theta_vals)
    ! stop
    call read_phi_vals(phi_vals)

    phi_dim = size(phi_vals,1)
    theta_dim = size(theta_vals,1)

    ! allocate
    ! allocate(theta_vals(1:theta_dim))
    ! allocate(phi_vals(1:phi_dim))
    allocate(x(1:phi_dim,1:theta_dim))
    allocate(y(1:phi_dim,1:theta_dim))
    allocate(z(1:phi_dim,1:theta_dim))

    ! stop
    
    ! convert theta vals to rad
    theta_vals = theta_vals*pi/180d0
    ! stop
    ! convert phi vals to rad
    phi_vals = phi_vals*pi/180d0

    ! numerical fixes due to divide by 0 in contour integral
    do i = 1, size(phi_vals)
        if(abs(phi_vals(i)*180/pi) .lt. 0.000001) phi_vals(i) = phi_vals(i) + 0.00001*pi/180
        if(abs(phi_vals(i)*180/pi - 360.0) .lt. 0.000001) phi_vals(i) = phi_vals(i) + 0.00001*pi/180
    end do

    ! stop
    
    ! convert 1-d arrays to 2-d arrays
    call meshgrid_real(theta_vals, phi_vals, theta_vals_mesh, phi_vals_mesh)

    ! DO NOT SET THIS TOO BIG - if waveno*r > 1e6, trig function precision fails vs matlab
    r = 1e4 ! distance to far-field
    
    ! convert to cartesian
    x = r*sin(theta_vals_mesh)*cos(phi_vals_mesh)
    y = r*sin(theta_vals_mesh)*sin(phi_vals_mesh)
    z = -r*cos(theta_vals_mesh) ! to do: check that this is ok to do (swaps forwards -> backwards)
    
    ! print*,'x(1,66)',x(1,66)
    ! print*,'x(181,66)',x(181,66)
    ! print*,'y(1,66)',y(1,66)
    ! print*,'y(181,66)',y(181,66)
    ! print*,'z(1,66)',z(1,66)
    ! print*,'z(181,66)',z(181,66)
    ! print*,'theta_vals_mesh(1,66)',theta_vals_mesh(1,66)
    ! print*,'theta_vals_mesh(181,66)',theta_vals_mesh(181,66)    
    ! print*,'phi_vals_mesh(1,66)',phi_vals_mesh(1,66)
    ! print*,'phi_vals_mesh(181,66)',phi_vals_mesh(181,66)
    ! print*,'phi_vals(1)',phi_vals(1)
    ! print*,'phi_vals(181)',phi_vals(181)
    ! stop

    end subroutine
        
    subroutine read_theta_vals(theta_vals)

        ! sr read_theta,vals reads and makes the phi values
        ! taken from sr anglesplit

        ! important values:
        !   k is the total number of bins
        !   anglesout are the centres of each angular bin - size (1:k)
        !   anglesteps are the sizes of each angular bin - size (1:k)
        real(8), dimension(:), allocatable, intent(out) :: theta_vals
        integer nlines, i, io, j, jmax
        real(8), dimension(:) ,allocatable :: anglesint ,splitsint
        real(8), dimension(:) ,allocatable :: anglesout, anglesteps
        integer kint

        if(allocated(anglesout)) deallocate(anglesout)
        if(allocated(anglesteps)) deallocate(anglesteps)
        allocate( anglesout(10000), anglesteps(10000) )    ! set up arrays for anglesint and splitsint
        nlines = 0  ! initialise line counter
        
        open(119, file = "theta_vals.txt", status = 'old', action = 'read')  ! open input file
        
        do  ! read in number of lines in input file
            read(119,*,iostat=io)
            if (io/=0) exit
            nlines = nlines + 1
        end do
        ! print*,'nlines = ',nlines
        allocate( anglesint(nlines) ,splitsint(nlines-1) )    ! set up arrays for anglesint and splitsint
        rewind(119)  ! rewind to top of input file
        
        do i = 1,nlines     ! read in anglesint and splitsint
            if(i .lt. nlines) then
                read(119,*) anglesint(i),splitsint(i)
            else
                read(119,*) anglesint(i)
            end if
        end do
        
        kint = 0
        do i = 1,nlines-1   ! for each line in the input file
            jmax = nint((anglesint(i+1)-anglesint(i))/splitsint(i))    ! compute the number of angle steps between one line and the next
            if(mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1)) .gt. 1e-3 .and. mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1)) .lt. 1-1e-3) then 
              print*,'bad angle choice, line:',i
              print*,'this should be close to integer:',(anglesint(i+1)-anglesint(i))/splitsint(i)
              print*,'should be less than 0.001',mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1))
              error stop     ! please check choice of bin splitting
            end if
            do j = 0,jmax
                if (i .ne. nlines-1 .and. j .eq. jmax) then
                    ! do nothign
                else
                ! if (j .eq. 0. .and. i .ne. 1 .and. i .ne. nlines-1) then ! skip duplicates where 2 lines meet
                ! else if (i .ne. nlines-1 .and. j .eq. jmax) then ! bodge
                ! else
                    kint = kint + 1
                    anglesout(kint) = anglesint(i) + j*splitsint(i)    ! central angle of each angular bin
                end if
            end do
        end do

        close(119)

        allocate(theta_vals(1:kint))

        do i = 1, kint
            theta_vals(i) = anglesout(i)
            ! print*,'i',i,' theta:', theta_vals(i)
        end do
        ! stop
        end subroutine

    subroutine read_phi_vals(phi_vals)

            ! sr read_phi_vals reads and makes the phi values
            ! taken from sr anglesplit

            ! important values:
            !   k is the total number of bins
            !   anglesout are the centres of each angular bin - size (1:k)
            !   anglesteps are the sizes of each angular bin - size (1:k)

            real(8), dimension(:), allocatable, intent(out) :: phi_vals
            integer nlines, i, io, j, jmax
            real(8), dimension(:) ,allocatable :: anglesint ,splitsint
            real(8), dimension(:) ,allocatable :: anglesout, anglesteps
            integer kint
    
            if(allocated(anglesout)) deallocate(anglesout)
            if(allocated(anglesteps)) deallocate(anglesteps)
            allocate( anglesout(10000), anglesteps(10000) )    ! set up arrays for anglesint and splitsint
            nlines = 0  ! initialise line counter
            
            open(119, file = "phi_vals.txt", status = 'old', action = 'read')  ! open input file
            
            do  ! read in number of lines in input file
                read(119,*,iostat=io)
                if (io/=0) exit
                nlines = nlines + 1
            end do
            ! print*,'nlines = ',nlines
            allocate( anglesint(nlines) ,splitsint(nlines-1) )    ! set up arrays for anglesint and splitsint
            rewind(119)  ! rewind to top of input file
            
            do i = 1,nlines     ! read in anglesint and splitsint
                if(i .lt. nlines) then
                    read(119,*) anglesint(i),splitsint(i)
                else
                    read(119,*) anglesint(i)
                end if
            end do
            
            kint = 0
            do i = 1,nlines-1   ! for each line in the input file
                jmax = nint((anglesint(i+1)-anglesint(i))/splitsint(i))    ! compute the number of angle steps between one line and the next
                if(mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1)) .gt. 1e-3 .and. mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1)) .lt. 1-1e-3) then 
                  print*,'bad angle choice, line:',i
                  print*,'this should be close to integer:',(anglesint(i+1)-anglesint(i))/splitsint(i)
                  print*,'should be less than 0.001',mod((anglesint(i+1)-anglesint(i))/splitsint(i),real(1))
                  error stop     ! please check choice of bin splitting
                end if
                do j = 0,jmax
                    if (j .eq. 0. .and. i .ne. 1 .and. i .ne. nlines-1) then ! skip duplicates where 2 lines meet
                    else
                        kint = kint + 1
                        anglesout(kint) = anglesint(i) + dble(j)*splitsint(i)    ! central angle of each angular bin
                    end if
                end do
            end do
    
            close(119)
    
            allocate(phi_vals(1:kint))
    
            do i = 1, kint
                phi_vals(i) = anglesout(i)
                ! print*,'i',i,' phi:', phi_vals(i)
            end do
    
            end subroutine        
        

    end module diff_mod