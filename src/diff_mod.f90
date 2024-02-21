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

    subroutine translate_far_field_bins(x,y,z,com,x1,y1,z1)
    
        ! translates and rotates bins to aperture system (loop method)
    
    real(8), dimension(:,:), allocatable, intent(in) :: x, y, z ! far-field bin positions
    real(8), intent(in) :: com(1:3) ! aperture centre of mass
    
    real(8), dimension(:,:), allocatable, intent(out) :: x1, y1, z1 ! far-field bin positions after translation to new com
    
    integer i, j
        
    ! allocate translated far-field bins
    allocate(x1(1:size(x,1),1:size(x,2)))
    allocate(y1(1:size(x,1),1:size(x,2)))
    allocate(z1(1:size(x,1),1:size(x,2)))

    ! translate far-field bins to new com
    do j = 1, size(x,2)
        do i = 1, size(x,1)
            x1(i,j) = x(i,j) - com(1)
            y1(i,j) = y(i,j) - com(2)
            z1(i,j) = z(i,j) - com(3)
        end do
    end do
    end subroutine

    subroutine rotate_far_field_bins(x,y,z,rot,x3,y3,z3,mapping,mapping2,theta_i)
    
        ! translates and rotates bins to aperture system (loop method)
    
    real(8), dimension(:,:), allocatable, intent(in) :: x, y, z ! far-field bin positions
    real(8), intent(in) :: rot(1:3,1:3) ! update to new rotation matrix which aligns with incidence
    real(8), dimension(:,:), allocatable, intent(out) :: x3, y3, z3 ! far-field bin positions
    integer(8), dimension(:,:), allocatable, intent(in) :: mapping
    integer(8), dimension(:), allocatable, intent(in) :: mapping2
    integer(8), intent(in), optional :: theta_i

    integer(8) i, j, ii
    
    ! allocate arrays to hold reshaped, rotated arrays
    allocate(x3(1:size(x,1),1:size(x,2)))
    allocate(y3(1:size(x,1),1:size(x,2)))
    allocate(z3(1:size(x,1),1:size(x,2)))
    
    ! rotate
    do j = 1, size(x,2)
        do ii = 1, mapping2(j)
        i = mapping(ii,j)
            x3(i,j) = rot(1,1)*x(i,j) + rot(1,2)*y(i,j) + rot(1,3)*z(i,j)
            y3(i,j) = rot(2,1)*x(i,j) + rot(2,2)*y(i,j) + rot(2,3)*z(i,j)
            z3(i,j) = rot(3,1)*x(i,j) + rot(3,2)*y(i,j) + rot(3,3)*z(i,j)
        end do
    end do
    
    if(present(theta_i)) then
        do j = 1, size(x,1)
            x3(j,theta_i) = rot(1,1)*x(j,theta_i) + rot(1,2)*y(j,theta_i) + rot(1,3)*z(j,theta_i)
            y3(j,theta_i) = rot(2,1)*x(j,theta_i) + rot(2,2)*y(j,theta_i) + rot(2,3)*z(j,theta_i)
            z3(j,theta_i) = rot(3,1)*x(j,theta_i) + rot(3,2)*y(j,theta_i) + rot(3,3)*z(j,theta_i)
        end do
    end if

    end subroutine
    
    subroutine karczewski(  diff_ampl, & ! polarisation matrix
                            m, & ! KW perp field vector
                            k, & ! KW par field vector
                            prop2, & ! outgoing propagation direction in aperture system
                            x3, y3, z3, & ! x, y, z coordinate of far-field bin
                            r) ! distance to bin vector
        
        ! sr karczewski computes the diff ampl matrix for polarisation of far-field diffraction at a far-field bin
        ! to do: add a debug mode with the numerical checks found in matlab code
    
    real(8), intent(out) :: diff_ampl(1:2,1:2)
    real(8), intent(out) :: m(1:3)
    real(8), intent(out) :: k(1:3)
    real(8), intent(in) :: prop2(1:3)
    real(8), intent(in) :: x3, y3, z3 ! far-field bin positions
    real(8), intent(in) :: r ! distance to bin vector
    real(8) bin_vec_size ! far-field bin distances
    real(8) big_kx, big_ky, big_kz ! outward propagation vector in aperture system
    real(8) frac
    real(8) a1m, b2m, a1e, b2e, b1m, a2e, a1em, a2em, b1em, b2em
    
    bin_vec_size = r ! distance to each far-field bin

    k(1) = x3/bin_vec_size ! propagation vector components for each bin vector in aperture system
    k(2) = y3/bin_vec_size
    k(3) = z3/bin_vec_size

    if(abs(k(2)) > 0.999999) k(2) = sign(dble(0.999999),k(2))
    
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

    ! ##### e-m theokry #####
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
    
    end subroutine
    
    subroutine get_rotation_matrix2(v0,rot)
    
        ! computes a rotation matrix for rotating into an aperture located in the xy plane
    
    real(8), dimension(:,:), allocatable, intent(in) :: v0 ! vertices of facet
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
    
    subroutine contour_integral(lambda,area_facs2,rot1,rot2,v0,prop2,x3,y3,z3,mapping,mapping2,nv)
    
        ! computes the area factor aka scalar fraunhofer diffraction pattern using contour integral method
    
    integer(8), dimension(:,:), allocatable, intent(in) :: mapping
    integer(8), dimension(:), allocatable, intent(in) :: mapping2
    real(8), intent(in) :: lambda ! wavelength
    complex(8), dimension(:,:), allocatable, intent(inout) :: area_facs2
    real(8) bin_vec_size_k ! distance to far-field bins - important for accurate phase
    real(8) kxx, kyy
    real(8) waveno ! wave number
    real(8), dimension(:,:), allocatable, intent(in) :: v0 ! vertices of facet after translating to com system
    real(8), intent(in) :: rot1(1:3,1:3) ! rotation matrix #1
    real(8), intent(in) :: rot2(1:3,1:3) ! rotation matrix (about x-axis)
    real(8) rot(1:3,1:3) ! update to new rotation matrix which aligns with incidence
    real(8), dimension(:,:), allocatable :: v1 ! vertices of facet after rotating
    integer(8) j, i1, i2, ii
    real(8), dimension(:,:), allocatable, intent(in) :: x3, y3, z3
    real(8), dimension(:,:), allocatable :: kxxs, kyys, bin_vec_size_ks
    real(8), intent(in) :: prop2(1:3)
    real(8) kInc(1:3) ! incoming propagation vector * wavenumber
    real(8), dimension(:), allocatable :: x, y ! x and y coordinates of aperture after rotating into x-y plane
    real(8), dimension(:), allocatable :: m, n ! gradient and inverse gradient between vertex j and j + 1
    real(8) mj, nj, xj, yj, xj_plus1, yj_plus1
    real(8) dx, dy
    real(8) alpha, beta, delta, omega1, omega2
    real(8) sumre, sumim
    real(8) delta1, delta2
    integer(8), intent(in) :: nv
    ! real(8) nf

    ! allocate
    allocate(x(1:nv))
    allocate(y(1:nv))
    allocate(m(1:nv))
    allocate(n(1:nv))
    allocate(v1(1:3,1:nv))
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
    
    x(1:nv) = v1(1,1:nv)
    y(1:nv) = v1(2,1:nv)

    do j = 1, nv ! loop over vertices, compute gradient
        if(j .eq. nv) then
            m(j) = (y(1) - y(j)) / (x(1) - x(j))
        else
            m(j) = (y(j+1) - y(j)) / (x(j+1) - x(j))
        end if
    end do
    
    n(:) = 1/m(:) ! get  inverse gradient
    area_facs2 = 0

    do i2 = 1, size(x3,2)
        do ii = 1, mapping2(i2) ! looping over the number of phi values to evaluate at, for this theta
            i1 = mapping(ii,i2)
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
    do j = 1, nv
        mj = m(j)
        nj = n(j)
        xj = x(j)
        yj = y(j)

        if(abs(mj) .gt. HUGE(1D0)) mj = 1e6 ! fix infinite gradient (these terms have little contribution becaese x_(j+1) = x(j) and they cancel)
        if(abs(nj) .gt. HUGE(1D0)) nj = 1e6 ! fix infinite 1/gradient
        if(abs(mj) .lt. 1d-9) nj = 1e6 ! fix infinite 1/gradient
        if(abs(nj) .lt. 1d-9) mj = 1e6 ! fix infinite 1/gradient

        if(j .eq. nv) then
            xj_plus1 = x(1)
            yj_plus1 = y(1)
        else
            xj_plus1 = x(j+1)
            yj_plus1 = y(j+1)
        end if

        dx = xj_plus1 - xj
        dy = yj_plus1 - yj

        do i2 = 1, size(x3,2)
            do ii = 1, mapping2(i2) ! looping over the number of phi values to evaluate at, for this theta
                i1 = mapping(ii,i2)
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

                area_facs2(i1,i2) = area_facs2(i1,i2) + cmplx(cos(bin_vec_size_k),sin(bin_vec_size_k),8) * cmplx(sumre,sumim,8) / lambda
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
                            amplC,      & ! the far-field diffraction amplitude matrix to be calculated
                            phi_vals,   & ! phi values
                            theta_vals, & ! theta values
                            fov,        &
                            job_params, &
                            outbeam, &
                            geometry)
    
        ! main diffraction subroutine
        ! takes in a few arguments and outputs the far-field scattering pattern
    
    complex(8), intent(inout) :: ampl(1:2,1:2) ! amplitude matrix to pass to diffraction sr
    real(8), intent(in) :: perp0(1:3) ! perp field direction to pass to diffraction sr
    real(8), intent(in) :: prop0(1:3) ! outgoing propagation direction to pass to diffraction sr
    real(8), intent(in) :: v_in(1:3,1:3) ! vertices of facet to pass to diffraction sr
    real(8), dimension(:,:), allocatable, intent(in) :: xfar, yfar, zfar ! far-field bin positions
    real(8), intent(in) :: lambda ! wavelength
    complex(8), dimension(:,:,:,:), allocatable, intent(inout) :: amplC
    real(8), dimension(:), allocatable, intent(in) :: phi_vals
    real(8), dimension(:), allocatable, intent(in) :: theta_vals
    real(8), intent(in) :: fov
    type(job_parameters_type), intent(in) :: job_params
    type(outbeamtype), intent(in) :: outbeam
    type(geometry_type), intent(in) :: geometry

    complex(8), dimension(:,:), allocatable :: area_facs2
    real(8) diff_ampl(1:2,1:2)
    real(8) rot(1:3,1:3), rot2(1:3,1:3)
    real(8) m(1:3)
    real(8) k(1:3)
    logical anti_parallel
    real(8) incidence(1:3), incidence2(1:3)
    real(8), dimension(:,:), allocatable :: v2 ! vertices of facet
    real(8) hc(1:3)
    real(8) evo2(1:3)
    real(8) rot4(1:2,1:2)
    integer(8) i, j, ii
    real(8) temp_rot1(1:2,1:2) 
    real(8) temp_vec3(1:3)
    complex(8) ampl_temp2(1:2,1:2)
    real(8), dimension(:,:), allocatable :: v20 ! vertices of facet after translating to com system
    real(8) com(1:3) ! aperture centre of mass
    real(8) prop1(1:3), perp1(1:3)
    real(8) angle
    real(8) prop2(1:3), perp2(1:3), e_par2(1:3)
    real(8), dimension(:,:), allocatable:: x3, y3, z3 ! translated, rotated far-field bin positions
    real(8), dimension(:,:), allocatable:: x1, y1, z1 ! translated, unrotated far-field bin positions
    ! real(8) cos_rot, sin_rot
    ! real(8) bin_vec_size ! distance to far-field bins - important for accurate phase
    real(8) phi
    real(8) my_k(1:3)
    real(8), dimension(:), allocatable :: my_phi_vals ! phi values in aperture system for this aperture only
    logical success
    integer(8) theta_i
    logical, dimension(:,:), allocatable :: is_in_fov ! whether or not each far-field bin is inside the field of view
    real(8) vec(1:3)
    real(8) cos_fov
    integer(8), dimension(:,:), allocatable :: mapping ! maps theta and phi index
    integer(8), dimension(:), allocatable :: mapping2 ! number of phi values at each theta value
    real(8), dimension(:,:), allocatable :: r1 ! far-field bin distances
    integer(8) nv ! number of vertices in outbeam face

    allocate(my_phi_vals(1:size(phi_vals,1)))
    allocate(area_facs2(1:size(xfar,1),1:size(xfar,2)))
    allocate(is_in_fov(1:size(xfar,1),1:size(xfar,2)))
    allocate(mapping(1:size(xfar,1),1:size(xfar,2)))
    allocate(mapping2(1:size(xfar,2)))
    allocate(r1(1:size(xfar,1),1:size(xfar,2)))

    is_in_fov(:,:) = .true. ! init
    amplC(:,:,:,:) = 0d0 ! init
    mapping(:,:) = 0 ! init
    mapping2(:) = 0 ! init
    cos_fov = cos(fov)

    nv = geometry%f(outbeam%fi)%nv ! get number of vertices in face

    allocate(v2(1:3,1:nv)) ! allocate arrays to hold vertex information
    allocate(v20(1:3,1:nv))

    do i = 1, nv ! for each vertex in face
        v2(:,i) = geometry%v(geometry%f(outbeam%fi)%vi(i),:) ! get information about vertices
    end do

    com = sum(v2,2)/nv ! get centre of mass
    do i = 1, nv
        v20(1,i) = v2(1,i) - com(1) ! translate aperture to centre of mass system
        v20(2,i) = v2(2,i) - com(2) ! translate aperture to centre of mass system
        v20(3,i) = v2(3,i) - com(3) ! translate aperture to centre of mass system
    end do

    call get_rotation_matrix2(v20,rot) ! get rotation matrix using first 2 vertices in face

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

    call cross(perp2,prop2,e_par2,.true.)
    
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

    call translate_far_field_bins(xfar,yfar,zfar,com,x1,y1,z1)

    r1(:,:) = sqrt(x1(:,:)**2 + y1(:,:)**2 + z1(:,:)**2) ! distance to each far-field bin

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

    if(job_params%is_fast_diff) then ! according to Jackson, most of the energy is confined to the region lambda/d
        do i = 1, size(xfar,1)
            do j = 1, size(xfar,2)
                vec(:) = (/x1(i,j),y1(i,j),z1(i,j)/) / r1(i,j) ! unit vector to far-field bin
                ! if(dot_product(vec,prop0) < cos_fov .and. j /= theta_i) then ! if in field of view, and not needed for numerical fix
                if(dot_product(vec,prop0) < cos_fov) then ! if in field of view, and not needed for numerical fix
                    is_in_fov(i,j) = .false. ! if dot product is less than threshold, bin is outside fov
                end if
            end do
        end do
    end if

    ! create mapping, to allow for vectorisation
    mapping2(:) = 0 ! init
    do j = 1, size(xfar,2) ! looping over theta
        ii = 0 ! reset counter for each value of theta
        do i = 1, size(xfar,1) ! looping over values of phi
            if(is_in_fov(i,j)) then ! if far-field is to be evaluated
                ii = ii + 1 ! update number of phi values to evaluate at, for this theta
                mapping(ii,j) = i ! save the position of this phi value
            end if
        end do ! end: looping over values of phi
        mapping2(j) = ii ! save the number of phi values to evaluate at
    end do

    call rotate_far_field_bins(x1,y1,z1,matmul(rot2,rot),x3,y3,z3,mapping,mapping2,theta_i)

    do j = 1, size(xfar,2) ! looping over theta
        do ii = 1, mapping2(j) ! looping over the number of phi values to evaluate at, for this theta
            i = mapping(ii,j) ! get the position of this phi value in the main arrays, and continue as normal

                phi = phi_vals(i) 

                ! karczewski theory
                call karczewski(diff_ampl,m,k,prop2,x3(i,j), y3(i,j),z3(i,j),r1(i,j))

                ! get the vector perpendicular to the scattering plane as viewed in the aperture system
                if(abs(dot_product(incidence2,k)) .lt. 0.999) then ! if bin is not in direct forwards
                    call cross(incidence2,k,hc)
                else ! rotate vector perp to scattering plane (y-z) into aperture system
                    my_k = (/x3(i,theta_i),y3(i,theta_i),z3(i,theta_i)/)/sqrt(x3(i,theta_i)**2 + y3(i,theta_i)**2 + z3(i,theta_i)**2)
                    call cross(incidence2,my_k,hc)
                end if

                call cross(k,m,evo2)

                ! make rotation matrix to rotate from KW plane to scattering plane
                rot4(1,1) = dot_product(hc,m)
                rot4(2,2) = dot_product(hc,m)
                rot4(1,2) = -dot_product(hc,evo2)
                rot4(2,1) = +dot_product(hc,evo2)          

                ! get rotation matrix to rotate from incidence plane (x-y) about incidence direction to scattering plane
                temp_vec3 = (/cos(phi),sin(phi),0d0/)
                call getRotationMatrix(temp_rot1,-temp_vec3(2),temp_vec3(1),0d0,1d0,0d0,0d0,0d0,0d0,-1d0)

                ! bodge (rotation was wrong way)
                temp_rot1 = transpose(temp_rot1)

                ! apply rotation matrices
                ampl_temp2 = matmul(rot4,matmul(diff_ampl,matmul(ampl,temp_rot1)))

                amplC(i,j,1,1) = ampl_temp2(1,1)
                amplC(i,j,1,2) = ampl_temp2(1,2)
                amplC(i,j,2,1) = ampl_temp2(2,1)
                amplC(i,j,2,2) = ampl_temp2(2,2)
        end do
    end do

    ! scalar fraunhofer integral, output contained in area_facs2
    call contour_integral(lambda,area_facs2,rot,rot2,v20,prop2,x3,y3,z3,mapping,mapping2,nv)

    do j = 1, size(xfar,2) ! looping over theta
        do ii = 1, mapping2(j) ! looping over the number of phi values to evaluate at, for this theta
            i = mapping(ii,j) ! get the position of this phi value in the main arrays, and continue as normal
            amplC(i,j,1,1) = amplC(i,j,1,1) * area_facs2(i,j)
            amplC(i,j,1,2) = amplC(i,j,1,2) * area_facs2(i,j)
            amplC(i,j,2,1) = amplC(i,j,2,1) * area_facs2(i,j)
            amplC(i,j,2,2) = amplC(i,j,2,2) * area_facs2(i,j)
        end do
    end do

    end subroutine
    
    subroutine diff_main(   beam_outbeam_tree,          & ! tree of outgoing beams
                            beam_outbeam_tree_counter,  & ! total number of outoing beams in beam tree
                            ampl_far_beam,            & ! far-field diffracted beam amplitude matrix
                            ext_diff_outbeam_tree,      & ! tree of outgoing external diffraction beams
                            ampl_far_ext_diff,        & ! far-field external diffraction amplitude matrix
                            job_params, &
                            geometry)
    
        ! sr diff_main is the main shell for diffraction of all beams + external diffraction at a fixed orientation

    type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
    integer(8), intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
    real(8) lambda ! wavelength
    real(8), dimension(:), allocatable :: theta_vals
    real(8), dimension(:), allocatable :: phi_vals
    real(8), dimension(:,:), allocatable :: xfar, yfar, zfar ! far-field bin positions
    logical is_multithreaded
    type(job_parameters_type), intent(in) :: job_params
    integer(8) j
    complex(8), dimension(:,:,:,:), allocatable, intent(out) :: ampl_far_beam
    complex(8), dimension(:,:,:,:), allocatable, intent(out) :: ampl_far_ext_diff
    complex(8), dimension(:,:,:,:), allocatable :: amplC ! summand
    complex(8) ampl(1:2,1:2) ! amplitude matrix to pass to diffraction sr
    real(8) perp0(1:3) ! perp field direction to pass to diffraction sr
    real(8) prop0(1:3) ! outgoing propagation direction to pass to diffraction sr
    real(8) v(1:3,1:3) ! vertices of facet to pass to diffraction sr
    complex(8), dimension(:,:), allocatable :: area_facs2
    type(outbeamtype), dimension(:), allocatable, intent(in) :: ext_diff_outbeam_tree
    real(8) start, finish, start1, finish1 ! cpu timing variables
    real(8) progressReal
    real work_done
    integer progressInt
    integer(8) num_threads
    real(8) fov
    type(outbeamtype) outbeam
    type(geometry_type), intent(in) :: geometry

    if(job_params%timing) then
        start = omp_get_wtime()
    end if

    lambda = job_params%la
    theta_vals = job_params%theta_vals
    phi_vals = job_params%phi_vals
    is_multithreaded = job_params%is_multithreaded

    call make_far_field_bins(xfar,yfar,zfar,theta_vals,phi_vals) ! get meshgrid-style far-field bins x, y, z

    ! allocate some far-field bins and init
    allocate(ampl_far_beam(1:size(xfar,1),1:size(xfar,2),1:2,1:2))
    allocate(ampl_far_ext_diff(1:size(xfar,1),1:size(xfar,2),1:2,1:2))
    allocate(area_facs2(1:size(xfar,1),1:size(xfar,2)))
    allocate(amplC(1:size(xfar,1),1:size(xfar,2),1:2,1:2))
    
    ! init
    area_facs2 = 0
    ampl_far_beam = 0
    ampl_far_ext_diff = 0
    progressInt = 0
    work_done = 0

    if(job_params%debug >= 1) then
        write(101,*)'start diff beam loop...'
    end if

    if(job_params%is_multithreaded) then
        num_threads = omp_get_max_threads()
    else
        num_threads = 1
    end if

    !$OMP PARALLEL num_threads(num_threads) PRIVATE(amplC,ampl,perp0,prop0,v,start1,finish1,fov,outbeam)
    if(job_params%timing) then
        start1 = omp_get_wtime()
    end if
    !$OMP DO
    do j = 1, beam_outbeam_tree_counter
        
        if(omp_get_thread_num() .eq. 0) then
            if(job_params%timing .and. job_params%debug >=1) then
                work_done = work_done + 1
                progressReal = work_done*100/beam_outbeam_tree_counter*omp_get_num_threads()          ! my thread percent completion
                if(int(progressReal) .gt. progressInt .and. mod(int(floor(progressReal)),10) .eq. 0) then  ! if at least 10% progress has been made
                    progressInt = int(progressReal)           ! update progress counter
                    call progress_bar(progressInt, 100)
                end if
            end if
        end if

        ampl(1:2,1:2) = beam_outbeam_tree(j)%ampl(1:2,1:2)
        perp0(1:3) = beam_outbeam_tree(j)%vk7(1:3)
        prop0(1:3) = beam_outbeam_tree(j)%prop_out(1:3)
        v(1:3,1:3) = beam_outbeam_tree(j)%verts(1:3,1:3)
        fov = beam_outbeam_tree(j)%fov
        outbeam = beam_outbeam_tree(j)
    
        ! get far-field contribution from this outbeam
        call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC,phi_vals,theta_vals,fov,job_params,outbeam,geometry)                                 

        !$OMP CRITICAL
        ampl_far_beam(:,:,:,:) = ampl_far_beam(:,:,:,:) + amplC(:,:,:,:)
        !$OMP END CRITICAL
    
    end do

    if(omp_get_thread_num() .eq. 0) then
        if(job_params%debug >= 1) then
            write(101,*)'end diff beam loop...'
            if(job_params%timing) then
                finish1 = omp_get_wtime()
                write(101,'(A,f16.8,A)')"beam diffraction took: ",finish1-start1," secs"
                start1 = omp_get_wtime()
            end if
            write(101,*)'start ext diff loop...'
        end if
    end if

    work_done = 0
    progressInt = 0
    !$OMP DO
    do j = 1, size(ext_diff_outbeam_tree,1)

        if(omp_get_thread_num() .eq. 0) then
            if(job_params%timing .and. job_params%debug >=1) then
                work_done = work_done + 1
                progressReal = work_done*100/size(ext_diff_outbeam_tree,1)*omp_get_num_threads()          ! my thread percent completion
                if(int(progressReal) .gt. progressInt .and. mod(int(floor(progressReal)),10) .eq. 0) then  ! if at least 10% progress has been made
                    progressInt = int(progressReal)           ! update progress counter
                    call progress_bar(progressInt, 100)
                end if
            end if
        end if

        ampl(1:2,1:2) = ext_diff_outbeam_tree(j)%ampl(1:2,1:2)
        perp0(1:3) = ext_diff_outbeam_tree(j)%vk7(1:3)
        prop0(1:3) = ext_diff_outbeam_tree(j)%prop_out(1:3)
        v(1:3,1:3) = ext_diff_outbeam_tree(j)%verts(1:3,1:3)
        fov = pi/2
        outbeam = ext_diff_outbeam_tree(j)
    
        call diffraction(ampl,v,prop0,perp0,xfar,yfar,zfar,lambda,amplC,phi_vals,theta_vals,fov,job_params,outbeam,geometry)

        !$OMP CRITICAL
        ampl_far_ext_diff(:,:,:,:) = ampl_far_ext_diff(:,:,:,:) + amplC(:,:,:,:)
        !$OMP END CRITICAL
            
    end do

    if(omp_get_thread_num() .eq. 0) then
        if(job_params%debug >= 1) then
            write(101,*)'end ext diff loop...'
            if(job_params%timing) then
                finish1 = omp_get_wtime()
                write(101,'(A,f16.8,A)')"external diffraction took: ",finish1-start1," secs"
            end if
        end if
    end if
    !$OMP END PARALLEL
    
    if(job_params%debug >= 1) then
        if(job_params%timing) then
            finish = omp_get_wtime()
            write(101,'(A,f16.8,A)')"diffraction took: ",finish-start," secs"
        end if
    end if

    end subroutine
    
    subroutine make_far_field_bins(x, y, z,theta_vals, phi_vals) ! optional args
    
        ! makes the far field x, y, z bins at which the diffracted field will be computed
    
    integer theta_dim, phi_dim ! theta and phi discretisation
    real(8), dimension(:,:), allocatable, intent(out) :: x, y, z ! far-field bin positions    
    real(8) r ! distance to far-field
    real(8), dimension(:), allocatable, intent(in) :: theta_vals
    real(8), dimension(:), allocatable, intent(in) :: phi_vals
    real(8), dimension(:,:), allocatable :: theta_vals_mesh, phi_vals_mesh

    ! call read_theta_vals(theta_vals)

    ! call read_phi_vals(phi_vals)

    phi_dim = size(phi_vals,1)
    theta_dim = size(theta_vals,1)

    ! allocate
    ! allocate(theta_vals(1:theta_dim))
    ! allocate(phi_vals(1:phi_dim))
    allocate(x(1:phi_dim,1:theta_dim))
    allocate(y(1:phi_dim,1:theta_dim))
    allocate(z(1:phi_dim,1:theta_dim))

    ! stop
    
    ! stop


    ! numerical fixes due to divide by 0 in contour integral
    ! do i = 1, size(phi_vals)
    !     if(abs(phi_vals(i)*180/pi) .lt. 0.000001) phi_vals(i) = phi_vals(i) + 0.00001*pi/180
    !     if(abs(phi_vals(i)*180/pi - 360.0) .lt. 0.000001) phi_vals(i) = phi_vals(i) + 0.00001*pi/180
    ! end do

    ! stop
    
    ! convert 1-d arrays to 2-d arrays
    call meshgrid_real(theta_vals, phi_vals, theta_vals_mesh, phi_vals_mesh)

    r = 1e6 ! distance to far-field
    
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