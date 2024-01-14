! beam_loop_mod.f90
! module for the beam tracing loop

module beam_loop_mod
    
    use misc_submod
    use types_mod
    use omp_lib
    
    implicit none
    
    contains
    
    subroutine print_beam_info(beam,job_params)

        ! print_beam_info
        ! prints useful debugging information about the propagation of a beam

        type(beam_type), intent(in) :: beam
        type(job_parameters_type), intent(in) :: job_params

        if(job_params%debug >= 3) then
            print*,'-----------------------------------------------'
            if(beam%is_int) then
                print'(a,i6,a)','beam ',beam%id,': is internally propagating'
            else
                print'(a,i6,a)','beam ',beam%id,': is externally propagating'
            end if
            print'(a,i6,a,i8)','beam ',beam%id,': originates from aperture: ',beam%ap
            print'(a,i6,a,i8)','beam ',beam%id,': number of facets in this beam: ',beam%nf_in
            print'(a,i6,a,i8)','beam ',beam%id,': number of facets illuminated: ',beam%nf_out
            print'(a,i6,a,f14.8)','beam ',beam%id,': scatt. cross section in: ',beam%scatt_in
            print'(a,i6,a,f14.8)','beam ',beam%id,': scatt. cross section out: ',beam%scatt_out
            print'(a,i6,a,f14.8)','beam ',beam%id,': geo cross section in: ',beam%proj_area_in
            print'(a,i6,a,f14.8)','beam ',beam%id,': geo cross section out: ',beam%proj_area_out
        end if


    end subroutine

    subroutine recursion_inc(beam_inc,geometry,job_params,beam_geometry,ext_diff_outbeam_tree)

        ! recursion_inc
        ! propagates the incident beam
        ! finds the visible facets as viewed along the initial propagation direction (-z)
        ! propagates the field from the incident beam to the visible facets
        ! initialises the array for the external diffraction
        
        type(geometry_type), intent(in) :: geometry ! geometry data structure
        type(job_parameters_type), intent(in) :: job_params ! job parameters
        type(beam_type), intent(inout) :: beam_inc ! current beam to be traced
        type(geometry_type), intent(in) :: beam_geometry ! geometry of the incident beam
        type(outbeamtype), dimension(:), allocatable, intent(out) :: ext_diff_outbeam_tree ! external diffraction tree
        
        integer(8) i, j ! some counters
        integer(8) ai ! an aperture id
        real(8) rot2(1:2,1:2) ! a 2x2 rotation matrix
        logical, dimension(:), allocatable :: in_beam ! whether each facet was within the beam
        real(8), dimension(:), allocatable :: dist_beam ! the distance from the beam to each facet
        integer(8), dimension(:), allocatable :: id_beam ! the facet id of the beam face which illuminated each facet
        logical, dimension(:), allocatable :: is_shad ! whether each facet was in the shadow of another
        logical, dimension(:), allocatable :: is_vis ! whether each facet was visible as viewed along the beam propagation direction
        real(8) prop(1:3) ! incoming propagation direction in unrotated system
        real(8) prop_int(1:3) ! reflected/refracted propagation vector
        real(8) prop_ext(1:3) ! reflected/refracted propagation vector
        integer(8) num_ill_facets ! a counter
        integer(8) n ! a counter
        real(8) an(1:3) ! an aperture normal
        real(8) ran(1:3) ! a rotated aperture normal
        real(8) vk7(1:3) ! a reflected/refracted e-perp vector
        real(8) e_perp(1:3) ! an incident e-perp vector
        integer(8) bfi ! an beam face id
        complex(8) ampl(1:2,1:2) ! an amplitude matrix
        complex(8) refl_ampl(1:2,1:2) ! a reflected amplitude matrix
        complex(8) trans_ampl(1:2,1:2) ! a transmitted amplitude matrix
        real(8) waveno ! wave number
        real(8) rbi_int ! internal real part refractive index
        real(8) ibi_int ! internal imaginary part refractive index
        complex(8) m_int ! internal complex refractive index
        real(8) dist ! distance travelled
        logical is_tir ! whether or not there is total internal reflection
        real(8) theta_i ! incident angle with the aperture
        real(8) theta_i_facet ! incident angle with the facet
        real(8) theta_t ! transmitted angle
        complex(8) fr(1:2,1:2) ! fresnel reflection matrix
        complex(8) ft(1:2,1:2) ! fresnel reflection matrix
        real(8) norm(1:3) ! a facet normal
        real(8) int_intensity ! internal intensitiy
        real(8) ext_intensity ! external intensitiy
        real(8) proj_area ! a projected area
        
        waveno = 2*pi/job_params%la ! wavenumber
        rbi_int = job_params%rbi ! real part refractive index
        ibi_int = job_params%ibi ! imaginary part refractive index
        m_int = cmplx(rbi_int,ibi_int,kind=8) ! refractive index
        
        ! no need to rotate as propagation is already along the -z direction
        
        ! find the facets illuminated by the beam
        call find_vis_incidence(in_beam, is_vis, dist_beam, id_beam, is_shad, beam_inc, geometry, beam_geometry)
        
        ! get the total number of facets illuminated by this beam
        num_ill_facets = 0
        do i = 1, geometry%na
            do j = 1, geometry%nf
                if(in_beam(j) .and. geometry%f(j)%ap == i) then
                    num_ill_facets = num_ill_facets + 1
                end if
            end do
        end do
        
        ! allocate space in the beam structure to hold the illuminated field
        beam_inc%nf_out = num_ill_facets
        allocate(beam_inc%field_out(1:beam_inc%nf_out))
        allocate(ext_diff_outbeam_tree(1:beam_inc%nf_out)) ! external diffraction
        
        ! add info about the illuminated field to the beam structure
        n = 0 ! init a counter to count the number of faces added to the beam struct
        do i = 1, geometry%nf ! for each face
            if(in_beam(i)) then ! if this face was illuminated by the beam
                n = n + 1 ! update counter
                ai = geometry%f(i)%ap ! get the aperture id (in unrotated system)
                an(:) = geometry%ap(ai)%n(:) ! get the aperture normal (in unrotated systen)
                ran(:) = geometry%ap(ai)%n(:) ! get the aperture normal (in rotated system)
                bfi = id_beam(i) ! get the beam facet id which illuminated this one
                ampl(:,:) = beam_inc%field_in(bfi)%ampl(:,:) ! get the amplitude matrix of the illuminating facet
                e_perp(:) = beam_inc%field_in(bfi)%e_perp(:) ! get the incident e-perp vector
                dist = dist_beam(i) ! distance travelled from illuminating to illuminated facet
                prop(:) = beam_inc%prop(:) ! incoming propagation direction in unrotated system
                
                ! get the reflected propagation vector in unrotated system
                prop_ext(1) =  0 + 2*geometry%ap(ai)%n(3)*geometry%ap(ai)%n(1) ! reflected propagation vector (in rotated system)
                prop_ext(2) =  0 + 2*geometry%ap(ai)%n(3)*geometry%ap(ai)%n(2) ! reflected propagation vector (in rotated system)
                prop_ext(3) = -1 + 2*geometry%ap(ai)%n(3)*geometry%ap(ai)%n(3) ! reflected propagation vector (in rotated system)
                
                ! get vk7, the new e-perp vector (normal x reflected prop vector)
                call cross(an,prop,vk7,.true.)
                
                ! get the rotation matrix to rotate about prop. vector into new scattering plane
                call get_rot_matrix(rot2,vk7,e_perp,prop)
                
                ! rotate the amplitude matrix into the new scattering plane
                ampl(:,:) = matmul(rot2(:,:),ampl(:,:))
                
                ! apply distance phase factor
                ampl(:,:) = ampl(:,:) * exp2cmplx(waveno*dist)

                ! compute fresnel matrices
                fr(:,:) = 0d0 ! init fresnel reflection matrix
                ft(:,:) = 0d0 ! init fresnel transmission matrix
                if(is_shad(i)) then ! if facet is in shadow (needs careful consideration)
                    theta_i = acos(an(3)) ! using incident angle as that of the aperture if in shadow
                else
                    theta_i = acos(geometry%n(geometry%f(i)%ni,3))
                end if
                is_tir = .false.
                theta_t = asin(sin(theta_i)/real(m_int))
                fr(2,2) = (cos(theta_i) - m_int*cos(theta_t))/(cos(theta_i) + m_int*cos(theta_t))
                ft(2,2) = (2*cos(theta_i))/(cos(theta_i) + m_int*cos(theta_t))
                fr(1,1) = (m_int*cos(theta_i) - cos(theta_t))/(cos(theta_t) + m_int*cos(theta_i))
                ft(1,1) = (2*cos(theta_i))/(cos(theta_t) + m_int*cos(theta_i))
                
                ! apply fresnel matrices to amplitude matrix
                refl_ampl(:,:) = matmul(fr(:,:),ampl(:,:))
                trans_ampl(:,:) = matmul(ft(:,:),ampl(:,:))
                
                ! compute internal intensity
                int_intensity = real(0.5*(  trans_ampl(1,1)*conjg(trans_ampl(1,1)) + &
                trans_ampl(1,2)*conjg(trans_ampl(1,2)) + &
                trans_ampl(2,1)*conjg(trans_ampl(2,1)) + &
                trans_ampl(2,2)*conjg(trans_ampl(2,2))))

                ! compute external intensity
                ext_intensity = real(0.5*(  refl_ampl(1,1)*conjg(refl_ampl(1,1)) + &
                refl_ampl(1,2)*conjg(refl_ampl(1,2)) + &
                refl_ampl(2,1)*conjg(refl_ampl(2,1)) + &
                refl_ampl(2,2)*conjg(refl_ampl(2,2))))

                ! possibly remove reflection from shadow facets here
                ! if(is_shad(i)) refl_ampl(:,:) = 0d0
                
                ! get projected facet area
                theta_i_facet = acos(geometry%n(geometry%f(i)%ni,3))
                proj_area = geometry%f(i)%area * cos(theta_i_facet)

                ! compute the transmitted propagation vector (from the whole aperture)
                norm(:) = -an(:) ! get aperture normal
                call get_trans_prop(theta_i,theta_t,prop,norm,prop_int) ! get the transmitted propagation vector
                        
                ! save stuff to the beam structure
                beam_inc%field_out(n)%ampl_ext(:,:) = refl_ampl(:,:) ! save the transmitted amplitude matrix
                beam_inc%field_out(n)%ampl_int(:,:) = trans_ampl(:,:) ! save the reflected amplitude matrix
                beam_inc%field_out(n)%e_perp(:) = vk7(:) ! save the perpendicular field vector
                beam_inc%field_out(n)%fi = i ! save the face id
                beam_inc%field_out(n)%ap = ai ! save the aperture
                beam_inc%field_out(n)%prop_int(:) = prop_int(:) ! save the internally reflected propagation vector
                beam_inc%field_out(n)%prop_ext(:) = prop_ext(:) ! save the externally transmitited propagation vector
                beam_inc%field_out(n)%is_tir = is_tir ! save whether or not this field was total internal reflection
                beam_inc%field_out(n)%scatt_int = int_intensity * proj_area * (rbi_int * (cos(theta_t)/cos(theta_i_facet))) ! save the scattering cross section contribution
                beam_inc%field_out(n)%scatt_ext = ext_intensity * proj_area ! save the scattering cross section contribution
                beam_inc%field_out(n)%proj_area = proj_area ! save the geometric cross section
                
                ! save stuff to external diffraction tree
                ext_diff_outbeam_tree(n)%ampl(:,:) = ampl(:,:)
                ext_diff_outbeam_tree(n)%vk7(:) = vk7(:)
                ext_diff_outbeam_tree(n)%prop_out(:) = prop(:)
                ext_diff_outbeam_tree(n)%prop_in(:) = prop(:)
                ext_diff_outbeam_tree(n)%fout = i ! face id from which the beam was emitted
                ext_diff_outbeam_tree(n)%verts(:,:) = transpose(geometry%v(geometry%f(i)%vi(:),:)) ! verts needs readjusting for arbitrary number of vertices
            end if
        end do
        
    end subroutine
    
    subroutine recursion_int(beam,geometry,job_params)

        ! recursion_int
        ! propagates an internal beam
        ! finds the visible facets as viewed along the initial propagation direction (-z)
        ! propagates the field from the incident beam to the visible facets      
        
        type(geometry_type), intent(in) :: geometry ! geometry data structure
        type(job_parameters_type), intent(in) :: job_params ! job parameters
        type(beam_type), intent(inout) :: beam ! current beam to be traced

        type(geometry_type) :: rot_geometry ! the geometry structure after rotating
        integer(8) i, j ! some counters
        integer(8) ai ! an aperture id
        real(8) rot(1:3,1:3) ! a 3x3 rotation matrix
        real(8) rot2(1:2,1:2) ! a 2x2 rotation matrix
        logical, dimension(:), allocatable :: is_shad ! whether each facet was in the shadow of another
        logical, dimension(:), allocatable :: in_beam ! whether each facet was within the beam
        integer(8), dimension(:), allocatable :: id_beam ! the facet id of the beam face which illuminated each facet
        real(8), dimension(:), allocatable :: dist_beam ! the distance from the beam to each facet
        real(8) prop(1:3) ! incoming propagation direction in unrotated system
        real(8) prop_int(1:3) ! reflected/refracted propagation vector
        real(8) prop_ext(1:3) ! reflected/refracted propagation vector
        integer(8) num_ill_facets    
        integer(8) n ! a counter
        real(8) an(1:3) ! an aperture normal
        real(8) ran(1:3) ! a rotated aperture normal
        real(8) vk7(1:3) ! a reflected/refracted e-perp vector
        real(8) e_perp(1:3) ! an incident e-perp vector
        integer(8) bfi ! an beam face id
        complex(8) ampl(1:2,1:2) ! an amplitude matrix
        complex(8) refl_ampl(1:2,1:2) ! a reflected amplitude matrix
        complex(8) trans_ampl(1:2,1:2) ! a transmitted amplitude matrix
        real(8) waveno ! wave number
        real(8) rbi_int ! internal real part refractive index
        real(8) ibi_int ! internal imaginary part refractive index
        complex(8) m_int ! internal complex refractive index
        real(8) dist ! distance travelled
        logical is_tir ! whether or not there is total internal reflection
        real(8) theta_i ! incident angle
        real(8) theta_i_facet ! incident angle with the facet
        real(8) theta_t ! transmitted angle
        complex(8) fr(1:2,1:2) ! fresnel reflection matrix
        complex(8) ft(1:2,1:2) ! fresnel reflection matrix
        real(8) abs_intensity ! absorbed intensity
        real(8) in_intensity ! absorbed intensity
        real(8) out_intensity ! intensity after absorption
        real(8) int_intensity ! internal intensitiy
        real(8) ext_intensity ! external intensitiy        
        real(8) norm(1:3) ! a facet normal
        real(8) proj_area ! a projected area
        
        waveno = 2*pi/job_params%la ! wavenumber
        rbi_int = job_params%rbi ! real part refractive index
        ibi_int = job_params%ibi ! imaginary part refractive index
        m_int = cmplx(rbi_int,ibi_int,kind=8) ! refractive index
        
        ! rotate geometry so that the propagation is along the z-axis, save the rotation matrix
        call rotate(beam,geometry,rot_geometry,rot)
        
        ! find the facets illuminated by the beam
        call find_vis_int(in_beam, dist_beam, id_beam, is_shad, beam, rot_geometry)
        
        ! get the total number of facets illuminated by this beam
        num_ill_facets = 0
        do i = 1, rot_geometry%na
            do j = 1, rot_geometry%nf
                if(in_beam(j) .and. rot_geometry%f(j)%ap == i) then
                    num_ill_facets = num_ill_facets + 1
                end if
            end do
        end do
        
        ! allocate space in the beam structure to hold the illuminated field
        beam%nf_out = num_ill_facets
        allocate(beam%field_out(1:beam%nf_out))
        
        ! add info about the illuminated field to the beam structure
        n = 0 ! init a counter to count the number of faces added to the beam struct
        do i = 1, rot_geometry%nf ! for each face
            if(in_beam(i)) then ! if this face was illuminated by the beam
                n = n + 1 ! update counter
                ai = geometry%f(i)%ap ! get the aperture id (in unrotated system)
                an(:) = geometry%ap(ai)%n(:) ! get the aperture normal (in unrotated systen)
                ran(:) = rot_geometry%ap(ai)%n(:) ! get the aperture normal (in rotated system)
                bfi = id_beam(i) ! get the beam facet id which illuminated this one
                ampl(:,:) = beam%field_in(bfi)%ampl(:,:) ! get the amplitude matrix of the illuminating facet
                e_perp(:) = beam%field_in(bfi)%e_perp(:) ! get the incident e-perp vector
                dist = dist_beam(i) ! distance travelled from illuminating to illuminated facet
                prop(:) = beam%prop(:) ! incoming propagation direction in unrotated system
                
                ! get the reflected propagation vector in unrotated system
                prop_int(1) =  0 + 2*rot_geometry%ap(ai)%n(3)*rot_geometry%ap(ai)%n(1) ! reflected propagation vector (in rotated system)
                prop_int(2) =  0 + 2*rot_geometry%ap(ai)%n(3)*rot_geometry%ap(ai)%n(2) ! reflected propagation vector (in rotated system)
                prop_int(3) = -1 + 2*rot_geometry%ap(ai)%n(3)*rot_geometry%ap(ai)%n(3) ! reflected propagation vector (in rotated system)
                prop_int(:) = matmul(transpose(rot),prop_int(:)) ! rotate reflected/refracted propagation vector back to original coordinate system
                
                ! get vk7, the new e-perp vector (normal x reflected prop vector)
                call cross(-an,prop,vk7,.true.)
                
                ! get the rotation matrix to rotate about prop. vector into new scattering plane
                call get_rot_matrix(rot2,vk7,e_perp,prop)
                
                ! rotate the amplitude matrix into the new scattering plane
                ampl(:,:) = matmul(rot2(:,:),ampl(:,:))
                
                ! apply distance phase factor
                ampl(:,:) = ampl(:,:) * exp2cmplx(waveno*rbi_int*dist)
                
                ! compute incident intensity (before applying the absorption factor)
                in_intensity = real(0.5*(  ampl(1,1)*conjg(ampl(1,1)) + &
                ampl(1,2)*conjg(ampl(1,2)) + &
                ampl(2,1)*conjg(ampl(2,1)) + &
                ampl(2,2)*conjg(ampl(2,2))))

                ! compute absorbed intensity (before applying the absorption factor)
                abs_intensity = in_intensity * (1d0-exp(-2*waveno*ibi_int*sqrt(dist))**2)
                
                ! apply absorption factor
                ampl(:,:) = ampl(:,:) * exp(-2*waveno*ibi_int*sqrt(dist))
                
                ! compute remaining intensity after absorption
                out_intensity = real(0.5*(  ampl(1,1)*conjg(ampl(1,1)) + &
                ampl(1,2)*conjg(ampl(1,2)) + &
                ampl(2,1)*conjg(ampl(2,1)) + &
                ampl(2,2)*conjg(ampl(2,2))))
                
                ! compute fresnel matrices
                fr(:,:) = 0d0 ! init fresnel reflection matrix
                ft(:,:) = 0d0 ! init fresnel transmission matrix
                theta_i = acos(-ran(3)) ! using incident angle as that of the aperture
                if(beam%is_int .and. theta_i >= asin(1/rbi_int)) then ! if tir
                    is_tir = .true. 
                    fr(2,2) = -1d0
                    ft(2,2) = 0d0
                    fr(1,1) = -1d0
                    ft(1,1) = 0d0
                else ! if not tir
                    is_tir = .false.
                    theta_t = asin(sin(theta_i)*rbi_int)
                    fr(2,2) = (m_int*cos(theta_i) - cos(theta_t))/(m_int*cos(theta_i) + cos(theta_t))
                    ft(2,2) = (2*m_int*cos(theta_i))/(m_int*cos(theta_i) + cos(theta_t))
                    fr(1,1) = (cos(theta_i) - m_int*cos(theta_t))/(m_int*cos(theta_t) + cos(theta_i))
                    ft(1,1) = (2*m_int*cos(theta_i)) / (m_int*cos(theta_t) + cos(theta_i))
                end if
                
                ! apply fresnel matrices to amplitude matrix
                refl_ampl(:,:) = matmul(fr(:,:),ampl(:,:))
                trans_ampl(:,:) = matmul(ft(:,:),ampl(:,:))
                
                ! compute internal intensity
                int_intensity = real(0.5*(  refl_ampl(1,1)*conjg(refl_ampl(1,1)) + &
                refl_ampl(1,2)*conjg(refl_ampl(1,2)) + &
                refl_ampl(2,1)*conjg(refl_ampl(2,1)) + &
                refl_ampl(2,2)*conjg(refl_ampl(2,2))))

                ! compute external intensity
                ext_intensity = real(0.5*(  trans_ampl(1,1)*conjg(trans_ampl(1,1)) + &
                trans_ampl(1,2)*conjg(trans_ampl(1,2)) + &
                trans_ampl(2,1)*conjg(trans_ampl(2,1)) + &
                trans_ampl(2,2)*conjg(trans_ampl(2,2))))

                ! possibly remove transmission from shadow facets here
                if(is_shad(i)) trans_ampl(:,:) = 0d0
                
                ! get projected facet area
                theta_i_facet = acos(-(rot_geometry%n(rot_geometry%f(i)%ni,3)))
                proj_area = geometry%f(i)%area * cos(theta_i_facet)

                ! add to the beam absorption cross section
                beam%abs = beam%abs + (abs_intensity * proj_area * rbi_int) ! sum the absorbed energy

                ! add to the beam output cross section
                beam%scatt_out = beam%scatt_out + (out_intensity * proj_area * rbi_int)

                ! compute the transmitted propagation vector
                if(.not. is_tir) then ! if no internal reflection
                    norm(:) = geometry%n(geometry%f(i)%ni,:)
                    call get_trans_prop(theta_i,theta_t,prop,norm,prop_ext) ! get the transmitted propagation vector
                end if
                
                ! save stuff to the beam structure
                beam%field_out(n)%ampl_ext(:,:) = trans_ampl(:,:) ! save the transmitted amplitude matrix
                beam%field_out(n)%ampl_int(:,:) = refl_ampl(:,:) ! save the reflected amplitude matrix
                beam%field_out(n)%e_perp(:) = vk7(:) ! save the perpendicular field vector
                beam%field_out(n)%fi = i ! save the face id
                beam%field_out(n)%ap = ai ! save the aperture
                beam%field_out(n)%prop_int(:) = prop_int(:) ! save the internally reflected propagation vector
                beam%field_out(n)%prop_ext(:) = prop_ext(:) ! save the externally transmitited propagation vector
                beam%field_out(n)%is_tir = is_tir ! save whether or not this field was total internal reflection
                beam%field_out(n)%scatt_int = int_intensity * proj_area ! save the scattering cross section contribution
                beam%field_out(n)%scatt_ext = ext_intensity * proj_area * (cos(theta_t)/cos(theta_i_facet)) ! save the scattering cross section contribution
                beam%field_out(n)%proj_area = proj_area ! save the geometric cross section

                beam%proj_area_out = beam%proj_area_out + proj_area ! sum the total projected area for this beam
            end if
        end do

    end subroutine
    
    subroutine get_trans_prop(theta_i,theta_t,prop_in,norm,prop_out)

        ! get_trans_prop
        ! gets a transmitted propagation vector from an incident angle,
        !   transmitted angle, incident propagation vector, and surface normal
        
        real(8), intent(in) :: theta_i
        real(8), intent(in) :: theta_t
        real(8), intent(in) :: prop_in(1:3)
        real(8), intent(out) :: prop_out(1:3)
        real(8), intent(in) :: norm(1:3)
        
        real(8) a, b, alpha, nf
        real(8), dimension(1:3) :: a_vec, b_vec
        
        alpha = pi - theta_t
        a = sin(theta_t-theta_i)/sin(theta_i)
        b = sin(alpha)/sin(theta_i)
        b_vec(:) = prop_in(:)
        a_vec(:) = norm(:)
        prop_out(:) = b*b_vec(:) - a*a_vec(:)
        nf = sqrt(prop_out(1)**2 + prop_out(2)**2 + prop_out(3)**2)
        prop_out(:) = prop_out(:) / nf
        
    end subroutine
    
    subroutine add_to_beam_tree_internal(beam_tree,beam,num_beams,geometry,job_params)
        
        ! add_to_beam_tree_internal
        ! for an internally propagating beam, adds the new beams to the beam tree
        ! after an internal beam has been propagated in sr recursion_int, it passed here to add
        !   the new beams to the beam tree, which may be then later traced
        
        type(beam_type), dimension(:), allocatable, intent(inout) :: beam_tree
        type(beam_type), intent(in) :: beam ! current beam to be traced
        integer(8), intent(inout) :: num_beams
        type(geometry_type), intent(in) :: geometry
        type(job_parameters_type), intent(in) :: job_params
        
        integer(8) i
        integer(8) num_new_beams ! number of beams to add to the beam tree
        logical, dimension(:), allocatable :: suff_ill ! whether each aperture was sufficiently illuminated
        integer(8), dimension(:), allocatable :: num_ill ! number of illuminated facets in each aperture
        integer(8), dimension(:), allocatable :: mapping ! a mapping
        integer(8) ai ! aperture id
        integer(8) fi ! aperture id
        integer(8) index ! index of a beam in the beam tree
        integer(8) nf ! a counter for the number of facets in a beam tree entry
        real(8) proj_area ! the face area projected along the new propagation direction

        ! determine which apertures were sufficiently illuminated to create new beams
        call get_suff_illuminated(geometry,beam,job_params,suff_ill,num_ill)
        
        allocate(mapping(1:geometry%na)) ! maps each aperture to a position in the beam tree (if suff. illuminated)
        mapping = 0
        num_new_beams = 0 ! init

        ! initialise the new beams to be added to the beam tree
        do i = 1, geometry%na ! for each aperture in the geometry
            if(suff_ill(i)) then ! if it was sufficiently illuminated by this beam
                num_new_beams = num_new_beams + 1 ! update counter
                num_beams = num_beams + 1 ! add an internally reflected beam to the beam tree
                if(num_beams > size(beam_tree,1)) then
                    print*,'not enough space in beam tree'
                    stop
                end if
                mapping(i) = num_beams ! save the position of this beam in the beam tree
                allocate(beam_tree(num_beams)%field_in(1:num_ill(i))) ! allocate space
                ! add stuff to the beam tree
                beam_tree(num_beams)%ap = i ! aperture id
                beam_tree(num_beams)%nf_in = 0 ! init
                beam_tree(num_beams)%is_int = .true. ! is internally propagating
                beam_tree(num_beams)%scatt_in = 0 ! init
                beam_tree(num_beams)%scatt_out = 0 ! init
                beam_tree(num_beams)%abs = 0 ! init
                beam_tree(num_beams)%proj_area_in = 0 ! init
                beam_tree(num_beams)%id = num_beams ! save position in beam tree
                beam_tree(num_beams)%rec = beam%rec + 1 ! save recursion number
            end if
        end do
        
        ! add information about the field of the new beams to the beam tree
        do i = 1, beam%nf_out ! for each illuminated facet
            fi = beam%field_out(i)%fi ! get the facet id
            ai = beam%field_out(i)%ap ! get the aperture id
            index = mapping(ai) ! get the position of new beam in the beam tree
            nf = beam_tree(index)%nf_in ! get (current) number of facets in this beam
            nf = nf + 1 ! update the number of facets
            proj_area = geometry%f(fi)%area * dot_product(geometry%n(geometry%f(fi)%ni,:),-beam%field_out(i)%prop_int(:))

            beam_tree(index)%nf_in = nf ! save updated number of facets
            beam_tree(index)%field_in(nf)%ampl(:,:) = beam%field_out(i)%ampl_int(:,:) ! the internally reflected field becomes the new beam field
            beam_tree(index)%field_in(nf)%e_perp(:) = beam%field_out(i)%e_perp(:)
            beam_tree(index)%field_in(nf)%fi = fi ! save the facet id
            beam_tree(index)%prop(:) = beam%field_out(i)%prop_int(:) ! save internally reflected propagation vector (overwrites each time)
            beam_tree(index)%scatt_in = beam_tree(index)%scatt_in + beam%field_out(i)%scatt_int ! sum the scattering cross section of this beam
            beam_tree(index)%proj_area_in = beam_tree(index)%proj_area_in + proj_area ! sum the projected area along the new propagation direction
        end do
        
        if(job_params%debug >= 3) print'(a,i6,a,i8,a,i8,a)','beam ',beam%id,': added ',num_new_beams,' beams to beam tree -----> ',num_beams,' total beams'
        
    end subroutine
    
    subroutine add_to_beam_tree_external(beam_tree,beam,num_beams,geometry,job_params)
        
        ! add_to_beam_tree_external
        ! for an externally propagating beam, adds the new beams to the beam tree
        ! after an external beam has been propagated in sr recursion_int, it passed here to add
        !   the new beams to the beam tree, which may be then later traced

        type(beam_type), dimension(:), allocatable, intent(inout) :: beam_tree
        type(beam_type), intent(in) :: beam ! current beam to be traced
        integer(8), intent(inout) :: num_beams
        type(geometry_type), intent(in) :: geometry
        type(job_parameters_type), intent(in) :: job_params
        
        integer(8) i
        integer(8) num_new_beams ! number of beams to add to the beam tree
        logical, dimension(:), allocatable :: suff_ill ! whether each aperture was sufficiently illuminated
        integer(8), dimension(:), allocatable :: num_ill ! number of illuminated facets in each aperture
        integer(8), dimension(:), allocatable :: mapping ! a mapping
        integer(8) ai ! aperture id
        integer(8) fi ! aperture id
        integer(8) index ! index of a beam in the beam tree
        integer(8) nf ! a counter for the number of facets in a beam tree entry
        real(8) proj_area ! the face area projected along the new propagation direction
        
        ! determine which apertures were sufficiently illuminated to create new beams
        call get_suff_illuminated(geometry,beam,job_params,suff_ill,num_ill)
        
        allocate(mapping(1:geometry%na)) ! maps each aperture to a position in the beam tree (if suff. illuminated)
        mapping = 0
        num_new_beams = 0 ! init

        ! initialise the new beams to be added to the beam tree
        do i = 1, geometry%na ! for each aperture in the geometry
            if(suff_ill(i)) then ! if it was sufficiently illuminated by this beam
                num_new_beams = num_new_beams + 1 ! update counter
                num_beams = num_beams + 1 ! add an internally reflected beam to the beam tree
                if(num_beams > size(beam_tree,1)) then
                    print*,'not enough space in beam tree'
                    stop
                end if
                mapping(i) = num_beams ! save the position of this beam in the beam tree
                allocate(beam_tree(num_beams)%field_in(1:num_ill(i))) ! allocate space
                ! add stuff to the beam tree
                beam_tree(num_beams)%ap = i ! aperture id
                beam_tree(num_beams)%nf_in = 0 ! init
                beam_tree(num_beams)%is_int = .true. ! is internally propagating
                beam_tree(num_beams)%scatt_in = 0 ! init
                beam_tree(num_beams)%scatt_out = 0 ! init
                beam_tree(num_beams)%id = num_beams ! save position in beam tree
                beam_tree(num_beams)%rec = beam%rec + 1 ! save recursion number
            end if
        end do
        
        ! add information about the field of the new beams to the beam tree
        do i = 1, beam%nf_out ! for each illuminated facet
            fi = beam%field_out(i)%fi ! get the facet id
            ai = beam%field_out(i)%ap ! get the aperture id
            index = mapping(ai) ! get the position of new beam in the beam tree
            nf = beam_tree(index)%nf_in ! get (current) number of facets in this beam
            nf = nf + 1 ! update the number of facets
            proj_area = geometry%f(fi)%area * dot_product(geometry%n(geometry%f(fi)%ni,:),-beam%field_out(i)%prop_int(:))

            beam_tree(index)%nf_in = nf ! save updated number of facets
            beam_tree(index)%field_in(nf)%ampl(:,:) = beam%field_out(i)%ampl_int(:,:) ! the internally reflected field becomes the new beam field
            beam_tree(index)%field_in(nf)%e_perp(:) = beam%field_out(i)%e_perp(:)
            beam_tree(index)%field_in(nf)%fi = fi ! save the facet id
            beam_tree(index)%prop(:) = beam%field_out(i)%prop_int(:) ! save internally reflected propagation vector (overwrites each time but they are all the same)
            beam_tree(index)%scatt_in = beam_tree(index)%scatt_in + beam%field_out(i)%scatt_int ! sum the scattering cross section of this beam
            beam_tree(index)%proj_area_in = beam_tree(index)%proj_area_in + proj_area ! sum the projected area along the new propagation direction
        end do
        
        if(job_params%debug >= 3) print'(a,i6,a,i8,a,i8,a)','beam ',beam%id,': added ',num_new_beams,' beams to beam tree -----> ',num_beams,' total beams'

    end subroutine
    
    subroutine get_suff_illuminated(geometry,beam,job_params,suff_ill,num_ill)
        
        ! get_suff_illuminated
        ! this subroutine examines the facets illuminated by a beam
        ! if the total area is greater than a threshold amount, the aperture is sufficiently illuminated
        ! sufficiently illuminated apertures create new beams that are propagated
        
        type(geometry_type), intent(in) :: geometry
        type(beam_type), intent(in) :: beam
        type(job_parameters_type), intent(in) :: job_params
        
        integer(8) i
        real(8), dimension(:), allocatable :: ap_areas ! the areas of each aperture that were illuminated
        logical, dimension(:), allocatable, intent(out) :: suff_ill ! whether each aperture was sufficiently illuminated
        integer(8), dimension(:), allocatable, intent(out) :: num_ill ! number of illuminated facets in each aperture
        integer(8) ai ! aperture id
        integer(8) fi ! face id
        
        allocate(ap_areas(1:geometry%na)) ! allocate
        allocate(suff_ill(1:geometry%na)) ! allocate
        allocate(num_ill(1:geometry%na)) ! allocate
        ap_areas = 0d0 ! init
        suff_ill = .false. ! init
        num_ill = 0 ! init
        
        do i = 1, beam%nf_out ! for each facet illuminated by the beam
            ai = beam%field_out(i)%ap ! get the aperture id of the illuminated face
            fi = beam%field_out(i)%fi ! get the facet id of the illuminated face
            ap_areas(ai) = ap_areas(ai) + geometry%f(fi)%area ! add the area to the total for this aperture
            num_ill(ai) = num_ill(ai) + 1 ! update the total number of illuminated facets for this aperture
        end do
        
        do i = 1, geometry%na ! for each aperture
            if(ap_areas(i) > job_params%threshold) suff_ill(i) = .true. ! if area large enough, set suff. illuminated
        end do
        
    end subroutine
    
    subroutine rotate(beam,geometry,rot_geometry,rot)

        ! rotate
        ! rotates a geometry
        ! the geometry is rotated so the the propagation direction of a beam is
        !   along the -z direction
        
        type(geometry_type), intent(in) :: geometry
        type(geometry_type), intent(out) :: rot_geometry
        real(8), intent(out) :: rot(1:3,1:3)
        type(beam_type), intent(in) :: beam ! current beam to be traced
        
        real(8) theta_1, theta_2
        real(8) rot1(1:3,1:3), rot2(1:3,1:3)
        real(8) temp_vector(1:3)
        integer(8) i
        
        rot1 = 0 ! initialise
        rot2 = 0 ! initialise
        rot = 0 ! initialise
        
        theta_1 = 2*pi - atan2(beam%prop(2),beam%prop(1)) ! angle to rotate about z axis into x-z plane in +ive x direction
        rot1(1,1) = cos(theta_1)
        rot1(1,2) = -sin(theta_1)
        rot1(2,1) = sin(theta_1)
        rot1(2,2) = cos(theta_1)
        rot1(3,3) = 1
        temp_vector = matmul(rot1,beam%prop(1:3)) ! to hold the propagation vector after rotating into x-z plane before we rotate onto z-axis
        call normalise_vec(temp_vector)
        theta_2 = acos(-temp_vector(3))
        rot2(1,1) = cos(theta_2)
        rot2(1,3) = sin(theta_2)
        rot2(3,1) = -sin(theta_2)
        rot2(3,3) = cos(theta_2)
        rot2(2,2) = 1
        rot = matmul(rot2,rot1)
        
        rot_geometry = geometry ! init
        
        do i = 1, geometry%nv
            rot_geometry%v(i,:) = matmul(rot,geometry%v(i,:)) ! rotate vertices
        end do
        do i = 1, geometry%nn
            rot_geometry%n(i,:) = matmul(rot,geometry%n(i,:)) ! rotate normals
        end do
        do i = 1, geometry%nf
            rot_geometry%f(i)%mid(:) = matmul(rot,geometry%f(i)%mid(:)) ! rotate midpoints
        end do
        do i = 1, geometry%na
            rot_geometry%ap(i)%mid(:) = matmul(rot,geometry%ap(i)%mid(:)) ! rotate aperture midpoints
            rot_geometry%ap(i)%n(:) = matmul(rot,geometry%ap(i)%n(:)) ! rotate aperture normals
        end do
        
    end subroutine
    
    subroutine energy_checks(   beam_outbeam_tree, &
        beam_outbeam_tree_counter, &
        output_parameters, &
        ext_diff_outbeam_tree, &
        job_params, &
        geometry)
        
        ! energy_checks
        ! performs some energy checks after the beam loop has finished

        type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
        integer(8), intent(in) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
        type(output_parameters_type), intent(inout) :: output_parameters
        real(8) :: energy_in
        real(8) :: energy_out_beam
        real(8) :: energy_out_ext_diff
        type(outbeamtype), dimension(:), allocatable, intent(in) :: ext_diff_outbeam_tree
        type(job_parameters_type), intent(in) ::  job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details
        type(geometry_type), intent(in) :: geometry
        
        integer(8) i
        integer(8) face_id
        real(8) prop(1:3)
        real(8) normal(1:3)
        real(8) cos_theta
        real(8) out_intensity
        complex(8) ampl(1:2,1:2)
        real(8) area
        
        if(job_params%debug >= 1) then
            print*,'========== start sr energy_checks'
        end if
        
        energy_in = output_parameters%geo_cross_sec
        energy_out_beam = 0
        
        do i = 1, beam_outbeam_tree_counter
            prop = beam_outbeam_tree(i)%prop_out
            face_id = beam_outbeam_tree(i)%fout
            ampl = beam_outbeam_tree(i)%ampl
            normal = geometry%n(geometry%f(face_id)%ni,:)
            area = geometry%f(face_id)%area
            cos_theta = dot_product(prop,normal)
            out_intensity = real(0.5*(   ampl(1,1)*conjg(ampl(1,1)) + &
            ampl(1,2)*conjg(ampl(1,2)) + &
            ampl(2,1)*conjg(ampl(2,1)) + &
            ampl(2,2)*conjg(ampl(2,2))))
            
            if (isnan(out_intensity)) then
                print*,'oh dear - nan ibeam = ',i
                stop
                beam_outbeam_tree(i)%ampl = 0
            else
                energy_out_beam = energy_out_beam + out_intensity*area*cos_theta
            end if
            
        end do
        
        energy_out_ext_diff = 0
        
        do i = 1, size(ext_diff_outbeam_tree,1)
            prop = ext_diff_outbeam_tree(i)%prop_out
            face_id = ext_diff_outbeam_tree(i)%fout
            ampl = ext_diff_outbeam_tree(i)%ampl
            normal = geometry%n(geometry%f(face_id)%ni,:)
            area = geometry%f(face_id)%area
            cos_theta = -dot_product(prop,normal)
            out_intensity = real(0.5*(   ampl(1,1)*conjg(ampl(1,1)) + &
            ampl(1,2)*conjg(ampl(1,2)) + &
            ampl(2,1)*conjg(ampl(2,1)) + &
            ampl(2,2)*conjg(ampl(2,2))))
            energy_out_ext_diff = energy_out_ext_diff + out_intensity*area*cos_theta
        end do
        
        output_parameters%beam_energy_out = energy_out_beam
        output_parameters%ext_energy_out = energy_out_ext_diff
        
        if(job_params%debug >= 1) then
            write(101,*)'------------------------------------------------------'
            write(101,'(a41,f16.8)')'energy in (ill. geom. cross sec.): ', energy_in
            write(101,'(a41,f16.8)')'beam energy out: ',energy_out_beam
            write(101,'(a41,f16.8)')'absorbed beam energy: ',output_parameters%abs
            write(101,'(a41,f16.8)')'ext diff energy out: ',energy_out_ext_diff
            write(101,'(a41,f16.8,a2)')'beam energy conservation: ',(energy_out_beam+output_parameters%abs)/energy_in*100,' %'
            write(101,'(a41,f16.8,a2)')'ext diff energy conservation: ',energy_out_ext_diff/energy_in*100,' %'
            
            print'(a40,f16.8)','energy in (ill. geom. cross sec.): ', energy_in
            print'(a40,f16.8)','beam energy out: ',energy_out_beam
            print'(a40,f16.8)','absorbed beam energy: ',output_parameters%abs
            print'(a40,f16.8)','ext diff energy out: ',energy_out_ext_diff
            print'(a40,f16.8,a2)','beam energy conservation: ',(energy_out_beam+output_parameters%abs)/energy_in*100,' %'
            print'(a40,f16.8,a2)','ext diff energy conservation: ',energy_out_ext_diff/energy_in*100,' %'
            
            print*,'========== end sr energy_checks'
        end if
        ! stop
        
    end subroutine
    
    subroutine beam_loop(   beam_outbeam_tree, &
        beam_outbeam_tree_counter, &
        ext_diff_outbeam_tree, &
        output_parameters, &
        job_params, &
        geometry, &
        beam_geometry, &
        beam_inc)
        
        ! beam_loop
        ! performs the main beam loop
        ! to be called from the main program
        ! for an incident beam and a geometry, this returns the near field
        !   in the form of the beam_outbeam_tree, as well as the externally
        !   diffracted near field, in the form of the ext_diff_outbeam_tree
        
        ! inputs
        ! character(100), intent(in) :: afn ! apertures filename
        type(outbeamtype), dimension(:), allocatable, intent(out) :: beam_outbeam_tree ! outgoing beams from the beam tracing
        type(outbeamtype), dimension(:), allocatable, intent(out) :: ext_diff_outbeam_tree
        integer(8), intent(out) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
        type(output_parameters_type), intent(inout) :: output_parameters 
        type(job_parameters_type), intent(in) :: job_params
        type(geometry_type), intent(in) :: geometry
        type(geometry_type), intent(in) :: beam_geometry
        type(beam_type), intent(inout) :: beam_inc
        
        real(8) start, finish, start1, finish1 ! cpu timing variables
        integer(8) i, j
        
        ! new beam tree
        type(beam_type), dimension(:), allocatable :: beam_tree
        integer(8) num_beams
        integer(8) i_start, i_end
        type(beam_type) beam ! current beam to be traced
        integer(8) num_threads
        
        ! init
        start = 0d0
        start1 = 0d0
        finish = 0d0
        finish1 = 0d0

        if(job_params%timing) then
            start = omp_get_wtime()
        endif
        
        if(job_params%debug >= 1) then
            print*,'start beam loop...'
        end if
        
        beam_outbeam_tree_counter = 0 ! counts the current number of beam outbeams
        allocate(beam_outbeam_tree(1:1000000)) ! set to 1000000 as guess for max outbeams
        
        num_beams = 0 ! total number of beams in the beam tree
        allocate(beam_tree(1:10000)) ! allocate some space to hold beams to be traced (might need to add more space later)
        
        if(job_params%debug >= 2) then
            print*,'incidence:'
            if(job_params%timing) then
                start1 = omp_get_wtime()
            end if
        end if

        call recursion_inc(beam_inc,geometry,job_params,beam_geometry,ext_diff_outbeam_tree) ! do the initial incidence

        call add_to_beam_tree_external(beam_tree,beam_inc,num_beams,geometry,job_params) ! add beams to be propagated to the tree
        
        call add_to_outbeam_tree(beam_outbeam_tree,beam_outbeam_tree_counter,beam_inc,job_params) ! add externally reflected beams to diffraction tree
        
        call get_geo_cross_section(geometry,beam_inc,output_parameters) ! for the first recursion, get the illuminated geometric cross section
        
        if(job_params%timing) then
            if(job_params%debug >= 2) then
                finish1 = omp_get_wtime()
                print*,'======================================'
                print'(a,f16.8,a)',"incidence - time taken: ",finish1-start1," secs"
                write(101,'(a,f16.8,a)')"incidence - time taken: ",finish1-start1," secs"
                print*,'======================================'
            end if
        end if

        ! main loop
        i_start = 1 ! entry in beam tree to start at
        do i = 1, job_params%rec ! for each recursion
            i_end = num_beams ! entry in beam tree to stop at
            
            if(job_params%debug >= 2) then
                print*,'internal recursion: ',i
                if(job_params%timing) then
                    start1 = omp_get_wtime()
                end if
            end if

            ! get number of threads required
            if(job_params%is_multithreaded) then ! if multithreading enabled
                num_threads = min(omp_get_max_threads(),i_end-i_start+1) ! set threads
            else ! if multithreading disabled
                num_threads = 1
            end if

            ! loop over each beam for this recursion
            !$omp parallel num_threads(num_threads) private(beam)
            !$omp do
            do j = i_start, i_end ! for each entry in the beam tree that belongs to this recursion
                beam = beam_tree(j) ! get a beam from the beam_tree
                ! propagate this beam
                if(beam%is_int) then ! if the beam is internally propagating
                    call recursion_int(beam,geometry,job_params) ! propagate the beam and populate the beam structure
                else ! if the beam is externally propagating
                    print*,'external propagation not supported yet'
                    cycle
                end if ! end: if the beam is internally propagating
                !$omp critical     
                beam_tree(beam%id) = beam ! update the information about the propagated beam in the beam tree
                if(job_params%debug >= 2) call print_beam_info(beam,job_params) ! print some info about the beam
                if(job_params%ibi > 0d0) output_parameters%abs = output_parameters%abs + beam%abs ! udpate absorption cross section
                ! add the new beams to the beam tree (need to also add the new information about the current beam)
                call add_to_beam_tree_internal(beam_tree,beam,num_beams,geometry,job_params)
                ! add the outgoing surface field to the diffraction structure
                call add_to_outbeam_tree(beam_outbeam_tree,beam_outbeam_tree_counter,beam,job_params)
                !$omp end critical
            end do ! end: for each entry in the beam tree that belongs to this recursion
            !$omp end parallel
            
            if(job_params%timing) then
                if(job_params%debug >= 2) then
                    finish1 = omp_get_wtime()
                    print*,'======================================'
                    print'(a,i3,a,f16.8,a)',"recursion",i," - time taken: ",finish1-start1," secs"
                    write(101,'(a,i3,a,f16.8,a)')"recursion",i," - time taken: ",finish1-start1," secs"
                    print*,'======================================'
                end if
            end if
            
            i_start = i_end + 1 ! update the starting index for the next recursion
        end do ! end: for each recursion
        
        if(job_params%debug >= 1) then
            if(job_params%timing) then
                finish = omp_get_wtime()
                print*,'=========='
                print'(a,f16.8,a)',"end beam loop - time elapsed: ",finish-start," secs"   
            end if
        end if
        
        call get_beamtree_vert(beam_outbeam_tree, beam_outbeam_tree_counter,geometry)    
        
        call energy_checks( beam_outbeam_tree, &
        beam_outbeam_tree_counter, &
        output_parameters, &
        ext_diff_outbeam_tree, &
        job_params, &
        geometry)
        
        if(job_params%debug >= 2) then
            print'(a)','memory usage breakdown (per mpi process):'
            print'(a,f8.2,a)','particle geometry: ',real(sizeof(geometry))/1048576d0,' mb'
            print'(a,f8.2,a)','beam tree: ',real(sizeof(beam_tree))/1048576d0,' mb'
            print'(a,f8.2,a)','ext. diffraction tree: ',real(sizeof(ext_diff_outbeam_tree))/1048576d0,' mb'
            print'(a,f8.2,a)','outbeam tree: ',real(sizeof(beam_outbeam_tree))/1048576d0,' mb'
            ! print'(a)','note: this is an underestimate of the total memory usage.'
            print'(a)',' =========='
        end if
        
        ! stop
        
    end subroutine
    
    subroutine get_beamtree_vert(beam_outbeam_tree, beam_outbeam_tree_counter, geometry)
        
        ! get_beamtree_vert
        ! adds information about the particle geometry to the diffraction tree
        ! to be removed at a later date
        
        type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
        integer(8), intent(in) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
        type(geometry_type), intent(in) :: geometry
        
        integer(8) i, fi
        
        do i = 1, beam_outbeam_tree_counter ! for each beam
            fi = beam_outbeam_tree(i)%fout ! get the face id
            beam_outbeam_tree(i)%verts = transpose(geometry%v(geometry%f(fi)%vi(:),:))
        end do
        
    end subroutine
    
    subroutine find_vis_int(in_beam, dist_beam, id_beam, is_shad, beam, rot_geometry)
        
        ! find_vis_int
        ! finds the facets illuminated by a an internally propagating beam
        ! for a given particle orientation with assumed propagation along the -z axis,
        ! computes the internally illuminated facets for a given illuminating aperture
        ! uses a z-buffer technique to increase speed, among a few other tricks
        
        logical, dimension(:), allocatable, intent(out) :: in_beam
        real(8), dimension(:), allocatable, intent(out) :: dist_beam
        integer(8), dimension(:), allocatable, intent(out) :: id_beam
        logical, dimension(:), allocatable, intent(out) :: is_shad
        type(beam_type), intent(in) :: beam
        type(geometry_type), intent(in) :: rot_geometry
        
        integer(8) i, k, m, fi, vi
        integer(8) j
        logical, dimension(:), allocatable :: is_beam
        logical, dimension(:), allocatable :: is_vis
        integer(8), dimension(:), allocatable :: f3 ! bounding box ids
        integer(8), dimension(:,:), allocatable :: f4 ! fuzzy bounding box ids
        real(8), dimension(:), allocatable :: dist_to_bb, dist_to_fbb
        integer(8) bb
        logical within_bounds
        real(8) vecb1, vecb2
        real(8) edge_norm1, edge_norm2
        real(8) edge_check
        real(8) beamxmax0, beamxmin0, beamymax0, beamymin0
        real(8) beamxmax, beamxmin, beamymax, beamymin
        logical, dimension(:,:), allocatable :: f5
        real(8) start, finish
        integer(8), dimension(:), allocatable :: mapping
        type(geometry_type) bb_geometry ! bounding box geometry
        
        ! ################################
        ! start new ray tracing algorithm
        
        call cpu_time(start)
        
        ! use the current crystal vertices to create some bounding boxes in x-y plane
        call beam_aligned_bounding_boxes(rot_geometry,bb_geometry)
        call compute_geometry_midpoints(bb_geometry)
        call compute_geometry_areas(bb_geometry)
        
        allocate(f3(1:rot_geometry%nf)) ! array to hold index of bounding box that each face belongs to
        allocate(f4(1:rot_geometry%nf,1:3)) ! array to hold index of fuzzy bounding box that each face belongs to
        allocate(dist_to_bb(1:bb_geometry%nf)) ! array to hold the distance of a given vertex to each bounding box
        allocate(dist_to_fbb(1:bb_geometry%nf)) ! array to hold the distance of a given vertex to each bounding box
        
        do i = 1, rot_geometry%nf ! for each face
            dist_to_bb(:) = sqrt((bb_geometry%f(:)%mid(1) - rot_geometry%f(i)%mid(1))**2 + (bb_geometry%f(:)%mid(2) - rot_geometry%f(i)%mid(2))**2) ! distance to each bb
            f3(i) = minloc(dist_to_bb,1) ! record which bounding box midpoint this facet was closest to
            do j = 1, 3 ! for each vertex
                dist_to_fbb(:) = sqrt((bb_geometry%f(:)%mid(1) - rot_geometry%v(rot_geometry%f(i)%vi(j),1))**2 + (bb_geometry%f(:)%mid(2) - rot_geometry%v(rot_geometry%f(i)%vi(j),2))**2) ! distance to each fuzzy bb
                f4(i,j) = minloc(dist_to_fbb,1) ! record which bounding box midpoint this facet vertex was closest to
            end do 
        end do
        
        ! more optimisation
        allocate(f5(1:rot_geometry%nf,1:bb_geometry%nf))
        f5 = .false. ! init
        do i = 1, rot_geometry%nf
            do j = 1, bb_geometry%nf
                if(f4(i,1) .eq. j .or. f4(i,2) .eq. j .or. f4(i,3) .eq. j) f5(i,j) = .true.
            end do
        end do
        
        ! allocate and init
        allocate(is_vis(1:rot_geometry%nf))
        allocate(is_shad(1:rot_geometry%nf))
        allocate(in_beam(1:rot_geometry%nf))
        allocate(is_beam(1:rot_geometry%nf))
        allocate(id_beam(1:rot_geometry%nf))
        allocate(dist_beam(1:rot_geometry%nf))
        allocate(mapping(1:rot_geometry%nf))
        
        is_vis = .true.
        is_shad = .false.
        in_beam = .false.
        is_beam = .false.
        id_beam = 0
        dist_beam = 0
        mapping = 0
        
        ! record which facets are part of the illuminating beam
        do i = 1, beam%nf_in
            is_beam(beam%field_in(i)%fi) = .true.
            mapping(beam%field_in(i)%fi) = i ! map each illuminating facet to its position in the beam data structure
        end do
        
        ! need to add a check to ignore facets outside beam max/min here
        ! init some stuff
        fi = beam%field_in(1)%fi ! get the first facet in the beam
        vi = rot_geometry%f(fi)%vi(1) ! get the vertex id of the first vertex in this facet
        beamxmax = rot_geometry%v(vi,1) ! save the x coordinate
        beamxmin = rot_geometry%v(vi,1) ! save the x coordinate
        beamymax = rot_geometry%v(vi,2) ! save the y coordinate
        beamymin = rot_geometry%v(vi,2) ! save the y coordinate
        
        do i = 1, beam%nf_in ! for each face in the incident beam
            fi = beam%field_in(i)%fi ! get the facet id
            do j = 1, rot_geometry%f(fi)%nv ! for each vertex in this face
                vi = rot_geometry%f(fi)%vi(j) ! get the vertex id
                beamxmax0 = rot_geometry%v(vi,1) ! get the x coordinate
                beamxmin0 = rot_geometry%v(vi,1) ! get the x coordinate
                beamymax0 = rot_geometry%v(vi,2) ! get the y coordinate
                beamymin0 = rot_geometry%v(vi,2) ! get the y coordinate
                if(beamxmax0 .gt. beamxmax) beamxmax = beamxmax0
                if(beamxmin0 .lt. beamxmin) beamxmin = beamxmin0
                if(beamymax0 .gt. beamymax) beamymax = beamymax0
                if(beamymin0 .lt. beamymin) beamymin = beamymin0 
            end do
        end do
        
        do i = 1, rot_geometry%nf
            if(rot_geometry%f(i)%mid(1) .lt. beamxmin) then
                is_vis(i) = .false.
            else if(rot_geometry%f(i)%mid(1) .gt. beamxmax) then
                is_vis(i) = .false.
            else if(rot_geometry%f(i)%mid(2) .lt. beamymin) then
                is_vis(i) = .false.
            else if(rot_geometry%f(i)%mid(2) .gt. beamymax) then
                is_vis(i) = .false.
            end if
        end do
        
        do m = 1, rot_geometry%nf ! for each facet m
            if(is_vis(m) .eqv. .false.) then ! if facet isnt visible
                ! do nothing
            else
                if(rot_geometry%ap(rot_geometry%f(m)%ap)%n(3) .gt. -0.01) then ! if aperture is downfacing
                    is_vis(m) = .false. ! set not visible
                else ! if aperture was facing towards incidence
                    bb = f3(m) ! get bounding box id
                    do j = 1, rot_geometry%nf ! for each potentially blocking facet j
                        ! if(rotatednorm(face2(j),3) .lt. 0.01) then ! sign flip here for internal
                        ! if(f4(j,1) .eq. bb .or. f4(j,2) .eq. bb .or. f4(j,3) .eq. bb) then ! 
                        if(f5(j,bb)) then ! 
                            ! if(any(f4(j,1:3)) .eq. bb) then ! if blocker was in fuzzy bounding box
                            if(j .ne. m) then ! ignore self-block
                                if(rot_geometry%n(rot_geometry%f(j)%ni,3) .lt. 0.01) then ! if down-facing, sign flip here for internal
                                    ! do nothing
                                else ! if up-facing
                                    if (rot_geometry%f(m)%mid(3) .gt. rot_geometry%f(j)%mid(3)) then ! if potential blocker was behind facet m
                                        ! do nothing
                                    else ! if potential blocker was in front of facet m
                                        ! do bounded surface edge check
                                        within_bounds = .true. ! assume centroid of facet m is within the bounded surface of potentially blocking facet j
                                        do k = 1, rot_geometry%f(j)%nv ! looping over vertices of potentially blocking facet j
                                            ! compute edge normal
                                            if(k .eq. rot_geometry%f(j)%nv) then
                                                edge_norm2 = -rot_geometry%v(rot_geometry%f(j)%vi(1),1) + rot_geometry%v(rot_geometry%f(j)%vi(rot_geometry%f(j)%nv),1) ! cross product of edge vector with reverse beam direction
                                                edge_norm1 = rot_geometry%v(rot_geometry%f(j)%vi(1),2) - rot_geometry%v(rot_geometry%f(j)%vi(rot_geometry%f(j)%nv),2)
                                            else
                                                edge_norm2 = -rot_geometry%v(rot_geometry%f(j)%vi(k+1),1) + rot_geometry%v(rot_geometry%f(j)%vi(k),1)
                                                edge_norm1 = rot_geometry%v(rot_geometry%f(j)%vi(k+1),2) - rot_geometry%v(rot_geometry%f(j)%vi(k),2)
                                            end if
                                            vecb1 = rot_geometry%f(m)%mid(1) - rot_geometry%v(rot_geometry%f(j)%vi(k),1) ! vector from vertex k of potential blocker j to centroid of facet m
                                            vecb2 = rot_geometry%f(m)%mid(2) - rot_geometry%v(rot_geometry%f(j)%vi(k),2)
                                            edge_check = vecb1*edge_norm1 + vecb2*edge_norm2 ! dot product of edge vector with vetor b
                                            if(edge_check .gt. 0) within_bounds = .false. ! if edge check fails, centroid of facet m is not within the bounded surface of facet j
                                        end do
                                        if(within_bounds .eqv. .false.) then ! if facet m is not within bounded surface of facet j
                                            !do nothing
                                        else
                                            if(is_beam(j)) then ! if facet j was part of the illuminating surface
                                                if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                    ! check to see which is the closest
                                                    if(rot_geometry%f(m)%mid(3) .gt. rot_geometry%f(id_beam(m))%mid(3)) then ! if facet j was closer than the previous blocker
                                                        in_beam(m) = .true. ! set to be within beam
                                                        id_beam(m) = mapping(j) ! record the blocking facet
                                                        dist_beam(m) = rot_geometry%f(j)%mid(3) - rot_geometry%f(m)%mid(3) ! record the distance from centroid of blocker to centroid of facet m
                                                    else
                                                        ! do nothing                                     
                                                    end if
                                                else ! if this is the first time finding a blocking beam facet
                                                    in_beam(m) = .true. ! set to be within beam
                                                    id_beam(m) = mapping(j) ! record the position of the blocking facet in the beam data struct
                                                    dist_beam(m) = rot_geometry%f(j)%mid(3) - rot_geometry%f(m)%mid(3) ! record the distance from centroid of blocker to centroid of facet m
                                                end if
                                            else
                                                if(rot_geometry%f(m)%ap .eq. rot_geometry%f(j)%ap) then ! if facet j and facet m belong to the same aperture
                                                    is_shad(m) = .true. ! set m as a shadow facet and continue to search for blockers
                                                else ! if facet j and facet m dont belong to the same aperture
                                                    if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                        if(rot_geometry%f(j)%mid(3) .gt. rot_geometry%f(id_beam(m))%mid(3)) then ! if facet j was behind than the blocking beam
                                                            ! do nothing
                                                        else
                                                            is_vis(m) = .false. ! set m as not in the shadow and has been blocked by non-illuminating facet
                                                        end if
                                                    end if
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end do
                end if
            end if
        end do
        
        call cpu_time(finish)
        
    end subroutine
    
    subroutine find_vis_incidence(in_beam, is_vis, dist_beam, id_beam, is_shad, beam, geometry, beam_geometry)
        
 
        
        logical, dimension(:), allocatable, intent(out) :: in_beam
        real(8), dimension(:), allocatable, intent(out) :: dist_beam
        integer(8), dimension(:), allocatable, intent(out) :: id_beam
        logical, dimension(:), allocatable, intent(out) :: is_shad
        type(beam_type), intent(in) :: beam
        type(geometry_type), intent(in) :: geometry
        type(geometry_type), intent(in) :: beam_geometry
        
        integer(8) i, k, m
        integer(8) j
        logical, dimension(:), allocatable :: is_beam
        logical, dimension(:), allocatable, intent(out) :: is_vis
        integer(8), dimension(:), allocatable :: f3 ! bounding box ids
        integer(8), dimension(:,:), allocatable :: f4 ! fuzzy bounding box ids
        real(8), dimension(:), allocatable :: dist_to_bb, dist_to_fbb
        integer(8) bb
        logical within_bounds
        real(8) vecb1, vecb2
        real(8) edge_norm1, edge_norm2
        real(8) edge_check
        logical, dimension(:,:), allocatable :: f5
        real(8) start, finish
        integer(8), dimension(:), allocatable :: mapping
        type(geometry_type) bb_geometry ! bounding box geometry
        
        ! ################################
        ! start new ray tracing algorithm
        
        call cpu_time(start)
        
        ! use the current crystal vertices to create some bounding boxes in x-y plane
        call beam_aligned_bounding_boxes(geometry,bb_geometry)
        call compute_geometry_midpoints(bb_geometry)
        call compute_geometry_areas(bb_geometry)
        
        allocate(f3(1:geometry%nf)) ! array to hold index of bounding box that each face belongs to
        allocate(f4(1:geometry%nf,1:3)) ! array to hold index of fuzzy bounding box that each face belongs to
        allocate(dist_to_bb(1:bb_geometry%nf)) ! array to hold the distance of a given vertex to each bounding box
        allocate(dist_to_fbb(1:bb_geometry%nf)) ! array to hold the distance of a given vertex to each bounding box
        
        do i = 1, geometry%nf ! for each face in the particle geometry
            dist_to_bb(:) = sqrt((bb_geometry%f(:)%mid(1) - geometry%f(i)%mid(1))**2 + (bb_geometry%f(:)%mid(2) - geometry%f(i)%mid(2))**2) ! distance to each bb
            f3(i) = minloc(dist_to_bb,1) ! record which bounding box midpoint this facet was closest to
            do j = 1, 3 ! for each vertex
                dist_to_fbb(:) = sqrt((bb_geometry%f(:)%mid(1) - geometry%v(geometry%f(i)%vi(j),1))**2 + (bb_geometry%f(:)%mid(2) - geometry%v(geometry%f(i)%vi(j),2))**2) ! distance to each fuzzy bb
                f4(i,j) = minloc(dist_to_fbb,1) ! record which bounding box midpoint this facet vertex was closest to
            end do 
        end do
        
        ! more optimisation
        allocate(f5(1:geometry%nf,1:bb_geometry%nf))
        f5 = .false. ! init
        do i = 1, geometry%nf
            do j = 1, bb_geometry%nf
                if(f4(i,1) .eq. j .or. f4(i,2) .eq. j .or. f4(i,3) .eq. j) f5(i,j) = .true.
            end do
        end do
        
        ! allocate and init
        allocate(is_vis(1:geometry%nf))
        allocate(is_shad(1:geometry%nf))
        allocate(in_beam(1:geometry%nf))
        allocate(is_beam(1:geometry%nf))
        allocate(id_beam(1:geometry%nf))
        allocate(dist_beam(1:geometry%nf))
        allocate(mapping(1:geometry%nf))
        
        is_vis = .true.
        is_shad = .false.
        in_beam = .true.
        is_beam = .false.
        id_beam = 0
        dist_beam = 0
        mapping = 0
        
        ! record which facets are part of the illuminating beam
        do i = 1, beam%nf_in
            is_beam(beam%field_in(i)%fi) = .true.
            mapping(beam%field_in(i)%fi) = i ! map each illuminating facet to its position in the beam data structure
        end do
        
        do m = 1, geometry%nf ! for each facet m
            if(geometry%ap(geometry%f(m)%ap)%n(3) .lt. 0.01) then ! if aperture is downfacing (flip for external)
                is_vis(m) = .false. ! set not visible
                in_beam(m) = .false. ! set not in beam
            else ! if aperture was facing towards incidence
                bb = f3(m) ! get bounding box id
                do j = 1, geometry%nf ! for each potentially blocking facet j
                    if(f5(j,bb)) then ! 
                        if(j .ne. m) then ! ignore self-block
                            if(geometry%n(geometry%f(j)%ni,3) .lt. 0.01) then ! if down-facing
                                ! do nothing
                            else ! if up-facing
                                if (geometry%f(m)%mid(3) .gt. geometry%f(j)%mid(3)) then ! if potential blocker was behind facet m
                                    ! do nothing
                                else ! if potential blocker was in front of facet m
                                    ! do bounded surface edge check
                                    within_bounds = .true. ! assume centroid of facet m is within the bounded surface of potentially blocking facet j
                                    do k = 1, geometry%f(j)%nv ! looping over vertices of potentially blocking facet j
                                        ! compute edge normal
                                        if(k .eq. geometry%f(j)%nv) then ! if its the last vertex in the facet, loop back to the first
                                            edge_norm2 = -geometry%v(geometry%f(j)%vi(1),1) + geometry%v(geometry%f(j)%vi(geometry%f(j)%nv),1) ! cross product of edge vector with reverse beam direction
                                            edge_norm1 = +geometry%v(geometry%f(j)%vi(1),2) - geometry%v(geometry%f(j)%vi(geometry%f(j)%nv),2)
                                        else
                                            edge_norm2 = -geometry%v(geometry%f(j)%vi(k+1),1) + geometry%v(geometry%f(j)%vi(k),1)
                                            edge_norm1 = +geometry%v(geometry%f(j)%vi(k+1),2) - geometry%v(geometry%f(j)%vi(k),2)
                                        end if
                                        vecb1 = geometry%f(m)%mid(1) - geometry%v(geometry%f(j)%vi(k),1) ! vector from vertex k of potential blocker j to centroid of facet m
                                        vecb2 = geometry%f(m)%mid(2) - geometry%v(geometry%f(j)%vi(k),2)
                                        edge_check = vecb1*edge_norm1 + vecb2*edge_norm2 ! dot product of edge vector with vetor b (sign flip)
                                        if(edge_check .gt. 0) within_bounds = .false. ! if edge check fails, centroid of facet m is not within the bounded surface of facet j
                                    end do
                                    if(within_bounds .eqv. .false.) then ! if facet m is not within bounded surface of facet j
                                        !do nothing
                                    else
                                        if(geometry%f(m)%ap .eq. geometry%f(j)%ap) then ! if facet j and facet m belong to the same aperture
                                            is_shad(m) = .true. ! set m as a shadow facet and continue to search for blockers
                                        else ! if facet j and facet m dont belong to the same aperture
                                            is_vis(m) = .false. ! set m as not in the shadow and has been blocked by non-illuminating facet
                                            in_beam(m) = .false. ! not in beam
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end do
            end if
        end do
        
        do m = 1, geometry%nf ! for each facet m
            if(in_beam(m)) then ! if in beam
                dist_beam(m) = beam_geometry%f(1)%mid(3) - geometry%f(m)%mid(3) ! compute the distance travelled ! cheat for simple incident beam with 1 facet
                id_beam(m) = 1 ! cheat for simple incident beam with 1 facet
            end if
        end do
        
        call cpu_time(finish)
        
    end subroutine
    
    subroutine get_rot_matrix(rot, vk7, ev1, ev3)
        
        ! get_rot_matrix
        ! get rotation matrix for rotation about vector ev3 from plane perp to ev1 to plane perp to vk7
        
        real(8), intent(in) :: vk7(1:3) ! reflected e-perp vector from each facet
        real(8), intent(in) :: ev1(1:3) ! incident e-perp vector
        real(8), intent(in) :: ev3(1:3) ! incident propagation vector
        real(8), intent(out) :: rot(1:2,1:2) ! rotation matrix
        
        real(8) myfunc_arg
        real(8) evo2(1:3)
        real(8) det ! normalisation factor
        
        myfunc_arg = dot_product(vk7,ev1)
        
        call cross(ev3,ev1,evo2,.true.)
        
        rot(1,1) = myfunc_arg
        rot(1,2) = -dot_product(vk7,evo2)
        rot(2,1) = +dot_product(vk7,evo2)
        rot(2,2) = myfunc_arg
        
        ! normalise (improvement over macke)
        det = rot(1,1)*rot(2,2) - rot(1,2)*rot(2,1)
        rot(:,:) = rot(:,:) / sqrt(abs(det))
        
    end subroutine
    
    subroutine get_geo_cross_section(geometry,inc_beam,output_parameters)

        ! get_geo_cross_section
        ! gets the illuminated geometric cross section for a beam
        ! currently adds this directly to the output parameters, but could be
        !   used more generally
        
        type(output_parameters_type), intent(inout) :: output_parameters
        type(geometry_type), intent(in) :: geometry
        type(beam_type), intent(in) :: inc_beam
        
        integer(8) i, fi
        real(8) geo_cross_sec
        real(8) area
        real(8) norm(1:3)
        
        geo_cross_sec = 0 ! init
        do i = 1, inc_beam%nf_out ! for each face illuminated by the incident beam
            fi = inc_beam%field_out(i)%fi
            area = geometry%f(fi)%area
            norm = geometry%n(geometry%f(fi)%ni,:)
            geo_cross_sec = geo_cross_sec + area*norm(3)
        end do 
        
        output_parameters%geo_cross_sec = geo_cross_sec
        
    end subroutine
    
    subroutine beam_aligned_bounding_boxes(geometry,bb_geometry)

    ! beam_aligned_bounding_boxes
    ! makes some 2-dimensional bounding boxes to help speed up some other subroutines
    ! this method was quick to implement but should be replaced in the future by a
    !   method which does not need to be re-called for every beam that is propagated
        
        type(geometry_type), intent(in) :: geometry
        type(geometry_type), intent(out) :: bb_geometry
        
        real(8) min_x, min_y, max_x, max_y, min_z, max_z
        integer(8), parameter :: bounding_box_x_dim = 8 ! bounding box x dimension
        integer(8), parameter :: bounding_box_y_dim = 8 ! bounding box y dimension
        real(8), parameter :: fac = 1.1
        real(8), dimension(1:bounding_box_x_dim+1, 1:bounding_box_y_dim+1) :: xvals, yvals ! bounding box x and y vertices
        integer(8) i,j, vert_counter, face_counter
        integer(8) num_verts, num_faces
        
        ! get the max xyz and min xy coordinates
        max_x = maxval(geometry%v(:,1))
        min_x = minval(geometry%v(:,1))
        max_y = maxval(geometry%v(:,2))
        min_y = minval(geometry%v(:,2))
        max_z = maxval(geometry%v(:,3))
        min_z = minval(geometry%v(:,3))
        
        !print*,'making bounding boxes'
        
        num_verts = (bounding_box_x_dim+1)*(bounding_box_y_dim+1)
        num_faces = (bounding_box_x_dim*bounding_box_y_dim)
        
        ! allocate arrays to hold outer bounding box
        allocate(bb_geometry%v(1:num_verts,1:3))
        allocate(bb_geometry%f(1:num_faces))
        bb_geometry%f(:)%nv = 4
        bb_geometry%nf = num_faces
        do i = 1, bb_geometry%nf
            allocate(bb_geometry%f(i)%vi(1:4))
        end do
        
        bb_geometry%nv = num_verts
        bb_geometry%nf = num_faces
        
        !print*,'size(boundingboxv,1)',size(boundingboxv,1)
        !print*,'boundingboxfsize',boundingboxfsize
        
        ! make the bounding box vertices (meshgrid style)
        bb_geometry%v(:,3) = 0 ! init
        vert_counter = 0
        do i = 1, bounding_box_x_dim + 1
            do j = 1, bounding_box_y_dim + 1
                xvals(i,j) = fac*(min_x + (max_x - min_x)/bounding_box_x_dim*(i-1))
                yvals(i,j) = fac*(min_y + (max_y - min_y)/bounding_box_x_dim*(j-1))
                !print*,'xvals(i,j)',xvals(i,j)
                !print*,'yvals(i,j)',yvals(i,j)
                vert_counter = vert_counter + 1
                bb_geometry%v(vert_counter,1) = xvals(i,j)
                bb_geometry%v(vert_counter,2) = yvals(i,j)
            end do
        end do
        
        !print*,'total vertices in vertex array',vert_counter
        
        ! make vertex and face arrays using the bottom left vertex as starting point (moving clockwise)
        face_counter = 0
        do i = 1, bounding_box_x_dim
            do j = 1, bounding_box_y_dim
                face_counter = face_counter + 1
                bb_geometry%f(face_counter)%vi(1) = (i-1)*(bounding_box_y_dim+1) + j
                bb_geometry%f(face_counter)%vi(2) = (i-1)*(bounding_box_y_dim+1) + j + 1
                bb_geometry%f(face_counter)%vi(3) = i*(bounding_box_y_dim+1) + j + 1
                bb_geometry%f(face_counter)%vi(4) = i*(bounding_box_y_dim+1) + j
            end do
        end do
        
    end subroutine
    
    subroutine add_to_outbeam_tree(beam_outbeam_tree,beam_outbeam_tree_counter,beam,job_params)
        
        ! add_to_outbeam_tree
        ! adds the external field after a beam has propagated to the outbeam tree, ready for diffraction

        integer(8), intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
        type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
        type(beam_type), intent(in) :: beam ! current beam to be traced
        type(job_parameters_type), intent(in) :: job_params

        integer(8) i, counter
        
        counter = 0 ! counts the number of beams added to the tree
        
        do i = 1, beam%nf_out ! for each illuminated facet in this beam structure
            if(.not. beam%field_out(i)%is_tir) then ! if it wasnt total internal reflection, add to outbeams
                beam_outbeam_tree_counter = beam_outbeam_tree_counter + 1 ! update outbeam tree counter
                counter = counter + 1
                if(beam_outbeam_tree_counter .gt. size(beam_outbeam_tree,1)) then
                    print*,'error: need more space in outbeam_tree. please increase in sr init'
                end if
                beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(:,:) = beam%field_out(i)%ampl_ext(:,:) ! external amplitude matrix
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(:) = beam%field_out(i)%prop_int(:) ! incident propagation direction
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(:) = beam%field_out(i)%prop_ext(:) ! outgoing propagation direction
                beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(:) = beam%field_out(i)%e_perp(:) ! perpendicular field vector
                beam_outbeam_tree(beam_outbeam_tree_counter)%fout = beam%field_out(i)%fi ! facet id
            end if
        end do

        if(job_params%debug >= 3) print'(a,i6,a,i8,a,i8,a)','beam ',beam%id,': added ',counter,' beams to outbeam tree -->',beam_outbeam_tree_counter,' total outbeams'

    end subroutine
    
end module beam_loop_mod
