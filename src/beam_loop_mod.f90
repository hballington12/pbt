! beam_loop_mod.f90
! module for the beam tracing loop

module beam_loop_mod

use misc_submod
use types_mod
use omp_lib

implicit none

real(8) ext_cross_section ! extinction cross section

contains

subroutine energy_checks(   beam_outbeam_tree, &
                            beam_outbeam_tree_counter, &
                            norm, &
                            face2, &
                            faceAreas, &
                            energy_out_beam, &
                            energy_in, &
                            ext_diff_outbeam_tree, &
                            energy_out_ext_diff, &
                            energy_abs_beam, &
                            job_params)

    type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
    integer, intent(in) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
    real(8), dimension(:,:), allocatable :: norm ! face normals
    integer, dimension(:), allocatable :: face2 ! face normal ID of each face
    real(8), dimension(:), allocatable :: faceAreas ! area of each facet
    real(8), intent(in) :: energy_in
    real(8), intent(out) :: energy_out_beam
    real(8), intent(out) :: energy_abs_beam
    real(8), intent(out) :: energy_out_ext_diff
    type(outbeamtype), dimension(:), allocatable, intent(in) :: ext_diff_outbeam_tree
    type(job_parameters_type), intent(in) ::  job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details

    integer i
    integer face_id
    real(8) prop(1:3)
    real(8) normal(1:3)
    real(8) cos_theta
    real(8) intensity_out
    complex(8) ampl(1:2,1:2)
    real(8) area

    if(job_params%debug >= 1) then
        print*,'========== start sr energy_checks'
    end if

    energy_out_beam = 0
    energy_abs_beam = ext_cross_section

    do i = 1, beam_outbeam_tree_counter
        prop = beam_outbeam_tree(i)%prop_out
        face_id = beam_outbeam_tree(i)%FOut
        ampl = beam_outbeam_tree(i)%ampl
        normal = norm(face2(face_id),1:3)
        area = faceAreas(face_id)
        cos_theta = dot_product(prop,normal)
        intensity_out = real(0.5*(   ampl(1,1)*conjg(ampl(1,1)) + &
                                ampl(1,2)*conjg(ampl(1,2)) + &
                                ampl(2,1)*conjg(ampl(2,1)) + &
                                ampl(2,2)*conjg(ampl(2,2))))
                                
        if (isnan(intensity_out)) then
            print*,'oh dear - nan ibeam = ',i
            beam_outbeam_tree(i)%ampl = 0
        else
            energy_out_beam = energy_out_beam + intensity_out*area*cos_theta
        end if
        
    end do

    energy_out_ext_diff = 0

    do i = 1, size(ext_diff_outbeam_tree,1)
        prop = ext_diff_outbeam_tree(i)%prop_out
        face_id = ext_diff_outbeam_tree(i)%FOut
        ampl = ext_diff_outbeam_tree(i)%ampl
        normal = norm(face2(face_id),1:3)
        area = faceAreas(face_id)
        cos_theta = -dot_product(prop,normal)
        intensity_out = real(0.5*(   ampl(1,1)*conjg(ampl(1,1)) + &
                                ampl(1,2)*conjg(ampl(1,2)) + &
                                ampl(2,1)*conjg(ampl(2,1)) + &
                                ampl(2,2)*conjg(ampl(2,2))))
        energy_out_ext_diff = energy_out_ext_diff + intensity_out*area*cos_theta
    end do

    if(job_params%debug >= 1) then
        write(101,*)'------------------------------------------------------'
        write(101,'(A41,f16.8)')'energy in (ill. geom. cross sec.): ', energy_in
        write(101,'(A41,f16.8)')'beam energy out: ',energy_out_beam
        write(101,'(A41,f16.8)')'absorbed beam energy: ',ext_cross_section
        write(101,'(A41,f16.8)')'ext diff energy out: ',energy_out_ext_diff
        write(101,'(A41,f16.8,A2)')'beam energy conservation: ',(energy_out_beam+ext_cross_section)/energy_in*100,' %'
        write(101,'(A41,f16.8,A2)')'ext diff energy conservation: ',energy_out_ext_diff/energy_in*100,' %'

        print'(A40,f16.8)','energy in (ill. geom. cross sec.): ', energy_in
        print'(A40,f16.8)','beam energy out: ',energy_out_beam
        print'(A40,f16.8)','absorbed beam energy: ',ext_cross_section
        print'(A40,f16.8)','ext diff energy out: ',energy_out_ext_diff
        print'(A40,f16.8,A2)','beam energy conservation: ',(energy_out_beam+ext_cross_section)/energy_in*100,' %'
        print'(A40,f16.8,A2)','ext diff energy conservation: ',energy_out_ext_diff/energy_in*100,' %'

        print*,'========== end sr energy_checks'
    end if
    ! stop

end subroutine

subroutine beam_loop(   beamV, &
                        beamF1, &
                        beamN, &
                        beamF2, &
                        beamMidpoints, &
                        ampl_beam, &
                        beam_outbeam_tree, &
                        beam_outbeam_tree_counter, &
                        ext_diff_outbeam_tree, &
                        energy_out_beam, &
                        energy_out_ext_diff, &
                        energy_abs_beam, &
                        output_parameters, &
                        job_params, &
                        geometry)

    ! main beam loop

    ! inputs
    integer(8), dimension(:,:), allocatable :: Face1 ! face vertex IDs
    real(8), dimension(:,:), allocatable :: verts ! unique vertices
    real(8) la ! wavelength
    real(8) rbi, ibi ! real part of the refractive index
    ! character(100), intent(in) :: afn ! apertures filename
    real(8), allocatable, dimension(:,:), intent(in) :: beamV ! beam vertices
    real(8), allocatable, dimension(:,:), intent(in) :: beamN ! beam normals
    real(8), allocatable, dimension(:,:), intent(in) :: beamMidpoints ! beam  midpoints
    integer(8), allocatable, dimension(:,:), intent(in) :: beamF1 ! beam face vertex indices
    integer(8), allocatable, dimension(:), intent(in) :: beamF2 ! beam face normal indices
    complex(8), allocatable, dimension(:,:,:), intent(in) :: ampl_beam ! amplitude matrix of incident beam
    integer rec ! max number of internal beam recursions
    type(outbeamtype), dimension(:), allocatable, intent(out) :: beam_outbeam_tree ! outgoing beams from the beam tracing
    type(outbeamtype), dimension(:), allocatable, intent(out) :: ext_diff_outbeam_tree
    integer, intent(out) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
    type(output_parameters_type), intent(inout) :: output_parameters 
    logical is_multithreaded ! whether or not code should use multithreading
    type(job_parameters_type), intent(in) :: job_params
    integer(8), dimension(:), allocatable :: num_face_vert ! number of vertices in each face
    type(geometry_type), intent(in) :: geometry

    real(8), dimension(:,:), allocatable :: Norm ! face normals
    integer, dimension(:), allocatable :: Face2 ! face normal ID of each face
    real(8) start, finish, start1, finish1 ! cpu timing variables
    integer(8) i, j

    ! sr makeAreas
    real(8), dimension(:), allocatable :: faceAreas ! area of each facet
    real(8), dimension(:,:), allocatable :: Midpoints ! face midpoints

    ! sr init
    logical, dimension(:), allocatable :: isVisible ! whether each facet is visible as viewed in -z direction
    logical, dimension(:), allocatable :: isVisiblePlusShadows ! whether each facet is visible, including down-facing facets of illuminated apertures
    logical, dimension(:), allocatable :: isWithinBeam ! whether each visible facet was within the illuminating beam
    logical, dimension(:), allocatable :: isWithinBeam_ps ! whether each visible facet was within the illuminating beam, including down-facing facets of illuminated apertures
    logical, dimension(:), allocatable :: isShadow ! whether the facet was in shadow (down-facing but within the illuminating beam and part of an illuminated aperture)
    real(8), dimension(:), allocatable :: distances ! distances to each illuminated face from the illuminating face ID given by beamIDs
    real(8), dimension(:), allocatable :: distances_ps ! distances to each illuminated face from the illuminating face ID given by beamIDs, including down-facing facets of illuminated apertures
    integer, dimension(:), allocatable :: beamIDs ! the beamF1 ID of the face which illuminated this facet
    integer, dimension(:), allocatable :: beamIDs_ps ! the beamF1 ID of the face which illuminated this facet, including down-facing facets of illuminated aperture
    complex(8), dimension(:,:,:), allocatable :: ampl_in ! amplitude matrix over the surface, for a specific recursion
    complex(8), dimension(:,:,:), allocatable :: ampl_in_ps ! amplitude matrix over the surface, for a specific recursion
    real(8) waveno ! wavenumber in vacuum
    complex(8), dimension(:), allocatable :: rperp ! Fresnel coefficient
    complex(8), dimension(:), allocatable :: rpar ! Fresnel coefficient
    complex(8), dimension(:), allocatable :: tperp ! Fresnel coefficient
    complex(8), dimension(:), allocatable :: tpar ! Fresnel coefficient
    real(8), dimension(:), allocatable :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
    real(8), dimension(:), allocatable :: vk91, vk92, vk93 ! reflected prop vector from each facet
    real(8), dimension(:,:,:), allocatable :: rot_ampl ! rotation matrix for beams incident on eahc facet
    complex(8), dimension(:,:,:), allocatable :: new_in_ampl ! amplitude matrix after rotating into new scattering plane
    complex(8), dimension(:,:,:), allocatable :: new_in_ampl_ps ! amplitude matrix after rotating into new scattering plane
    complex(8), dimension(:,:,:), allocatable :: trans_ampl_ps ! transmitted amplitude matrices, including shadowed facets
    complex(8), dimension(:,:,:), allocatable :: trans_ampl ! transmitted amplitude matrices
    complex(8), dimension(:,:,:), allocatable :: refl_ampl ! reflected amplitude matrices
    complex(8), dimension(:,:,:), allocatable :: refl_ampl_ps ! reflected amplitude matrices
    complex(8), dimension(:,:,:), allocatable :: ampl_diff ! external diffraction amplitude matrices
    integer interactionCounter ! counts the current number of interactions
    
    ! sr findVisibleFacets
    
    ! sr readApertures
    integer(8), dimension(:), allocatable :: apertures ! the aperture which each facet belongs to
    
    ! sr initApertures
    real(8), dimension(:,:), allocatable :: apertureNormals ! the normal of each aperture
    real(8), dimension(:,:), allocatable :: apertureMidpoints ! the midpoint of each aperture
    real(8), dimension(:,:), allocatable :: aperturePropagationVectors ! propagation vector of beam emitted from each aperture
    real(8), dimension(:), allocatable :: apertureAreas ! the total area of each aperture
    real(8), dimension(:), allocatable :: illuminatedApertureAreas ! the llluminated area of each aperture
    logical, dimension(:), allocatable :: sufficientlyIlluminated
    
    ! sr getIlluminatedGeoCrossSection
    real(8) illuminatedGeoCrossSection ! the total illuminated geometric cross section
    
    ! sr getThreshold
    real(8) threshold ! minimum area of illumination per aperture to create new beam
    
    ! beam_loop -> diffraction
    
    ! main loop
    integer, dimension(:,:), allocatable :: F1Mapping ! the face indices of each row of variables to go into main loop
    complex(8), dimension(:,:), allocatable :: refl_ampl_out11Int ! the amplitude matrix that goes into the main loop
    complex(8), dimension(:,:), allocatable :: refl_ampl_out12Int ! needs renaming
    complex(8), dimension(:,:), allocatable :: refl_ampl_out21Int
    complex(8), dimension(:,:), allocatable :: refl_ampl_out22Int
    real(8), dimension(:,:,:), allocatable :: propagationVectors ! propagation vector of beam emitted from each aperture
    real(8), dimension(:,:), allocatable :: vk71Int, vk72Int, vk73Int ! reflected e-perp vector
    
    ! energy_checks
    real(8), intent(out) :: energy_out_beam
    real(8), intent(out) :: energy_abs_beam
    real(8), intent(out) :: energy_out_ext_diff
    
    ! storage checks
    real(8) total_storage

    if(job_params%timing) then
        start = omp_get_wtime()
    endif

    if(job_params%debug >= 1) then
        print*,'start beam loop...'
    end if

    la = job_params%la
    rbi = job_params%rbi
    ibi = job_params%ibi
    rec = job_params%rec
    is_multithreaded = job_params%is_multithreaded
    ! is_multithreaded = .false.
    ! bodge conversion while i refactor
    allocate(verts(1:geometry%nv,1:3))
    verts(:,:) = geometry%v(:,:)
    allocate(num_face_vert(1:geometry%nf))
    num_face_vert(:) = geometry%f(:)%nv
    allocate(apertures(1:geometry%nf))
    apertures(:) = geometry%f(:)%ap
    allocate(Face1(1:geometry%nf,1:maxval(num_face_vert)))
    do i = 1, geometry%nf
        do j = 1, num_face_vert(i)
            Face1(i,j) = geometry%f(i)%vi(j)
        end do
    end do

    call make_normals(Face1, verts, Face2, Norm) ! recalculate normals

    call midPointsAndAreas(Face1, verts, Midpoints, faceAreas, num_face_vert) ! calculate particle facet areas
    
    call init(  Face1, isVisible, isVisiblePlusShadows, isWithinBeam, distances, beamIDs, &
                isWithinBeam_ps, distances_ps, beamIDs_ps, isShadow, ampl_in, ampl_in_ps, la, waveno, &
                rperp, rpar, tperp, tpar, vk71, vk72, vk73, vk91, vk92, vk93, rot_ampl, new_in_ampl, &
                new_in_ampl_ps, trans_ampl_ps, trans_ampl, refl_ampl, refl_ampl_ps, &
                beam_outbeam_tree, beam_outbeam_tree_counter, interactionCounter) ! initialise some variables
    ! move beam_loop_mod variables into separate beam_loop_init subroutine and out of main program later...

    ! ############# beam_loop_mod #############
    
    ! call readApertures(afn, apertures, Face1) ! read aperture assignments from file (now moved to pdal2 and passed into beam_loop)
    
    call initApertures(apertures, Norm, Face2, Midpoints, apertureNormals, apertureMidpoints, apertureAreas, illuminatedApertureAreas, sufficientlyIlluminated, &
                             aperturePropagationVectors) ! initialise the aperture variables
    
    call findVisibleFacets(verts, Face1, Norm, Face2, midPoints, isVisible, isVisiblePlusShadows, apertures, num_face_vert) ! find external visible facets in this orientation
    ! stop
    call findWithinBeam(Face1, Midpoints, isVisible, beamV, beamF1, beamN, beamF2, beamMidpoints, isWithinBeam, distances, beamIDs) ! find visible facets within the incident beam
    
    call findWithinBeam(Face1, Midpoints, isVisiblePlusShadows, beamV, beamF1, beamN, beamF2, beamMidpoints, isWithinBeam_ps, distances_ps, beamIDs_ps) ! find visible facets (including shadow) within the incident beam
    
    call getShadow(isWithinBeam,isWithinBeam_ps,isShadow)
    
    call getApertureAreas(apertures, isWithinBeam, apertureAreas, illuminatedApertureAreas, apertureNormals, faceAreas)
    
    call getIlluminatedGeoCrossSection(faceAreas, isWithinBeam, Norm, Face2, illuminatedGeoCrossSection) ! get illuminated geometric cross section using particle in current orientation

    call getThreshold(la,threshold,job_params) ! get minimum illuminatino required to create a new beam
    
    call getSufficientlyIlluminated(illuminatedApertureAreas,threshold,sufficientlyIlluminated) ! get logical array containing which apertures were sufficiently illuminated
    
    call propagateBeam(ampl_in, ampl_in_ps, ampl_beam, isWithinBeam, isWithinBeam_ps, beamIDs, beamIDs_ps, waveno, distances, distances_ps, Face1) ! use distances to propagate fields to the points of intersection
    
    call getFresnel(Face1, isShadow, Norm, Face2, apertureNormals, apertures, rbi, rperp, rpar, tperp, tpar, ibi) ! calculate fresnel coefficients based on facet normals (or aperture normals for shadow facets)
    
    call getPropagationVectors(aperturePropagationVectors,sufficientlyIlluminated,apertureNormals,rbi) ! get propagation vectors based on aperture normals
    
    call getReflectionVectors(Norm, Face2, isShadow, apertureNormals, apertures, &
                                    vk71, vk72, vk73, vk91, vk92, vk93) ! get reflection vectors based on facet normals (WILL NOT WORK FOR RE-ENTRY)
    
    call getRotationMatrices(rot_ampl, vk71, vk72, vk73, 1D0, 0D0, 0D0, 0D0, 0D0, -1D0) ! get rotation matrix for rotating into new scattering plane

    call rotateAmplitudeMatrices(rot_ampl,ampl_in,new_in_ampl) ! rotate amplitude matrices into new scattering plane

    call rotateAmplitudeMatrices(rot_ampl,ampl_in_ps,new_in_ampl_ps) ! rotate amplitude matrices, including shadowed facets, into new scattering plane

    call applyFresnelMatrices(new_in_ampl, rperp, rpar, tperp, tpar, refl_ampl, trans_ampl) ! apply fresnel matrices
    
    call applyFresnelMatrices(new_in_ampl_ps, rperp, rpar, tperp, tpar, refl_ampl_ps, trans_ampl_ps) ! apply fresnel matrices, including shadow facets
    
    call make_external_diff_outbeam_tree(new_in_ampl, ampl_diff, isVisible, ext_diff_outbeam_tree, vk71, vk72, vk73, verts, Face1) ! save amplitude matrix after rotation into scattering plane for later
    
    call add_refl_to_outbeam_tree(beam_outbeam_tree, beam_outbeam_tree_counter, isWithinBeam, &
                                        refl_ampl, vk91, vk92, vk93, vk71, vk72, vk73, apertures, interactionCounter)
    
    ! this sr needs optimising (currently has facets with 0 ampl. matrix being passed to main loop)
    call init_for_main_loop(      Face1, trans_ampl_ps, F1Mapping, &
                                  refl_ampl_out11Int, refl_ampl_out12Int, refl_ampl_out21Int, refl_ampl_out22Int, &
                                  aperturePropagationVectors, propagationVectors, &
                                  vk71, vk72, vk73, vk71Int, vk72Int, vk73Int)
    
    ! main loop
    do i = 1, rec ! for each beam recursion
        if(job_params%debug >= 2) then
            print*,'internal recursion #',i
        end if
        if(job_params%timing) then
            if(job_params%debug >= 1) then
                start1 = omp_get_wtime()
            end if
        end if
       call internal_recursion_outer(   F1Mapping, &
                                        propagationVectors, &
                                        verts, Norm, midPoints, Face1, Face2, &
                                        apertureMidpoints, apertureNormals, apertures, &
                                        faceAreas, &
                                        threshold, &
                                        rbi, ibi, waveno, &
                                        vk71Int, vk72Int, vk73Int, &
                                        beam_outbeam_tree, beam_outbeam_tree_counter, &
                                        refl_ampl_out11Int, refl_ampl_out12Int, refl_ampl_out21Int, refl_ampl_out22Int, &
                                        interactionCounter, &
                                        is_multithreaded)
        if(job_params%timing) then
            if(job_params%debug >= 2) then
                finish1 = omp_get_wtime()
                print'(A,I3,A,f16.8,A)',"recursion",i," - time elapsed: ",finish1-start1," secs"
                write(101,'(A,I3,A,f16.8,A)')"recursion",i," - time elapsed: ",finish1-start1," secs"
            end if
        end if
    end do

    if(job_params%debug >= 1) then
        if(job_params%timing) then
            finish = omp_get_wtime()
            print*,'=========='
            print'(A,f16.8,A)',"end beam loop - time elapsed: ",finish-start," secs"
            ! print*,'=========='
            write(101,*)'=========='
            write(101,'(A,f16.8,A)')"end beam loop - time elapsed: ",finish-start," secs"
            ! write(101,*)'=========='    
        end if
    end if

    call get_beamtree_vert(beam_outbeam_tree, beam_outbeam_tree_counter, verts, Face1)    

    call energy_checks( beam_outbeam_tree, &
                        beam_outbeam_tree_counter, &
                        norm, &
                        face2, &
                        faceAreas, &
                        energy_out_beam, &
                        illuminatedGeoCrossSection, &
                        ext_diff_outbeam_tree, &
                        energy_out_ext_diff, &
                        energy_abs_beam, &
                        job_params)

    output_parameters%geo_cross_sec = illuminatedGeoCrossSection

    ! get storage size:
    ! double precision is 8 bytes (64 bits)
    if(job_params%debug >= 1) then
        total_storage = real((  sizeof(verts) + &
                                sizeof(F1Mapping) + &
                                sizeof(propagationVectors) + &
                                sizeof(Norm) + &
                                sizeof(midPoints) + &
                                sizeof(faceAreas) + &
                                sizeof(vk71Int) + &
                                sizeof(beam_outbeam_tree) + &
                                sizeof(refl_ampl_out11Int) + &
                                sizeof(ext_diff_outbeam_tree)))/1048576D0
        print'(A,f10.2,A)','beam loop estimated memory usage (per mpi process):',total_storage,' MB'
    end if

    if(job_params%debug >= 2) then
        print'(A)','memory usage breakdown (per mpi process):'
        print'(A,f10.2,A)','particle vertices: ',real(sizeof(verts))/1048576D0,' MB'
        print'(A,f10.2,A)','particle illuminated faces: ',real(sizeof(F1Mapping))/1048576D0,' MB'
        print'(A,f10.2,A)','beam propagation vecs: ',real(sizeof(propagationVectors))/1048576D0,' MB'
        print'(A,f10.2,A)','particle normals: ',real(sizeof(Norm))/1048576D0,' MB'
        print'(A,f10.2,A)','particle midpoints: ',real(sizeof(midPoints))/1048576D0,' MB'
        print'(A,f10.2,A)','particle face areas: ',real(sizeof(faceAreas))/1048576D0,' MB'
        print'(A,f10.2,A)','perp. field vecs: ',real(sizeof(vk71Int))/1048576D0,' MB'
        print'(A,f10.2,A)','outbeam tree: ',real(sizeof(beam_outbeam_tree))/1048576D0,' MB'
        print'(A,f10.2,A)','amplitude matrices: ',real(sizeof(refl_ampl_out11Int))/1048576D0,' MB'
        print'(A,f10.2,A)','external diff beamtree: ',real(sizeof(ext_diff_outbeam_tree))/1048576D0,' MB'
        print'(A)','note: this is an underestimate of the total memory usage.'
        print'(A)',' =========='
    end if

    ! stop

end subroutine

subroutine get_beamtree_vert(beam_outbeam_tree, beam_outbeam_tree_counter, verts, Face1)

    ! uses the FOut data to add vertex information to a beamtree

type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
integer, intent(in) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs

integer i, face_id

do i = 1, beam_outbeam_tree_counter ! for each beam
    face_id = beam_outbeam_tree(i)%FOut ! get the face id
    beam_outbeam_tree(i)%verts = transpose(verts(Face1(face_id,1:3),1:3)) ! retrieve vertex information
end do

end subroutine

subroutine find_vis(rotatedVert, Face1, Face2, &
                    apertures, &
                    rotatedMidpoints, rotatedNorm, &
                    in_beam, dist_beam, id_beam, is_shad, &
                    rotatedapertureNormals, &
                    is_multithreaded, &
                    beam)

    ! for subroutine beam_scan
    ! for a given particle orientation with assumed propagation along the -z axis,
    ! computes the internally illuminated facets for a given illuminating aperture
    ! uses a z-buffer technique to increase speed, among a few other tricks

integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), dimension(:,:), allocatable, intent(in) :: rotatedVert
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
real(8), dimension(:,:), allocatable, intent(in) :: rotatedMidpoints
real(8), dimension(:,:), allocatable, intent(in) :: rotatedNorm
logical, dimension(:), allocatable, intent(out) :: in_beam
real(8), dimension(:), allocatable, intent(out) :: dist_beam
integer, dimension(:), allocatable, intent(out) :: id_beam
logical, dimension(:), allocatable, intent(out) :: is_shad
real(8), dimension(:,:), allocatable, intent(in) :: rotatedapertureNormals
logical, intent(in) :: is_multithreaded
type(beam_type), intent(in) :: beam

logical am_i_multithreaded

integer i, k, m
integer j
logical, dimension(:), allocatable :: is_beam
logical, dimension(:), allocatable :: is_vis
real(8), allocatable, dimension(:,:) :: boundingBoxV
integer(8), allocatable, dimension(:,:) :: boundingBoxF
integer(8) boundingBoxFSize
real(8), dimension(:,:), allocatable :: boundingBoxMidpoints ! unique vertices, face vertex IDs, face normals, face midpoints
real(8), dimension(:), allocatable :: boundingBoxFaceAreas
integer(8), dimension(:), allocatable :: boundingBoxNumFaceVert
integer, dimension(:), allocatable :: F3 ! bounding box IDs
integer, dimension(:,:), allocatable :: F4 ! fuzzy bounding box IDs
real(8), dimension(:), allocatable :: distanceToBB, distanceToFuzzyBB
integer BB
logical within_bounds
real(8) vecb1, vecb2
real(8) edge_norm1, edge_norm2
real(8) edge_check
real(8) beamXmax0, beamXmin0, beamYmax0, beamYmin0
real(8) beamXmax, beamXmin, beamYmax, beamYmin
logical, dimension(:,:), allocatable :: F5
real(8) start, finish
integer, dimension(:), allocatable :: mapping

! ################################
! start new ray tracing algorithm

call CPU_TIME(start)

! use the current crystal vertices to create some bounding boxes in x-y plane
call beam_aligned_bounding_boxes(rotatedVert, boundingBoxV, boundingBoxF)

allocate(boundingBoxNumFaceVert(1:size(boundingBoxF,1)))
boundingBoxNumFaceVert = 4

! compute bounding box midpoints
call midPointsAndAreas(boundingBoxF, boundingBoxV, boundingBoxMidpoints, boundingBoxFaceAreas, boundingBoxNumFaceVert)

allocate(F3(1:size(Face1,1))) ! array to hold index of bounding box that each face belongs to
allocate(F4(1:size(Face1,1),1:3)) ! array to hold index of fuzzy bounding box that each face belongs to
boundingBoxFSize = size(boundingBoxF,1) ! number of bounding box faces
allocate(distanceToBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box
allocate(distanceToFuzzyBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box

am_i_multithreaded = is_multithreaded
am_i_multithreaded = .false. ! disable multithreading for now...

! find which bounding box each vertex belongs to
if(am_i_multithreaded) then
    !$OMP PARALLEL DEFAULT(SHARED) num_threads(omp_get_max_threads()) PRIVATE(i,j,distanceToBB,distanceToFuzzyBB)
    !$OMP DO
    do i = 1, size(Face1,1) ! for each face
        distanceToBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - rotatedMidpoints(i,1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - rotatedMidpoints(i,2))**2) ! distance to each bb
        F3(i) = minloc(distanceToBB,1) ! record which bounding box midpoint this facet was closest to
        do j = 1, 3 ! for each vertex
            distanceToFuzzyBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - rotatedVert(Face1(i,j),1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - rotatedVert(Face1(i,j),2))**2) ! distance to each bb
            F4(i,j) = minloc(distanceToFuzzyBB,1) ! record which bounding box midpoint this facet vertex was closest to
        end do 
    end do
    !$OMP END PARALLEL
else
    do i = 1, size(Face1,1) ! for each face
        distanceToBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - rotatedMidpoints(i,1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - rotatedMidpoints(i,2))**2) ! distance to each bb
        F3(i) = minloc(distanceToBB,1) ! record which bounding box midpoint this facet was closest to
        do j = 1, 3 ! for each vertex
            distanceToFuzzyBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - rotatedVert(Face1(i,j),1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - rotatedVert(Face1(i,j),2))**2) ! distance to each bb
            F4(i,j) = minloc(distanceToFuzzyBB,1) ! record which bounding box midpoint this facet vertex was closest to
        end do 
    end do
end if

! print*,'F3(500)',F3(500)
! print*,'F4(500,3)',F4(500,3)

! stop

! more optimisation
allocate(F5(1:size(Face2,1),1:boundingBoxFSize))
F5 = .false. ! init
do i = 1, size(Face2,1)
    do j = 1, boundingBoxFSize
        if(F4(i,1) .eq. j .or. F4(i,2) .eq. j .or. F4(i,3) .eq. j) F5(i,j) = .true.
    end do
end do

! allocate and init
allocate(is_vis(1:size(Face2,1)))
allocate(is_shad(1:size(Face2,1)))
allocate(in_beam(1:size(Face2,1)))
allocate(is_beam(1:size(Face2,1)))
allocate(id_beam(1:size(Face2,1)))
allocate(dist_beam(1:size(Face2,1)))
allocate(mapping(1:size(Face2,1)))

is_vis = .true.
is_shad = .false.
in_beam = .false.
is_beam = .false.
id_beam = 0
dist_beam = 0
mapping = 0

! record which facets are part of the illuminating beam
do i = 1, beam%num_facets
    is_beam(beam%field(i)%face_id) = .true.
    mapping(beam%field(i)%face_id) = i ! map each illuminating facet to its position in the beam data structure
end do

! need to add a check to ignore facets outside beam max/min here
! ok well this gave 5x speedup
do i = 1, 3 ! bodge
    beamXmax0 = maxval(rotatedVert(Face1(beam%field(1:beam%num_facets)%face_id,i),1))
    beamXmin0 = minval(rotatedVert(Face1(beam%field(1:beam%num_facets)%face_id,i),1))
    beamYmax0 = maxval(rotatedVert(Face1(beam%field(1:beam%num_facets)%face_id,i),2))
    beamYmin0 = minval(rotatedVert(Face1(beam%field(1:beam%num_facets)%face_id,i),2))
    if(i .eq. 1) then
        beamXmax = beamXmax0
        beamXmin = beamXmin0
        beamYmax = beamYmax0
        beamYmin = beamYmin0
    else
        if(beamXmax0 .gt. beamXmax) beamXmax = beamXmax0
        if(beamXmin0 .lt. beamXmin) beamXmin = beamXmin0
        if(beamYmax0 .gt. beamYmax) beamYmax = beamYmax0
        if(beamYmin0 .lt. beamYmin) beamYmin = beamYmin0        
    end if
end do
do i = 1, size(Face2,1)
    if(rotatedMidpoints(i,1) .lt. beamXmin) then
        is_vis(i) = .false.
    else if(rotatedMidpoints(i,1) .gt. beamXmax) then
        is_vis(i) = .false.
    else if(rotatedMidpoints(i,2) .lt. beamYmin) then
        is_vis(i) = .false.
    else if(rotatedMidpoints(i,2) .gt. beamYmax) then
        is_vis(i) = .false.
    end if
end do

if(am_i_multithreaded) then
    !$OMP PARALLEL DEFAULT(SHARED) num_threads(omp_get_max_threads()) PRIVATE(j,k,m,BB,within_bounds,edge_norm1,edge_norm2,vecb1,vecb2,edge_check)
    !$OMP DO
    do m = 1, size(Face2,1) ! for each facet m
        if(is_vis(m) .eqv. .false.) then ! if facet isnt visible
            ! do nothing
        else
            if(rotatedapertureNormals(apertures(m),3) .gt. -0.01) then ! if aperture is downfacing
                is_vis(m) = .false. ! set not visible
            else ! if aperture was facing towards incidence
                BB = F3(m) ! get bounding box ID
                do j = 1, size(Face2,1) ! for each potentially blocking facet j
                    ! if(rotatedNorm(Face2(j),3) .lt. 0.01) then ! sign flip here for internal
                    ! if(F4(j,1) .eq. BB .or. F4(j,2) .eq. BB .or. F4(j,3) .eq. BB) then ! 
                    if(F5(j,BB)) then ! 
                    ! if(any(F4(j,1:3)) .eq. BB) then ! if blocker was in fuzzy bounding box
                        if(j .ne. m) then ! ignore self-block
                            if(rotatedNorm(Face2(j),3) .lt. 0.01) then ! if down-facing, sign flip here for internal
                                ! do nothing
                            else ! if up-facing
                                if (rotatedMidpoints(m,3) .gt. rotatedMidpoints(j,3)) then ! if potential blocker was behind facet m
                                    ! do nothing
                                else ! if potential blocker was in front of facet m
                                    ! do bounded surface edge check
                                    within_bounds = .true. ! assume centroid of facet m is within the bounded surface of potentially blocking facet j
                                    do k = 1, 3 ! looping over vertices of potentially blocking facet j
                                        ! compute edge normal
                                        if(k .eq. 3) then
                                            edge_norm2 = -rotatedVert(Face1(j,1),1) + rotatedVert(Face1(j,3),1) ! cross product of edge vector with reverse beam direction
                                            edge_norm1 = rotatedVert(Face1(j,1),2) - rotatedVert(Face1(j,3),2)
                                        else
                                            edge_norm2 = -rotatedVert(Face1(j,k+1),1) + rotatedVert(Face1(j,k),1)
                                            edge_norm1 = rotatedVert(Face1(j,k+1),2) - rotatedVert(Face1(j,k),2)
                                        end if
                                        vecb1 = rotatedMidpoints(m,1) - rotatedVert(Face1(j,k),1) ! vector from vertex k of potential blocker j to centroid of facet m
                                        vecb2 = rotatedMidpoints(m,2) - rotatedVert(Face1(j,k),2)
                                        ! edge_norm1 = edge_vec2
                                        ! edge_norm2 = -edge_vec1
                                        ! nf = sqrt(edge_norm1**2 + edge_norm2**2) ! probably dont even need this step
                                        ! edge_norm1 = edge_norm1 / nf ! can confirm dont need this, saves 16% time of beam loop
                                        ! edge_norm2 = edge_norm2 / nf
                                        edge_check = vecb1*edge_norm1 + vecb2*edge_norm2 ! dot product of edge vector with vetor B
                                        if(edge_check .gt. 0) within_bounds = .false. ! if edge check fails, centroid of facet m is not within the bounded surface of facet j
                                    end do
                                    if(within_bounds .eqv. .false.) then ! if facet m is not within bounded surface of facet j
                                        !do nothing
                                    else
                                        if(is_beam(j)) then ! if facet j was part of the illuminating surface
                                            if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                ! check to see which is the closest
                                                if(rotatedMidpoints(j,3) .gt. rotatedMidpoints(id_beam(m),3)) then ! if facet j was closer than the previous blocker
                                                    in_beam(m) = .true. ! set to be within beam
                                                    id_beam(m) = j ! record the blocking facet
                                                    dist_beam(m) = rotatedMidpoints(j,3) - rotatedMidpoints(m,3) ! record the distance from centroid of blocker to centroid of facet m
                                                else
                                                    ! do nothing                                     
                                                end if
                                            else ! if this is the first time finding a blocking beam facet
                                                in_beam(m) = .true. ! set to be within beam
                                                ! id_beam(m) = j ! record the blocking facet
                                                id_beam(m) = mapping(j) 
                                                dist_beam(m) = rotatedMidpoints(j,3) - rotatedMidpoints(m,3) ! record the distance from centroid of blocker to centroid of facet m
                                            end if
                                        else
                                            if(apertures(m) .eq. apertures(j)) then ! if facet j and facet m belong to the same aperture
                                                is_shad(m) = .true. ! set m as a shadow facet and continue to search Dfor blockers
                                            else ! if facet j and facet m dont belong to the same aperture
                                                if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                    if(rotatedMidpoints(j,3) .gt. rotatedMidpoints(id_beam(m),3)) then ! if facet j was behind than the blocking beam
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
    !$OMP END PARALLEL
else
    do m = 1, size(Face2,1) ! for each facet m
        if(is_vis(m) .eqv. .false.) then ! if facet isnt visible
            ! do nothing
        else
            if(rotatedapertureNormals(apertures(m),3) .gt. -0.01) then ! if aperture is downfacing
                is_vis(m) = .false. ! set not visible
            else ! if aperture was facing towards incidence
                BB = F3(m) ! get bounding box ID
                do j = 1, size(Face2,1) ! for each potentially blocking facet j
                    ! if(rotatedNorm(Face2(j),3) .lt. 0.01) then ! sign flip here for internal
                    ! if(F4(j,1) .eq. BB .or. F4(j,2) .eq. BB .or. F4(j,3) .eq. BB) then ! 
                    if(F5(j,BB)) then ! 
                    ! if(any(F4(j,1:3)) .eq. BB) then ! if blocker was in fuzzy bounding box
                        if(j .ne. m) then ! ignore self-block
                            if(rotatedNorm(Face2(j),3) .lt. 0.01) then ! if down-facing, sign flip here for internal
                                ! do nothing
                            else ! if up-facing
                                if (rotatedMidpoints(m,3) .gt. rotatedMidpoints(j,3)) then ! if potential blocker was behind facet m
                                    ! do nothing
                                else ! if potential blocker was in front of facet m
                                    ! do bounded surface edge check
                                    within_bounds = .true. ! assume centroid of facet m is within the bounded surface of potentially blocking facet j
                                    do k = 1, 3 ! looping over vertices of potentially blocking facet j
                                        ! compute edge normal
                                        if(k .eq. 3) then
                                            edge_norm2 = -rotatedVert(Face1(j,1),1) + rotatedVert(Face1(j,3),1) ! cross product of edge vector with reverse beam direction
                                            edge_norm1 = rotatedVert(Face1(j,1),2) - rotatedVert(Face1(j,3),2)
                                        else
                                            edge_norm2 = -rotatedVert(Face1(j,k+1),1) + rotatedVert(Face1(j,k),1)
                                            edge_norm1 = rotatedVert(Face1(j,k+1),2) - rotatedVert(Face1(j,k),2)
                                        end if
                                        vecb1 = rotatedMidpoints(m,1) - rotatedVert(Face1(j,k),1) ! vector from vertex k of potential blocker j to centroid of facet m
                                        vecb2 = rotatedMidpoints(m,2) - rotatedVert(Face1(j,k),2)
                                        ! edge_norm1 = edge_vec2
                                        ! edge_norm2 = -edge_vec1
                                        ! nf = sqrt(edge_norm1**2 + edge_norm2**2) ! probably dont even need this step
                                        ! edge_norm1 = edge_norm1 / nf ! can confirm dont need this, saves 16% time of beam loop
                                        ! edge_norm2 = edge_norm2 / nf
                                        edge_check = vecb1*edge_norm1 + vecb2*edge_norm2 ! dot product of edge vector with vetor B
                                        if(edge_check .gt. 0) within_bounds = .false. ! if edge check fails, centroid of facet m is not within the bounded surface of facet j
                                    end do
                                    if(within_bounds .eqv. .false.) then ! if facet m is not within bounded surface of facet j
                                        !do nothing
                                    else
                                        if(is_beam(j)) then ! if facet j was part of the illuminating surface
                                            if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                ! check to see which is the closest
                                                if(rotatedMidpoints(j,3) .gt. rotatedMidpoints(id_beam(m),3)) then ! if facet j was closer than the previous blocker
                                                    in_beam(m) = .true. ! set to be within beam
                                                    id_beam(m) = j ! record the blocking facet
                                                    dist_beam(m) = rotatedMidpoints(j,3) - rotatedMidpoints(m,3) ! record the distance from centroid of blocker to centroid of facet m
                                                else
                                                    ! do nothing                                     
                                                end if
                                            else ! if this is the first time finding a blocking beam facet
                                                in_beam(m) = .true. ! set to be within beam
                                                id_beam(m) = mapping(j) ! record the position of the blocking facet in the beam data struct
                                                dist_beam(m) = rotatedMidpoints(j,3) - rotatedMidpoints(m,3) ! record the distance from centroid of blocker to centroid of facet m
                                            end if
                                        else
                                            if(apertures(m) .eq. apertures(j)) then ! if facet j and facet m belong to the same aperture
                                                is_shad(m) = .true. ! set m as a shadow facet and continue to search Dfor blockers
                                            else ! if facet j and facet m dont belong to the same aperture
                                                if(in_beam(m)) then ! if a blocking beam facet had already been found
                                                    if(rotatedMidpoints(j,3) .gt. rotatedMidpoints(id_beam(m),3)) then ! if facet j was behind than the blocking beam
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
end if

call CPU_TIME(finish)

end subroutine

subroutine rotate_back_propagation_vectors(rotationMatrices, propagationVectors2, sufficientlyIlluminated2, aperture_id, counter)

! rotates propagation vectors back into original coordinate system

real(8), dimension(:,:,:), allocatable, intent(in) :: rotationMatrices
real(8), dimension(:,:,:), allocatable, intent(inout) :: propagationVectors2
logical, dimension(:), allocatable, intent(in) :: sufficientlyIlluminated2 ! whether each aperture was sufficiently illuminated by the previous beam
integer, intent(in) :: counter
integer, intent(in) :: aperture_id

real(8), dimension(:,:,:), allocatable :: propagationVectors2_temp
integer i

! allocate and save a temporary variable to avoid risk of overwriting
allocate(propagationVectors2_temp(1:size(propagationVectors2,1),1:size(propagationVectors2,2),1:size(propagationVectors2,3)))
propagationVectors2_temp(1:size(propagationVectors2,1),1:size(propagationVectors2,2),1:size(propagationVectors2,3)) = & 
    propagationVectors2(1:size(propagationVectors2,1),1:size(propagationVectors2,2),1:size(propagationVectors2,3))

do i = 1, size(sufficientlyIlluminated2,1) ! for each illuminating aperture
    if(sufficientlyIlluminated2(i)) then
        propagationVectors2(i,1:3,counter) = matmul(transpose(rotationMatrices(1:3,1:3,aperture_id)),propagationVectors2_temp(i,1:3,counter))
    end if
end do


end subroutine

subroutine get_reflection_vectors(propagationVectors2, rotatedapertureNormals, sufficientlyIlluminated2, counter)

integer, intent(in) :: counter
real(8), dimension(:,:,:), allocatable, intent(inout) :: propagationVectors2
real(8), dimension(:,:), allocatable, intent(in) :: rotatedapertureNormals
logical, dimension(:), allocatable, intent(in) :: sufficientlyIlluminated2 ! whether each aperture was sufficiently illuminated by the previous beam

integer i

do i = 1, size(sufficientlyIlluminated2,1)
    if(sufficientlyIlluminated2(i)) then ! for each sufficiently illuminated aperture
        propagationVectors2(i,1,counter) =  0 + 2*rotatedapertureNormals(i,3)*rotatedapertureNormals(i,1)
        propagationVectors2(i,2,counter) =  0 + 2*rotatedapertureNormals(i,3)*rotatedapertureNormals(i,2)
        propagationVectors2(i,3,counter) = -1 + 2*rotatedapertureNormals(i,3)*rotatedapertureNormals(i,3)
    end if
end do


end subroutine

subroutine rotate_into_aperture_system( aperturePropagationVectors, apertureMidpoints, apertureNormals, &
                                        verts, Norm, midPoints, &
                                        rotatedapertureMidpoints, &
                                        rotatedapertureNormals, &
                                        rotatedaperturePropagationVectors, &
                                        rotatedVert, &
                                        rotatedNorm, &
                                        rotatedMidpoints, &
                                        rotationMatrices, &
                                        beam)

real(8), dimension(:,:), allocatable, intent(in) :: aperturePropagationVectors
real(8), dimension(:,:), allocatable, intent(in) :: apertureMidpoints ! the midpoint of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints ! face midpoints
real(8), dimension(:,:), allocatable, intent(out) :: rotatedaperturePropagationVectors
real(8), dimension(:,:), allocatable, intent(out) :: rotatedapertureMidpoints
real(8), dimension(:,:), allocatable, intent(out) :: rotatedapertureNormals
real(8), dimension(:,:), allocatable, intent(out) :: rotatedVert
real(8), dimension(:,:), allocatable, intent(out) :: rotatedNorm
real(8), dimension(:,:), allocatable, intent(out) :: rotatedMidpoints
real(8), dimension(:,:,:), allocatable, intent(inout) :: rotationMatrices
type(beam_type), intent(in) :: beam

real(8) theta_1, theta_2
real(8) rot1(1:3,1:3), rot2(1:3,1:3), rot(1:3,1:3)
real(8) temp_vector(1:3)
integer i, j

! print*,'rotating into aperture system for aperture:',aperture_id
! print*,'propagation direction is',aperturePropagationVectors(aperture_id,1:3)

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

do i = 1, 3
    do j = 1, 3
        rotationMatrices(i,j,beam%aperture_id) = rot(i,j)
    end do
end do

! allocations (there might be a more elegant way to do this later)
allocate(rotatedaperturePropagationVectors(1:size(aperturePropagationVectors,1),1:size(aperturePropagationVectors,2)))
allocate(rotatedapertureMidpoints(1:size(apertureMidpoints,1),1:size(apertureMidpoints,2)))
allocate(rotatedapertureNormals(1:size(apertureNormals,1),1:size(apertureNormals,2)))
allocate(rotatedVert(1:size(verts,1),1:size(verts,2)))
allocate(rotatedNorm(1:size(Norm,1),1:size(Norm,2)))
allocate(rotatedMidpoints(1:size(Midpoints,1),1:size(Midpoints,2)))

do j = 1, size(verts,1) ! rotate vertices
    rotatedVert(j,1:3) = matmul(rot,verts(j,1:3))
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated vertices: ',rotatedVert(j,1), rotatedVert(j,2), rotatedVert(j,3)
end do
do j = 1, size(Norm,1) ! rotate normals
    rotatedNorm(j,1:3) = matmul(rot,Norm(j,1:3))
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated normals: ',rotatedNorm(j,1), rotatedNorm(j,2), rotatedNorm(j,3)
end do
do j = 1, size(Midpoints,1) ! rotate midpoints
    rotatedMidpoints(j,1:3) = matmul(rot,Midpoints(j,1:3))
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated midpoints: ',rotatedMidpoints(j,1), rotatedMidpoints(j,2), rotatedMidpoints(j,3)
end do
do j = 1, size(aperturePropagationVectors,1) ! rotate aperture variables
    rotatedapertureMidpoints(j,1:3) = matmul(rot,apertureMidpoints(j,1:3))
    rotatedapertureNormals(j,1:3) = matmul(rot,apertureNormals(j,1:3))
    rotatedaperturePropagationVectors(j,1:3) = matmul(rot,aperturePropagationVectors(j,1:3))
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated aperture midpoints: ',rotatedapertureMidpoints(j,1), rotatedapertureMidpoints(j,2), rotatedapertureMidpoints(j,3)
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated aperture normals: ',rotatedapertureNormals(j,1), rotatedapertureNormals(j,2), rotatedapertureNormals(j,3)
    ! print'(A,I,A,f12.6,f12.6,f12.6)','j',j,' rotated aperture propagation: ',rotatedaperturePropagationVectors(j,1), rotatedaperturePropagationVectors(j,2), rotatedaperturePropagationVectors(j,3)
end do

! some renormalisation might be required here

! print*,'temp_vector',temp_vector
! print*,'theta_2',theta_2

! do i = 1, 3
!     do j = 1,3
!         ! print'(A,I,A,I,A,f12.6)','rot1(',i,',',j,')',rot1(i,j)
!         ! print'(A,I,A,I,A,f12.6)','rot2(',i,',',j,')',rot2(i,j)
!         print'(A,I,A,I,A,f12.6)','rot(',i,',',j,')',rot(i,j)
!     end do
! end do


end subroutine

subroutine beam_scan(   aperturePropagationVectors, apertureMidpoints, apertureNormals, &
                        verts, Norm, midPoints, Face1, Face2, &
                        apertures, &
                        faceAreas, &
                        threshold, &
                        sufficientlyIlluminated2, &
                        in_beam, &
                        rotatedapertureNormals, &
                        rotationMatrices, &
                        id_beam, dist_beam, &
                        is_shad, &
                        is_multithreaded, &
                        beam)

real(8), dimension(:,:), allocatable, intent(in) :: aperturePropagationVectors
real(8), dimension(:,:), allocatable, intent(in) :: apertureMidpoints ! the midpoint of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints ! face midpoints
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), dimension(:), allocatable, intent(in) :: faceAreas ! area of each facet
real(8), intent(in) :: threshold ! minimum area of illumination per aperture to create new beam
logical, dimension(:), allocatable, intent(inout) :: sufficientlyIlluminated2 ! whether each aperture was sufficiently illuminated by the previous beam
real(8), dimension(:,:), allocatable, intent(out) :: rotatedapertureNormals
real(8), dimension(:,:,:), allocatable, intent(inout) :: rotationMatrices
real(8), dimension(:,:), allocatable :: rotatedVert
real(8), dimension(:,:), allocatable :: rotatedNorm
! new ray tracing algorithm
logical, dimension(:), allocatable, intent(out) :: is_shad
logical, dimension(:), allocatable, intent(out) :: in_beam
integer, dimension(:), allocatable, intent(out) :: id_beam
real(8), dimension(:), allocatable :: dist_beam
logical, intent(in) :: is_multithreaded
type(beam_type), intent(in) :: beam

integer k

! rotated variables
real(8), dimension(:,:), allocatable :: rotatedaperturePropagationVectors
real(8), dimension(:,:), allocatable :: rotatedapertureMidpoints
real(8), dimension(:,:), allocatable :: rotatedMidpoints
! beam variables
! findVisible and findWithinBeam
integer i, counter2
! real(8), dimension(:), allocatable :: distances
! integer, dimension(:), allocatable :: beamIDs
! illuminations
integer totalIlluminated_ps
integer j
real(8), dimension(:), allocatable :: illuminatedApertureAreas2
real(8), dimension(:), allocatable :: illuminatedApertureAreas2_ps

allocate(illuminatedApertureAreas2(1:size(apertureMidpoints,1)))
allocate(illuminatedApertureAreas2_ps(1:size(apertureMidpoints,1)))

! rotate into the aperture
call rotate_into_aperture_system(   aperturePropagationVectors, &
                                    apertureMidpoints, apertureNormals, &
                                    verts, Norm, midPoints, &
                                    rotatedapertureMidpoints, &
                                    rotatedapertureNormals, &
                                    rotatedaperturePropagationVectors, & ! unused
                                    rotatedVert, & 
                                    rotatedNorm, &
                                    rotatedMidpoints, &
                                    rotationMatrices, &
                                    beam)

! find the facets illuminated by the beam
call find_vis(  rotatedVert, Face1, Face2, &
                apertures, &
                rotatedMidpoints, rotatedNorm, &
                in_beam, dist_beam, id_beam, is_shad, &
                rotatedapertureNormals, &
                is_multithreaded, &
                beam)

! init
! totalIlluminated = 0
totalIlluminated_ps = 0
! illuminatedApertureAreas2 = 0
illuminatedApertureAreas2_ps = 0

do i = 1, size(apertureMidpoints,1) ! for each aperture
    ! counter = 0
    counter2 = 0
    do j = 1, size(face2,1) ! for each face
        ! if(isWithinBeam2(j) .and. apertures(j) .eq. i) then
        !     counter = counter + 1
        !     illuminatedApertureAreas2(i) = illuminatedApertureAreas2(i) + faceAreas(j)            
        ! end if
        if(in_beam(j) .and. apertures(j) .eq. i) then
            counter2 = counter2 + 1
            illuminatedApertureAreas2_ps(i) = illuminatedApertureAreas2_ps(i) + faceAreas(j)
        end if
    end do
    ! if(counter .gt. 0) then
    !     if(prescan) then
    !         print'(A,I5,A,I5,A40,I5,A,f10.6)','aperture ',aperture_id,' illuminated ',counter,' facets in aperture ',i,' with a total area of ',illuminatedApertureAreas2(i)
    !     end if
    !     totalIlluminated = totalIlluminated + counter
    ! end if
    if(counter2 .gt. 0) then
        ! print'(A,I5,A,I5,A40,I5,A,f10.6)','aperture ',aperture_id,' illuminated ',counter2,' facets (shadow) in aperture ',i,' with a total area of ',illuminatedApertureAreas2_ps(i)
        totalIlluminated_ps = totalIlluminated_ps + counter2
    end if    
end do
! print*,'checkpoint #2'

! if(totalIlluminated .gt. maxIlluminatedFacets) maxIlluminatedFacets = totalIlluminated
! if(totalIlluminated_ps .gt. maxIlluminatedFacets_ps) maxIlluminatedFacets_ps = totalIlluminated_ps

! init
sufficientlyIlluminated2 = .false.

do i = 1, size(apertureMidpoints,1) ! for each aperture
    if(illuminatedApertureAreas2_ps(i) .gt. threshold) then
        sufficientlyIlluminated2(i) = .true.
    end if
end do
! print*,'checkpoint #3'

! implementation for rough crystals - need to loop through apertures and disregard ones which face the wrong way
do i = 1, size(apertureMidpoints,1) ! for each aperture
    if(sufficientlyIlluminated2(i)) then
        if(rotatedapertureNormals(i,3) > -0.01) then ! if normal faces upwards, it is unphysical for internal interactions
            sufficientlyIlluminated2(i) = .false.
            ! print'(A,I2,A)','aperture ',i,' was suff. illuminated but set not visible because it was facing the wrong way'
        end if
    end if
end do
! print*,'checkpoint #4'

! do i = 1, size(apertureMidpoints,1) ! for each aperture
!     if(sufficientlyIlluminated2(i)) totalIlluminatedApertures = totalIlluminatedApertures + 1
! end do
! print*,'checkpoint #5'

! due to the above 2 do loops, the total ill. facets is probably an overestimate as some are ignored due to facing wrong way
! fix this later but for now we will just have a slightly too big array (i think)

! stop

end subroutine

subroutine beam_recursion(  sufficientlyIlluminated, &
                            aperturePropagationVectors, apertureMidpoints, apertureNormals, &
                            verts, Norm, midPoints, Face1, Face2, &
                            isWithinBeam, &
                            apertures, &
                            faceAreas, &
                            threshold, &
                            rbi, ibi, waveno, &
                            vk71, vk72, vk73, &
                            ampl, &
                            propagationVectors2, &
                            vk91Int, vk92Int, vk93Int, &
                            trans_ampl_out11_2, trans_ampl_out12_2, trans_ampl_out21_2, trans_ampl_out22_2, &
                            refl_ampl_out11_2, refl_ampl_out12_2, refl_ampl_out21_2, refl_ampl_out22_2, &
                            vk71Int2 ,vk72Int2, vk73Int2, &
                            vk121Int ,vk122Int, vk123Int, &
                            FInt, &
                            is_multithreaded)

! monster subroutine for the inner part of each beam recursion

logical, dimension(:), allocatable, intent(in) :: sufficientlyIlluminated
real(8), dimension(:,:), allocatable, intent(in) :: aperturePropagationVectors
real(8), dimension(:,:), allocatable, intent(in) :: apertureMidpoints ! the midpoint of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints ! face midpoints
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
logical, dimension(:), allocatable, intent(in) :: isWithinBeam
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), dimension(:), allocatable, intent(in) :: faceAreas ! area of each facet
real(8), intent(in) :: threshold ! minimum area of illumination per aperture to create new beam
real(8), intent(in) :: rbi ! real part of the refractive index
real(8), intent(in) :: ibi ! imaginary part of the refractive index
real(8), intent(in) :: waveno ! wavenumber in vacuum
real(8), dimension(:), allocatable, intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
complex(8), dimension(:,:,:), allocatable, intent(in) :: ampl
real(8), dimension(:,:,:), allocatable, intent(out) :: propagationVectors2
real(8), dimension(:,:), allocatable, intent(out) :: vk91Int, vk92Int, vk93Int
complex(8), dimension(:,:), allocatable, intent(out) :: trans_ampl_out11_2, trans_ampl_out12_2, trans_ampl_out21_2, trans_ampl_out22_2
complex(8), dimension(:,:), allocatable, intent(out) :: refl_ampl_out11_2, refl_ampl_out12_2, refl_ampl_out21_2, refl_ampl_out22_2
real(8), dimension(:,:), allocatable, intent(out) :: vk71Int2, vk72Int2, vk73Int2
real(8), dimension(:,:), allocatable, intent(out) :: vk121Int ,vk122Int, vk123Int
integer, dimension(:,:), allocatable, intent(out) :: FInt
logical, intent(in) :: is_multithreaded

! integer maxIlluminatedFacets ! counts the maximum number of facets illuminated across all illuminating apertures
integer maxIlluminatedFacets_ps ! counts the maximum number of facets illuminated across all illuminating apertures (including shadow)
integer totalIlluminatedApertures ! total number of new apertures illuminated
integer i, counter, num_sufficiently_illuminated_apertures, j, k, l
integer aperture_id
integer, dimension(:), allocatable :: sufficiently_illuminated_indices ! indices of sufficiently illuminated apertures
logical, dimension(:), allocatable :: sufficientlyIlluminated2 ! whether each aperture was sufficiently illuminated by the previous beam
integer sum_suff_ill ! total number of new sufficiently illuminated apertures
logical prescan
logical, dimension(:), allocatable :: isShadow
logical, dimension(:), allocatable :: in_beam
real(8), dimension(:,:), allocatable :: rotatedapertureNormals
real(8), dimension(:,:,:), allocatable :: rotationMatrices
real(8) vk7a(1:3)
integer aperture
real(8) rot(1:2,1:2)
complex(8) rot_ampl(1:2,1:2)
integer, dimension(:), allocatable :: id_beam
real(8), dimension(:), allocatable :: dist_beam
real(8) theta_i, theta_t, theta_i_aperture
logical tir
complex(8) r_perp, t_perp, r_par, t_par ! fresnel coefficients
complex(8) fresnel_matrix_refl(1:2,1:2), fresnel_matrix_trans(1:2,1:2)
complex(8) refl_ampl(1:2,1:2), trans_ampl(1:2,1:2)
real(8) alpha, A, B
real(8) a_vec(1:3), b_vec(1:3), c_vec(1:3)
real(8) intensity_in, intensity_abs, intensity_out
! real(8) start, finish

type(beam_type) beam ! current beam to be traced



! allocations
allocate(sufficientlyIlluminated2(1:size(sufficientlyIlluminated,1)))
allocate(rotationMatrices(1:3,1:3,1:size(apertureMidpoints,1)))
allocate(beam%field(1:size(Face1,1))) ! reduce this to just the size of each beam to preserve memory

! make a quick array containing the sufficiently illuminated aperture numbers
num_sufficiently_illuminated_apertures = 0
do i = 1, size(sufficientlyIlluminated,1)
    if(sufficientlyIlluminated(i)) num_sufficiently_illuminated_apertures = num_sufficiently_illuminated_apertures + 1
end do
allocate(sufficiently_illuminated_indices(1:num_sufficiently_illuminated_apertures))
counter = 0
do i = 1, size(sufficientlyIlluminated,1)
    if(sufficientlyIlluminated(i)) then
        counter = counter + 1
        sufficiently_illuminated_indices(counter) = i
    end if
end do
! print*,'sufficientlyIlluminated',sufficientlyIlluminated
! print*,'sufficiently_illuminated_indices',sufficiently_illuminated_indices(i)

! init
! maxIlluminatedFacets = 0
maxIlluminatedFacets_ps = 0
totalIlluminatedApertures = 0
rotationMatrices = 0

! start = omp_get_wtime()

! perform prescan to determine shape of arrays (removed)
! prescan = .true.
! do i = 1, num_sufficiently_illuminated_apertures
!     aperture_id = sufficiently_illuminated_indices(i)
!     call beam_scan( aperturePropagationVectors, apertureMidpoints, apertureNormals, &
!                     verts, Norm, midPoints, Face1, Face2, &
!                     aperture_id, &
!                     isWithinBeam, &
!                     apertures, &
!                     faceAreas, &
!                     maxIlluminatedFacets_ps, &
!                     threshold, &
!                     sufficientlyIlluminated2, &
!                     totalIlluminatedApertures, &
!                     prescan, &
!                     beamFaceAreas, & ! dummy
!                     beamvk71, beamvk72, beamvk73, & ! dummy
!                     beam_ampl, & ! dummy
!                     isWithinBeam2, isWithinBeam2_ps, & 
!                     rotatedapertureNormals, &
!                     rotationMatrices, &
!                     beamIDs_ps, distances_ps, &
!                     vk71, vk72, vk73, ampl, &
!                     rotatedVert, rotatedNorm, &
!                     isShadow, &
!                     is_multithreaded) 
! end do

! print*,'total new illuminated apertures: ',totalIlluminatedApertures
! print*,'max number of illuminated facets: ',maxIlluminatedFacets
! print*,'max number of illuminated facets (shadow): ',maxIlluminatedFacets_ps

! prescan complete, now we know the sizes to allocate for arrays

maxIlluminatedFacets_ps = size(Face1,1) ! override prescan

! get total number of illuminating apertures
sum_suff_ill = 0 
do i = 1,size(sufficientlyIlluminated,1)
    if(sufficientlyIlluminated(i)) then
        sum_suff_ill = sum_suff_ill + 1
    end if
end do

! print*,'sufficientlyIlluminated',sufficientlyIlluminated
! stop


! allocations
allocate(vk91Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk92Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk93Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk71Int2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk72Int2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk73Int2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk121Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk122Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(vk123Int(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(trans_ampl_out11_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(trans_ampl_out12_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(trans_ampl_out21_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(trans_ampl_out22_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(refl_ampl_out11_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(refl_ampl_out12_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(refl_ampl_out21_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(refl_ampl_out22_2(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(FInt(1:maxIlluminatedFacets_ps,1:sum_suff_ill))
allocate(propagationVectors2(1:size(apertureMidpoints,1),1:3,1:sum_suff_ill))

! init
vk91Int = 0
vk92Int = 0
vk93Int = 0
vk71Int2 = 0
vk72Int2 = 0
vk73Int2 = 0
vk121Int = 0
vk122Int = 0
vk123Int = 0
trans_ampl_out11_2 = 0
trans_ampl_out12_2 = 0
trans_ampl_out21_2 = 0
trans_ampl_out22_2 = 0
refl_ampl_out11_2 = 0
refl_ampl_out12_2 = 0
refl_ampl_out21_2 = 0
refl_ampl_out22_2 = 0
FInt = 0
propagationVectors2 = 0

l = 0 ! counting variable to count through apertures
prescan = .false.
do i = 1, num_sufficiently_illuminated_apertures
    l = l + 1 ! update counter
    aperture_id = sufficiently_illuminated_indices(i) ! replaces i in matlab code



    ! populate the beam data structure
    counter = 0 ! initialise a counter
    do j = 1, size(Face2,1) ! for each facet
        if(isWithinBeam(j) .and. apertures(j) .eq. aperture_id) then ! if the facet was within the beam of the previous recursion, and belongs to this aperture
            counter = counter + 1
            beam%field(counter)%face_id = j
            beam%field(counter)%e_perp = (/vk71(j),vk72(j),vk73(j)/)
            beam%field(counter)%ampl = ampl(1:2,1:2,j)
        end if
    end do

    beam%prop(1:3) = aperturePropagationVectors(aperture_id,1:3)
    beam%num_facets = counter
    beam%aperture_id = aperture_id

    call beam_scan( aperturePropagationVectors, &
                    apertureMidpoints, apertureNormals, &
                    verts, Norm, midPoints, Face1, Face2, &
                    apertures, &
                    faceAreas, &
                    threshold, &
                    sufficientlyIlluminated2, &
                    in_beam, &
                    rotatedapertureNormals, &
                    rotationMatrices, &
                    id_beam, dist_beam, &
                    isShadow, &
                    is_multithreaded, &
                    beam)

    ! get the reflected propagation vectors for the current illuminating beam
    call get_reflection_vectors(propagationVectors2, rotatedapertureNormals, sufficientlyIlluminated2, l)

    ! now rotate the propagation vectors back to the original coordinate system for later use
    call rotate_back_propagation_vectors(rotationMatrices, propagationVectors2, sufficientlyIlluminated2, aperture_id, l)

    k = 0 ! counting variable for storing e-field components later on

    do j = 1, size(Face2,1) ! for each facet
        if(in_beam(j)) then
            k = k + 1
            ! if(isShadow(j) .eqv. .false.) then
            !     call cross(-Norm(Face2(j),1:3),aperturePropagationVectors(aperture_id,1:3),vk7a(1:3))
            ! else
                aperture = apertures(j)
                call cross(-apertureNormals(aperture,1:3),aperturePropagationVectors(aperture_id,1:3),vk7a(1:3))
            ! end if

            call getRotationMatrix( rot,vk7a(1),vk7a(2),vk7a(3), &
                                    ! beamvk71(id_beam(j)),beamvk72(id_beam(j)),beamvk73(id_beam(j)), &
                                    beam%field(id_beam(j))%e_perp(1),beam%field(id_beam(j))%e_perp(2),beam%field(id_beam(j))%e_perp(3), &
                                    aperturePropagationVectors(aperture_id,1),aperturePropagationVectors(aperture_id,2),aperturePropagationVectors(aperture_id,3))

            rot_ampl(1:2,1:2) = matmul(rot(1:2,1:2),beam%field(id_beam(j))%ampl(:,:))

            intensity_in = real(0.5*(   rot_ampl(1,1)*conjg(rot_ampl(1,1)) + &
                                        rot_ampl(1,2)*conjg(rot_ampl(1,2)) + &
                                        rot_ampl(2,1)*conjg(rot_ampl(2,1)) + &
                                        rot_ampl(2,2)*conjg(rot_ampl(2,2))))

            rot_ampl = rot_ampl * exp2cmplx(waveno*rbi*dist_beam(j)) * exp(-2*waveno*ibi*sqrt(dist_beam(j))) ! absorption

            intensity_out = real(0.5*(  rot_ampl(1,1)*conjg(rot_ampl(1,1)) + &
                                        rot_ampl(1,2)*conjg(rot_ampl(1,2)) + &
                                        rot_ampl(2,1)*conjg(rot_ampl(2,1)) + &
                                        rot_ampl(2,2)*conjg(rot_ampl(2,2))))

            intensity_abs = real(0.5*(  beam%field(id_beam(j))%ampl(1,1)*conjg(beam%field(id_beam(j))%ampl(1,1)) + &
                                        beam%field(id_beam(j))%ampl(1,2)*conjg(beam%field(id_beam(j))%ampl(1,2)) + &
                                        beam%field(id_beam(j))%ampl(2,1)*conjg(beam%field(id_beam(j))%ampl(2,1)) + &
                                        beam%field(id_beam(j))%ampl(2,2)*conjg(beam%field(id_beam(j))%ampl(2,2))))*(1-exp(-2*waveno*ibi*sqrt(dist_beam(j)))**2)

            ! print*,'intensity conservation: ',(intensity_out + intensity_abs)/intensity_in * 100,'%'

            ! print'(A16,f10.4,A,f10.4,A,f10.4,A,f10.4,A)','     rot ampl: (',real(rot_ampl(1,1)),' + ',imag(rot_ampl(1,1)),'i, ',real(rot_ampl(1,2)),' + ',imag(rot_ampl(1,2)),'i)'
            ! print'(A16,f10.4,A,f10.4,A,f10.4,A,f10.4,A)','               (',real(rot_ampl(2,1)),' + ',imag(rot_ampl(2,1)),'i, ',real(rot_ampl(2,2)),' + ',imag(rot_ampl(2,2)),'i)'

            ! FRESNEL
            ! if(isShadow(j) .eqv. .false.) then
            !     theta_i = acos(-rotatedNorm(Face2(j),3)) ! get incident angle
            ! else
                aperture = apertures(j) ! probably dont need this line (see above)
                theta_i = acos(-rotatedapertureNormals(aperture,3)) ! get incident angle                
            ! end if
            if(theta_i .gt. asin(1/rbi)) then ! if tir
                tir = .true.
                r_perp = -1
                t_perp = 0
                r_par = -1
                t_par = 0
            else
                tir = .false.
                theta_t = asin(sin(theta_i)*rbi)
                r_perp = (rbi*cos(theta_i) - cos(theta_t))/(rbi*cos(theta_i) + cos(theta_t))
                t_perp = (2*rbi*cos(theta_i))/(rbi*cos(theta_i) + cos(theta_t))
                r_par = (cos(theta_i) - rbi*cos(theta_t))/(rbi*cos(theta_t) + cos(theta_i))
                t_par = (2*rbi*cos(theta_i)) / (rbi*cos(theta_t) + cos(theta_i))
            end if
            ! force TIR if angle with aperture normal is > critical angle
            aperture = apertures(j)
            theta_i_aperture = acos(-rotatedapertureNormals(aperture,3))
            if(theta_i_aperture .gt. asin(1/rbi)) then
                r_perp = -1
                t_perp = 0
                r_par = -1
                t_par = 0
            end if
            ! make fresnel matrix
            fresnel_matrix_refl = 0
            fresnel_matrix_trans = 0
            fresnel_matrix_refl(1,1) = r_par
            fresnel_matrix_refl(2,2) = r_perp
            fresnel_matrix_trans(1,1) = t_par
            fresnel_matrix_trans(2,2) = t_perp
            ! rotate fresnel matrix into new scattering plane
            trans_ampl = matmul(fresnel_matrix_trans,rot_ampl)
            refl_ampl = matmul(fresnel_matrix_refl,rot_ampl)
            ! remove transmission for shadowed facets (probably doesnt matter)
            if(isShadow(j)) trans_ampl = trans_ampl*0

            ! print*,'illuminated facet ID: ',j
            ! print*,'illuminated facet cos(theta_i): ',cos(theta_i)
            ! print*,'illuminated face area: ',faceAreas(j)
            ! print*,'illuminated projected face area: ',cos(theta_i)*faceAreas(j)
            ! print*,'illuminating facet ID: ',id_beam(j)
            ! print*,'illuminating facet area: ',faceAreas(id_beam(j))
            ! print*,'illuminating facet cos(theta_i)',-rotatedNorm(Face2(id_beam(j)),3)
            ! print*,'illuminating projected face area: ',-rotatedNorm(Face2(id_beam(j)),3)*faceAreas(id_beam(j))
            ! print*,'absorption factor: ',(1-exp(-2*waveno*ibi*dist_beam(j)))
            ! print*,'illuminated extinction cross section: ',(1-exp(-2*waveno*ibi*dist_beam(j)))*cos(theta_i)*faceAreas(j)

            ! ext_cross_section = ext_cross_section + sqrt(1-exp(-2*waveno*ibi*dist_beam(j)))*cos(theta_i)*faceAreas(j)*intensity_in*rbi ! this might need a check
            
            ext_cross_section = ext_cross_section + intensity_abs*cos(theta_i)*faceAreas(j)*rbi ! this might need a check, rbi because energy is more concentrated inside particle

            ! ext_cross_section = ext_cross_section + intensity_abs*cos(theta_i)*faceAreas(j)*rbi/(-rotatedNorm(Face2(id_beam(j)),3)) ! this might need a check
                                        ! stop

            ! output variables
            vk91Int(k,l) = aperturePropagationVectors(aperture_id,1)
            vk92Int(k,l) = aperturePropagationVectors(aperture_id,2)
            vk93Int(k,l) = aperturePropagationVectors(aperture_id,3)
            vk71Int2(k,l) = vk7a(1)
            vk72Int2(k,l) = vk7a(2)
            vk73Int2(k,l) = vk7a(3)
            FInt(k,l) = j
            trans_ampl_out11_2(k,l) = trans_ampl(1,1)
            trans_ampl_out12_2(k,l) = trans_ampl(1,2)
            trans_ampl_out21_2(k,l) = trans_ampl(2,1)
            trans_ampl_out22_2(k,l) = trans_ampl(2,2)
            refl_ampl_out11_2(k,l) = refl_ampl(1,1)
            refl_ampl_out12_2(k,l) = refl_ampl(1,2)
            refl_ampl_out21_2(k,l) = refl_ampl(2,1)
            refl_ampl_out22_2(k,l) = refl_ampl(2,2)

            ! calc transmission propagation vector for vector diffraction
            if(tir) then
                vk121Int(k,l) = 0
                vk122Int(k,l) = 0
                vk123Int(k,l) = 1
            else
                alpha = pi - theta_t
                A = sin(theta_t-theta_i)/sin(theta_i)
                B = sin(alpha)/sin(theta_i)
                b_vec(1) = aperturePropagationVectors(aperture_id,1)
                b_vec(2) = aperturePropagationVectors(aperture_id,2)
                b_vec(3) = aperturePropagationVectors(aperture_id,3)
                a_vec(1) = Norm(Face2(j),1)
                a_vec(2) = Norm(Face2(j),2)
                a_vec(3) = Norm(Face2(j),3)
                c_vec(1:3) = B*b_vec(1:3) - A*a_vec(1:3)
                c_vec(1:3) = c_vec(1:3) / sqrt(c_vec(1)**2 + c_vec(2)**2 + c_vec(3)**2)
                vk121Int(k,l) = c_vec(1)
                vk122Int(k,l) = c_vec(2)
                vk123Int(k,l) = c_vec(3)
                ! debugging
                ! if(k .eq. 268 .and. l .eq. 3) then
                !     print*,'alpha',alpha
                !     print*,'theta_i',theta_i
                !     print*,'theta_t',theta_t
                !     print*,'A',A
                !     print*,'B',B
                !     print*,'a_vec',a_vec
                !     print*,'b_vec',b_vec
                !     print*,'c_vec',c_vec
                !     print*,'vk121Int(k,l)',vk121Int(k,l)
                !     print*,'aperturePropagationVectors(aperture_id,1:3)',aperturePropagationVectors(aperture_id,1:3)
                !     stop
                ! end if

            end if


            ! stop


        end if
    end do


end do

! finish = omp_get_wtime()

! print'(A,f16.8,A)',"total beam time elapsed: ",finish-start," secs"

! stop

end subroutine

subroutine internal_recursion_outer(F1Mapping, &
                                    propagationVectors, &
                                    verts, Norm, midPoints, Face1, Face2, &
                                    apertureMidpoints, apertureNormals, apertures, &
                                    faceAreas, &
                                    threshold, &
                                    rbi, ibi, waveno, &
                                    vk71Int, vk72Int, vk73Int, &
                                    beam_outbeam_tree, beam_outbeam_tree_counter, &
                                    refl_ampl_out11Int, refl_ampl_out12Int, refl_ampl_out21Int, refl_ampl_out22Int, &
                                    interactionCounter, &
                                    is_multithreaded)

! monster subroutine for the outer part of each beam recursion
! vk71, vk72, and vk73 appear to have errors and have been incorrectly implemented (also in Matlab code)

integer, dimension(:,:), allocatable, intent(inout) :: F1Mapping ! the face indices of each row of variables to go into main loop
real(8), dimension(:,:,:), allocatable, intent(inout) :: propagationVectors ! propagation vector of beam emitted from each aperture
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints ! face midpoints
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
real(8), dimension(:,:), allocatable, intent(in) :: apertureMidpoints ! the midpoint of each aperture
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), dimension(:), allocatable, intent(in) :: faceAreas ! area of each facet
real(8), intent(in) :: threshold ! minimum area of illumination per aperture to create new beam
real(8), intent(in) :: rbi ! real part of the refractive index
real(8), intent(in) :: ibi ! imaginary part of the refractive index
real(8), intent(in) :: waveno ! wavenumber in vacuum
real(8), dimension(:,:), allocatable, intent(inout) :: vk71Int, vk72Int, vk73Int ! reflected e-perp vector
type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
integer, intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
complex(8), dimension(:,:), allocatable, intent(inout) :: refl_ampl_out11Int ! the amplitude matrix that goes into the main loop
complex(8), dimension(:,:), allocatable, intent(inout) :: refl_ampl_out12Int ! needs renaming
complex(8), dimension(:,:), allocatable, intent(inout) :: refl_ampl_out21Int ! needs renaming
complex(8), dimension(:,:), allocatable, intent(inout) :: refl_ampl_out22Int ! needs renaming
integer, intent(inout) :: interactionCounter ! counts the current number of interactions
logical, intent(in) :: is_multithreaded

real(8), dimension(:,:,:), allocatable :: propagationVectors2
real(8), dimension(:,:,:), allocatable :: propagationVectors3
real(8), dimension(:,:), allocatable :: vk91Int2, vk92Int2, vk93Int2
! real(8), dimension(:), allocatable :: vk91Int1, vk92Int1, vk93Int1
complex(8), dimension(:,:), allocatable :: trans_ampl_out11_2, trans_ampl_out12_2, trans_ampl_out21_2, trans_ampl_out22_2
! complex(8), dimension(:), allocatable :: trans_ampl_out11_1, trans_ampl_out12_1, trans_ampl_out21_1, trans_ampl_out22_1
complex(8), dimension(:,:), allocatable :: refl_ampl_out11_2, refl_ampl_out12_2, refl_ampl_out21_2, refl_ampl_out22_2
complex(8), dimension(:,:), allocatable :: refl_ampl_out11_3, refl_ampl_out12_3, refl_ampl_out21_3, refl_ampl_out22_3
real(8), dimension(:,:), allocatable :: vk71Int2, vk72Int2, vk73Int2
real(8), dimension(:,:), allocatable :: vk71Int3, vk72Int3, vk73Int3
! real(8), dimension(:), allocatable :: vk71Int1, vk72Int1, vk73Int1
real(8), dimension(:,:), allocatable :: vk121Int2 ,vk122Int2, vk123Int2
! real(8), dimension(:), allocatable :: vk121Int1 ,vk122Int1, vk123Int1
integer, dimension(:,:), allocatable :: FInt2
integer, dimension(:,:), allocatable :: FInt3
! integer, dimension(:), allocatable :: FInt1
integer, dimension(:,:), allocatable :: InteractionInt2
! integer, dimension(:), allocatable :: InteractionInt1

integer i, j, k, m
integer, dimension(:), allocatable :: illuminatedFaceIDs
logical, dimension(:), allocatable :: isWithinBeam
real(8), dimension(:,:), allocatable :: aperturePropagationVectors
logical, dimension(:), allocatable :: sufficientlyIlluminated
complex(8), dimension(:,:,:), allocatable :: ampl
real(8), dimension(:), allocatable :: vk71, vk72, vk73
! integer FTemp, counter
! integer, dimension(:), allocatable :: illuminated_apertures_temp
! integer, dimension(:), allocatable :: unique_illuminated_apertures
! logical, dimension(:), allocatable :: FInt1_mask
! real(8), dimension(:), allocatable :: remaining_aperture_energy
! integer(8) my_aperture_id
real(8) start, finish
logical is_first_beam_back ! whether or not this is the first loop iteration back
integer num_req_threads ! number of required threads


! start = omp_get_wtime()

! allocations
allocate(illuminatedFaceIDs(1:size(F1Mapping,1))) ! allocate array to hold each column of the Face1 Mappings
allocate(isWithinBeam(1:size(Face2,1))) ! which faces were within beam
allocate(ampl(1:2,1:2,1:size(Face2,1))) ! amplitude matrix
allocate(vk71(1:size(Face2,1))) ! e-perp direction
allocate(vk72(1:size(Face2,1))) ! e-perp direction
allocate(vk73(1:size(Face2,1))) ! e-perp direction
allocate(aperturePropagationVectors(1:size(propagationVectors,1),1:3)) ! propagation vectors for each internal field
allocate(sufficientlyIlluminated(1:size(propagationVectors,1)))

if (is_multithreaded) then ! multi-threaded diffraction

    ! get desired number of threads
    num_req_threads = min(omp_get_max_threads(),size(F1Mapping,2))
    ! num_req_threads = 1
    !$OMP PARALLEL num_threads(num_req_threads) ,PRIVATE( &
    !$OMP   i,is_first_beam_back,isWithinBeam,ampl,vk71,vk72,vk73, &
    !$OMP   aperturePropagationVectors,sufficientlyIlluminated,illuminatedFaceIDs, &
    !$OMP   propagationVectors2,vk91Int2, vk92Int2, vk93Int2, &
    !$OMP   trans_ampl_out11_2,trans_ampl_out12_2,trans_ampl_out21_2,trans_ampl_out22_2, &
    !$OMP   refl_ampl_out11_2,refl_ampl_out12_2,refl_ampl_out21_2,refl_ampl_out22_2, &
    !$OMP   vk71Int2,vk72Int2,vk73Int2,vk121Int2 ,vk122Int2, vk123Int2,FInt2 &
    !$OMP )
    !$OMP DO
    do i = 1, size(F1Mapping,2) ! looping over internal fields created by the previous recursion

        is_first_beam_back = .false. ! init

        ! step 1: retrieve the parameters needed to propagate the next set of beams

        call get_beam_params(   F1Mapping,                      & ! <-  the face IDs of all illuminated faces from all beams of the previous recursion
                                isWithinBeam,                   & ! <-> logical array of the faces illuminated by beam i from the previous recursion
                                ampl,                           & ! <-> amplitude matrix of the faces illuminated by beam i from the previous recursion
                                vk71,                           & ! <-> e-perp x component of the faces illuminated by beam i from the previous recursion
                                vk72,                           & ! <-> e-perp y component of the faces illuminated by beam i from the previous recursion
                                vk73,                           & ! <-> e-perp z component of the faces illuminated by beam i from the previous recursion
                                aperturePropagationVectors,     & ! <-> propagation vectors of all apertures illuminated by beam i from the previous recursion
                                sufficientlyIlluminated,        & ! <-> logical array of the apertures which were sufficiently illuminated by beam i from the previous recursion
                                propagationVectors,             & ! <-  the propagation vectors of all apertures for all beams from the previous recursion
                                i,                              & ! <-  the beam number from the previous recursion that we wish to propagate the reflected beams of
                                illuminatedFaceIDs,             & ! <-> the face IDs of all illuminated faces from beam i of the previous recursion
                                refl_ampl_out11Int,             & ! <-  the amplitude matrix (1,1) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out12Int,             & ! <-  the amplitude matrix (1,2) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out21Int,             & ! <-  the amplitude matrix (2,1) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out22Int,             & ! <-  the amplitude matrix (2,2) of all illuminated faces from all beams of the previous recursion
                                vk71Int, vk72Int, vk73Int)        ! <-  the e-perp components of all illuminated faces from all beams of the previous recursion

        ! step 2: propagate the next set of beams (all apertures illuminated by a beam from the previous recursion)

        call beam_recursion(    sufficientlyIlluminated,                & ! <-  logical array of which apertures were sufficiently illuminated by beam i from the previous recursion
                                aperturePropagationVectors,             & ! <-  propagation vectors of all apertures illuminated by beam i from the previous recursion
                                apertureMidpoints,                      & ! <-  midpoints of apertures
                                apertureNormals,                        & ! <-  average normals of apertures
                                verts, Norm, midPoints, Face1, Face2,   & ! <-  particle vertices, normals, midpoints, face ids, normal ids
                                isWithinBeam,                           & ! <-  logical array of the faces illuminated by beam i from the previous recursion
                                apertures,                              & ! <-  which aperture each face belongs to
                                faceAreas,                              & ! <-  area of each face
                                threshold,                              & ! <-  threshold illuminated area, under which new beams will not propagate
                                rbi, ibi, waveno,                       & ! <-  input parameters
                                vk71, vk72, vk73,                       & ! <-  e-perp components of the faces illuminated by beam i from the previous recursion
                                ampl,                                   & ! <-  amplitude matrix of the faces illuminated by beam i from the previous recursion
                                propagationVectors2,                    & !  -> the propagation vectors of all apertures illuminated from this recursion
                                vk91Int2, vk92Int2, vk93Int2,           & !  -> the reflected propagation direction components at all faces illuminated from this recursion
                                trans_ampl_out11_2,                     & !  -> the transmitted amplitude matrix (1,1) at all faces illuminated from this recursion
                                trans_ampl_out12_2,                     & !  -> the transmitted amplitude matrix (1,2) at all faces illuminated from this recursion
                                trans_ampl_out21_2,                     & !  -> the transmitted amplitude matrix (2,1) at all faces illuminated from this recursion
                                trans_ampl_out22_2,                     & !  -> the transmitted amplitude matrix (2,2) at all faces illuminated from this recursion
                                refl_ampl_out11_2,                      & !  -> the reflected amplitude matrix (1,1) at all faces illuminated from this recursion
                                refl_ampl_out12_2,                      & !  -> the reflected amplitude matrix (1,2) at all faces illuminated from this recursion
                                refl_ampl_out21_2,                      & !  -> the reflected amplitude matrix (2,1) at all faces illuminated from this recursion
                                refl_ampl_out22_2,                      & !  -> the reflected amplitude matrix (2,2) at all faces illuminated from this recursion
                                vk71Int2 ,vk72Int2, vk73Int2,           & !  -> the e-perp components of all faces illuminated from this recursion
                                vk121Int2 ,vk122Int2, vk123Int2,        & !  -> the transmitted propagation direction components at all faces illuminated from this recursion
                                FInt2,                                  & !  -> the face IDs of all faces illuminated from this recursion
                                is_multithreaded)                         ! <-  whether or not the beam recursion should use multithreaded operations (currently disabled)

        ! step 3: keep track of the different interaction numbers (optional)
        !$OMP CRITICAL
        call get_interaction(   InteractionInt2,    & ! the interaction number of each illuminated face at this recursion
                                FInt2,              & ! the face IDs of all faces illuminated from this recursion
                                apertures,          & ! which aperture each face belongs to
                                interactionCounter)   ! total number of interactions so far

        ! step 4: add outgoing rays to the outbeam tree

        call add_to_outbeam_tree(   FInt2,                              & ! <-  the face IDs of all faces illuminated from this recursion
                                    beam_outbeam_tree_counter,          & ! <-> total number of entries in outbeam tree
                                    beam_outbeam_tree,                  & ! <-> outbeam tree
                                    trans_ampl_out11_2,                 & ! <-  the transmitted amplitude matrix (1,1) at all faces illuminated from this recursion
                                    trans_ampl_out12_2,                 & ! <-  the transmitted amplitude matrix (1,2) at all faces illuminated from this recursion
                                    trans_ampl_out21_2,                 & ! <-  the transmitted amplitude matrix (2,1) at all faces illuminated from this recursion
                                    trans_ampl_out22_2,                 & ! <-  the transmitted amplitude matrix (2,2) at all faces illuminated from this recursion
                                    vk71Int2, vk72Int2, vk73Int2,       & ! <-  the e-perp components of all faces illuminated from this recursion
                                    vk121Int2, vk122Int2, vk123Int2,    & ! <-  the transmitted propagation direction components at all faces illuminated from this recursion
                                    vk91Int2, vk92Int2, vk93Int2,       & ! <-  the reflected propagation direction components at all faces illuminated from this recursion
                                    InteractionInt2)                      ! <-  the interaction number of each illuminated face at this recursion

        ! step 5: stitch together some arrays for use in the next recursion of the beam loop

        if(.not. allocated(vk71Int3)) is_first_beam_back = .true. ! determine whether this is the first beam back (multithreading support)

        ! ! init v3 variables
        if(is_first_beam_back) then ! if its the first beam back, we need to allocate first
            allocate(vk71Int3(1:size(FInt2,1),size(FInt2,2)))
            allocate(vk72Int3(1:size(FInt2,1),size(FInt2,2)))
            allocate(vk73Int3(1:size(FInt2,1),size(FInt2,2)))  
            allocate(refl_ampl_out11_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out12_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out21_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out22_3(1:size(FInt2,1),size(FInt2,2)))  
            allocate(FInt3(1:size(FInt2,1),size(FInt2,2)))
            FInt3 = FInt2 ! keep track, ready to be concatenated to on the next loop
            vk71Int3 = vk71Int2
            vk72Int3 = vk72Int2
            vk73Int3 = vk73Int2
            refl_ampl_out11_3 = refl_ampl_out11_2
            refl_ampl_out12_3 = refl_ampl_out12_2
            refl_ampl_out21_3 = refl_ampl_out21_2
            refl_ampl_out22_3 = refl_ampl_out22_2
            propagationVectors3 = propagationVectors2       
        else
            call cat_int_var(FInt2, FInt3)
            call cat_real_var(vk71Int2, vk71Int3)
            call cat_real_var(vk72Int2, vk72Int3)
            call cat_real_var(vk73Int2, vk73Int3)
            call cat_complex_var(refl_ampl_out11_2, refl_ampl_out11_3)
            call cat_complex_var(refl_ampl_out12_2, refl_ampl_out12_3)
            call cat_complex_var(refl_ampl_out21_2, refl_ampl_out21_3)
            call cat_complex_var(refl_ampl_out22_2, refl_ampl_out22_3)
            call cat_prop(propagationVectors2, propagationVectors3) ! concatenate propagation vectors
            ! stop
        end if
        !$OMP END CRITICAL
    end do
    !$OMP END PARALLEL
else ! single-threaded
    do i = 1, size(F1Mapping,2) ! looping over internal fields created by the previous recursion

        is_first_beam_back = .false. ! init

        ! step 1: retrieve the parameters needed to propagate the next set of beams

        call get_beam_params(   F1Mapping,                      & ! <-  the face IDs of all illuminated faces from all beams of the previous recursion
                                isWithinBeam,                   & ! <-> logical array of the faces illuminated by beam i from the previous recursion
                                ampl,                           & ! <-> amplitude matrix of the faces illuminated by beam i from the previous recursion
                                vk71,                           & ! <-> e-perp x component of the faces illuminated by beam i from the previous recursion
                                vk72,                           & ! <-> e-perp y component of the faces illuminated by beam i from the previous recursion
                                vk73,                           & ! <-> e-perp z component of the faces illuminated by beam i from the previous recursion
                                aperturePropagationVectors,     & ! <-> propagation vectors of all apertures illuminated by beam i from the previous recursion
                                sufficientlyIlluminated,        & ! <-> logical array of the apertures which were sufficiently illuminated by beam i from the previous recursion
                                propagationVectors,             & ! <-  the propagation vectors of all apertures for all beams from the previous recursion
                                i,                              & ! <-  the beam number from the previous recursion that we wish to propagate the reflected beams of
                                illuminatedFaceIDs,             & ! <-> the face IDs of all illuminated faces from beam i of the previous recursion
                                refl_ampl_out11Int,             & ! <-  the amplitude matrix (1,1) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out12Int,             & ! <-  the amplitude matrix (1,2) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out21Int,             & ! <-  the amplitude matrix (2,1) of all illuminated faces from all beams of the previous recursion
                                refl_ampl_out22Int,             & ! <-  the amplitude matrix (2,2) of all illuminated faces from all beams of the previous recursion
                                vk71Int, vk72Int, vk73Int)        ! <-  the e-perp components of all illuminated faces from all beams of the previous recursion

        ! step 2: propagate the next set of beams (all apertures illuminated by a beam from the previous recursion)

        call beam_recursion(    sufficientlyIlluminated,                & ! <-  logical array of which apertures were sufficiently illuminated by beam i from the previous recursion
                                aperturePropagationVectors,             & ! <-  propagation vectors of all apertures illuminated by beam i from the previous recursion
                                apertureMidpoints,                      & ! <-  midpoints of apertures
                                apertureNormals,                        & ! <-  average normals of apertures
                                verts, Norm, midPoints, Face1, Face2,   & ! <-  particle vertices, normals, midpoints, face ids, normal ids
                                isWithinBeam,                           & ! <-  logical array of the faces illuminated by beam i from the previous recursion
                                apertures,                              & ! <-  which aperture each face belongs to
                                faceAreas,                              & ! <-  area of each face
                                threshold,                              & ! <-  threshold illuminated area, under which new beams will not propagate
                                rbi, ibi, waveno,                       & ! <-  input parameters
                                vk71, vk72, vk73,                       & ! <-  e-perp components of the faces illuminated by beam i from the previous recursion
                                ampl,                                   & ! <-  amplitude matrix of the faces illuminated by beam i from the previous recursion
                                propagationVectors2,                    & !  -> the propagation vectors of all apertures illuminated from this recursion
                                vk91Int2, vk92Int2, vk93Int2,           & !  -> the reflected propagation direction components at all faces illuminated from this recursion
                                trans_ampl_out11_2,                     & !  -> the transmitted amplitude matrix (1,1) at all faces illuminated from this recursion
                                trans_ampl_out12_2,                     & !  -> the transmitted amplitude matrix (1,2) at all faces illuminated from this recursion
                                trans_ampl_out21_2,                     & !  -> the transmitted amplitude matrix (2,1) at all faces illuminated from this recursion
                                trans_ampl_out22_2,                     & !  -> the transmitted amplitude matrix (2,2) at all faces illuminated from this recursion
                                refl_ampl_out11_2,                      & !  -> the reflected amplitude matrix (1,1) at all faces illuminated from this recursion
                                refl_ampl_out12_2,                      & !  -> the reflected amplitude matrix (1,2) at all faces illuminated from this recursion
                                refl_ampl_out21_2,                      & !  -> the reflected amplitude matrix (2,1) at all faces illuminated from this recursion
                                refl_ampl_out22_2,                      & !  -> the reflected amplitude matrix (2,2) at all faces illuminated from this recursion
                                vk71Int2 ,vk72Int2, vk73Int2,           & !  -> the e-perp components of all faces illuminated from this recursion
                                vk121Int2 ,vk122Int2, vk123Int2,        & !  -> the transmitted propagation direction components at all faces illuminated from this recursion
                                FInt2,                                  & !  -> the face IDs of all faces illuminated from this recursion
                                is_multithreaded)                         ! <-  whether or not the beam recursion should use multithreaded operations (currently disabled)

        ! step 3: keep track of the different interaction numbers (optional)

        call get_interaction(   InteractionInt2,    & ! the interaction number of each illuminated face at this recursion
                                FInt2,              & ! the face IDs of all faces illuminated from this recursion
                                apertures,          & ! which aperture each face belongs to
                                interactionCounter)   ! total number of interactions so far

        ! step 4: add outgoing rays to the outbeam tree

        call add_to_outbeam_tree(   FInt2,                              & ! <-  the face IDs of all faces illuminated from this recursion
                                    beam_outbeam_tree_counter,          & ! <-> total number of entries in outbeam tree
                                    beam_outbeam_tree,                  & ! <-> outbeam tree
                                    trans_ampl_out11_2,                 & ! <-  the transmitted amplitude matrix (1,1) at all faces illuminated from this recursion
                                    trans_ampl_out12_2,                 & ! <-  the transmitted amplitude matrix (1,2) at all faces illuminated from this recursion
                                    trans_ampl_out21_2,                 & ! <-  the transmitted amplitude matrix (2,1) at all faces illuminated from this recursion
                                    trans_ampl_out22_2,                 & ! <-  the transmitted amplitude matrix (2,2) at all faces illuminated from this recursion
                                    vk71Int2, vk72Int2, vk73Int2,       & ! <-  the e-perp components of all faces illuminated from this recursion
                                    vk121Int2, vk122Int2, vk123Int2,    & ! <-  the transmitted propagation direction components at all faces illuminated from this recursion
                                    vk91Int2, vk92Int2, vk93Int2,       & ! <-  the reflected propagation direction components at all faces illuminated from this recursion
                                    InteractionInt2)                      ! <-  the interaction number of each illuminated face at this recursion

        ! step 5: stitch together some arrays for use in the next recursion of the beam loop

        if(.not. allocated(vk71Int3)) is_first_beam_back = .true. ! determine whether this is the first beam back (multithreading support)

        ! ! init v3 variables
        if(i .eq. 1) then ! if its the first beam back, we need to allocate first
            allocate(vk71Int3(1:size(FInt2,1),size(FInt2,2)))
            allocate(vk72Int3(1:size(FInt2,1),size(FInt2,2)))
            allocate(vk73Int3(1:size(FInt2,1),size(FInt2,2)))  
            allocate(refl_ampl_out11_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out12_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out21_3(1:size(FInt2,1),size(FInt2,2)))
            allocate(refl_ampl_out22_3(1:size(FInt2,1),size(FInt2,2)))  
            allocate(FInt3(1:size(FInt2,1),size(FInt2,2)))
            FInt3 = FInt2 ! keep track, ready to be concatenated to on the next loop
            vk71Int3 = vk71Int2
            vk72Int3 = vk72Int2
            vk73Int3 = vk73Int2
            refl_ampl_out11_3 = refl_ampl_out11_2
            refl_ampl_out12_3 = refl_ampl_out12_2
            refl_ampl_out21_3 = refl_ampl_out21_2
            refl_ampl_out22_3 = refl_ampl_out22_2
            propagationVectors3 = propagationVectors2       
        else
            call cat_int_var(FInt2, FInt3)
            call cat_real_var(vk71Int2, vk71Int3)
            call cat_real_var(vk72Int2, vk72Int3)
            call cat_real_var(vk73Int2, vk73Int3)
            call cat_complex_var(refl_ampl_out11_2, refl_ampl_out11_3)
            call cat_complex_var(refl_ampl_out12_2, refl_ampl_out12_3)
            call cat_complex_var(refl_ampl_out21_2, refl_ampl_out21_3)
            call cat_complex_var(refl_ampl_out22_2, refl_ampl_out22_3)
            call cat_prop(propagationVectors2, propagationVectors3) ! concatenate propagation vectors
            ! stop
        end if
    end do
end if

! stop

deallocate(F1Mapping)
deallocate(vk71Int)
deallocate(vk72Int)
deallocate(vk73Int)
deallocate(refl_ampl_out11Int)
deallocate(refl_ampl_out12Int)
deallocate(refl_ampl_out21Int)
deallocate(refl_ampl_out22Int)
deallocate(propagationVectors)

allocate(F1Mapping(1:size(FInt3,1),1:size(FInt3,2)))
allocate(vk71Int(1:size(vk71Int3,1),1:size(vk71Int3,2)))
allocate(vk72Int(1:size(vk72Int3,1),1:size(vk72Int3,2)))
allocate(vk73Int(1:size(vk73Int3,1),1:size(vk73Int3,2)))
allocate(refl_ampl_out11Int(1:size(refl_ampl_out11_3,1),1:size(refl_ampl_out11_3,2)))
allocate(refl_ampl_out12Int(1:size(refl_ampl_out12_3,1),1:size(refl_ampl_out12_3,2)))
allocate(refl_ampl_out21Int(1:size(refl_ampl_out21_3,1),1:size(refl_ampl_out21_3,2)))
allocate(refl_ampl_out22Int(1:size(refl_ampl_out22_3,1),1:size(refl_ampl_out22_3,2)))
allocate(propagationVectors(1:size(propagationVectors3,1),1:size(propagationVectors3,2),1:size(propagationVectors3,3)))

! init
F1Mapping = 0
vk71Int = 0
vk72Int = 0
vk73Int = 0
refl_ampl_out11Int = 0
refl_ampl_out12Int = 0
refl_ampl_out21Int = 0
refl_ampl_out22Int = 0
propagationVectors = 0

F1Mapping(1:size(FInt3,1),1:size(FInt3,2)) = FInt3(1:size(FInt3,1),1:size(FInt3,2))
vk71Int(1:size(vk71Int3,1),1:size(vk71Int3,2)) = vk71Int3(1:size(vk71Int3,1),1:size(vk71Int3,2))
vk72Int(1:size(vk72Int3,1),1:size(vk72Int3,2)) = vk72Int3(1:size(vk72Int3,1),1:size(vk72Int3,2))
vk73Int(1:size(vk73Int3,1),1:size(vk73Int3,2)) = vk73Int3(1:size(vk73Int3,1),1:size(vk73Int3,2))
refl_ampl_out11Int(1:size(refl_ampl_out11_3,1),1:size(refl_ampl_out11_3,2)) = refl_ampl_out11_3(1:size(refl_ampl_out11_3,1),1:size(refl_ampl_out11_3,2))
refl_ampl_out12Int(1:size(refl_ampl_out12_3,1),1:size(refl_ampl_out12_3,2)) = refl_ampl_out12_3(1:size(refl_ampl_out12_3,1),1:size(refl_ampl_out12_3,2))
refl_ampl_out21Int(1:size(refl_ampl_out21_3,1),1:size(refl_ampl_out21_3,2)) = refl_ampl_out21_3(1:size(refl_ampl_out21_3,1),1:size(refl_ampl_out21_3,2))
refl_ampl_out22Int(1:size(refl_ampl_out22_3,1),1:size(refl_ampl_out22_3,2)) = refl_ampl_out22_3(1:size(refl_ampl_out22_3,1),1:size(refl_ampl_out22_3,2))
propagationVectors(1:size(propagationVectors3,1),1:size(propagationVectors3,2),1:size(propagationVectors3,3)) = &
    propagationVectors3(1:size(propagationVectors3,1),1:size(propagationVectors3,2),1:size(propagationVectors3,3))

! stop

! finish = omp_get_wtime()

! print*,'=========='
! print'(A,f16.8,A)',"total time elapsed: ",finish-start," secs"

! stop

end subroutine

subroutine init_for_main_loop(Face1, trans_ampl_ps, F1Mapping, &
                              refl_ampl_out11Int, refl_ampl_out12Int, refl_ampl_out21Int, refl_ampl_out22Int, &
                              aperturePropagationVectors, propagationVectors, &
                              vk71, vk72, vk73, vk71Int, vk72Int, vk73Int)

! initialises the variables that go into the main beam loop

integer, dimension(:,:), allocatable, intent(out) :: F1Mapping ! the face indices of each row of variables to go into main loop
complex(8), dimension(:,:), allocatable, intent(out) :: refl_ampl_out11Int ! the amplitude matrix that goes into the main loop
complex(8), dimension(:,:), allocatable, intent(out) :: refl_ampl_out12Int ! needs renaming
complex(8), dimension(:,:), allocatable, intent(out) :: refl_ampl_out21Int
complex(8), dimension(:,:), allocatable, intent(out) :: refl_ampl_out22Int
complex(8), dimension(:,:,:), allocatable, intent(in) :: trans_ampl_ps ! transmitted amplitude matrices, including shadowed facets
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
real(8), dimension(:,:,:), allocatable, intent(out) :: propagationVectors ! propagation vector of beam emitted from each aperture
real(8), dimension(:,:), allocatable, intent(in) :: aperturePropagationVectors ! propagation vector of beam emitted from each aperture
real(8), dimension(:), allocatable, intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
real(8), dimension(:,:), allocatable, intent(out) :: vk71Int, vk72Int, vk73Int ! reflected e-perp vector

integer i

! allocations (over-allocated at the moment, can reduce later)
allocate(F1Mapping(1:size(Face1,1),1:1))
allocate(refl_ampl_out11Int(1:size(Face1,1),1:1))
allocate(refl_ampl_out12Int(1:size(Face1,1),1:1))
allocate(refl_ampl_out21Int(1:size(Face1,1),1:1))
allocate(refl_ampl_out22Int(1:size(Face1,1),1:1))
allocate(vk71Int(1:size(Face1,1),1:1))
allocate(vk72Int(1:size(Face1,1),1:1))
allocate(vk73Int(1:size(Face1,1),1:1))
allocate(propagationVectors(1:size(aperturePropagationVectors,1),1:3,1:1))

do i = 1, size(Face1,1) ! for each facet
    F1Mapping(i,1) = i ! record facet ID
    vk71Int(i,1) = vk71(i) ! record e-perp direction
    vk72Int(i,1) = vk72(i) ! record e-perp direction
    vk73Int(i,1) = vk73(i) ! record e-perp direction
    refl_ampl_out11Int(i,1) = trans_ampl_ps(1,1,i) ! record transmitted amplitude matrix (including shadow)
    refl_ampl_out12Int(i,1) = trans_ampl_ps(1,2,i)
    refl_ampl_out21Int(i,1) = trans_ampl_ps(2,1,i)
    refl_ampl_out22Int(i,1) = trans_ampl_ps(2,2,i)
end do

propagationVectors = 0 ! initialise
do i = 1, size(aperturePropagationVectors,1) ! for each aperture
    propagationVectors(i,1,1) = aperturePropagationVectors(i,1)
    propagationVectors(i,2,1) = aperturePropagationVectors(i,2)
    propagationVectors(i,3,1) = aperturePropagationVectors(i,3)
end do

end subroutine


subroutine add_refl_to_outbeam_tree(beam_outbeam_tree, beam_outbeam_tree_counter, isWithinBeam, &
                                    refl_ampl, vk91, vk92, vk93, vk71, vk72, vk73, apertures, interactionCounter)

! adds external reflection to outbeam tree

type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
integer, intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
logical, dimension(:), allocatable, intent(in) :: isWithinBeam ! whether each visible facet was within the illuminating beam
complex(8), dimension(:,:,:), allocatable, intent(in) :: refl_ampl ! reflected amplitude matrices
real(8), dimension(:), allocatable, intent(in) :: vk91, vk92, vk93 ! reflected prop vector from each facet
real(8), dimension(:), allocatable, intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
integer, intent(inout) :: interactionCounter ! counts the current number of interactions

integer i, j
integer, dimension(:), allocatable :: outbeam_apertures ! the aperture which each facet belongs to
integer, dimension(:), allocatable :: unique_apertures !

! add stuff to outbeam tree for the first external reflection
do i = 1, size(isWithinBeam,1) ! for each facet
    if(isWithinBeam(i)) then ! if the facet was within the illuminating beam
        beam_outbeam_tree_counter = beam_outbeam_tree_counter + 1
        beam_outbeam_tree(beam_outbeam_tree_counter)%FOut = i
        beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(1) = vk71(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(2) = vk72(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(3) = vk73(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(1) = vk91(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(2) = vk92(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(3) = vk93(i)
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(1) = 0
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(2) = 0
        beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(3) = -1
        beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(1:2,1:2) = refl_ampl(1:2,1:2,i)
    end if
end do

! also sort the interaction counter (ugh...)
! make an array to hold the apertures of outbeams
allocate(outbeam_apertures(1:beam_outbeam_tree_counter))
do i = 1, beam_outbeam_tree_counter
    outbeam_apertures(i) = apertures(beam_outbeam_tree(i)%FOut)
end do

! get the unique apertures
call unique_int(outbeam_apertures,unique_apertures)

! set interaction counter based on which of the unqiue illuminated apertures each outbeam facet belonged to
do i = 1, size(unique_apertures,1) ! for each illuminated aperture
    interactionCounter = interactionCounter + 1 ! update interaction counter
    do j = 1, beam_outbeam_tree_counter ! for each outbeam
        if(apertures(beam_outbeam_tree(j)%FOut) .eq. unique_apertures(i)) beam_outbeam_tree(j)%interactionOut = interactionCounter ! if facet of outbeam belonged to this aperture, assign interaction counter
    end do
end do


end subroutine

subroutine make_external_diff_outbeam_tree(new_in_ampl, ampl_diff, isVisible, ext_diff_outbeam_tree, vk71, vk72, vk73, verts, Face1)

! makes the outbeam tree using the visible facets from the first illumination
! can safely remove ampl_diff from the entire code at a later date

complex(8), dimension(:,:,:), allocatable, intent(in) :: new_in_ampl ! amplitude matrix after rotating into new scattering plane
complex(8), dimension(:,:,:), allocatable, intent(inout) :: ampl_diff ! external diffraction amplitude matrices
logical, dimension(:), allocatable, intent(in) :: isVisible ! whether each facet is visible as viewed in -z direction
type(outbeamtype), dimension(:), allocatable, intent(out) :: ext_diff_outbeam_tree
real(8), dimension(:), allocatable, intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs

integer i, num_outbeams_counter

num_outbeams_counter = 0 ! initialise
do i = 1, size(isVisible) ! for each facet
    if(isVisible(i)) num_outbeams_counter = num_outbeams_counter + 1 ! count the number of facets for external diffraction
end do

allocate(ext_diff_outbeam_tree(1:num_outbeams_counter)) ! allocate outbeam tree

num_outbeams_counter = 0 ! reset counter
do i = 1, size(isVisible) ! for each facet
    if(isVisible(i)) then 
        num_outbeams_counter = num_outbeams_counter + 1 ! count the number of facets for external diffraction
        ! save variables to outbeam tree
        ext_diff_outbeam_tree(num_outbeams_counter)%ampl(1:2,1:2) = new_in_ampl(1:2,1:2,i) ! amplitude matrix
        ext_diff_outbeam_tree(num_outbeams_counter)%vk7(1) = vk71(i) ! outgoing perpendicular e-field direction
        ext_diff_outbeam_tree(num_outbeams_counter)%vk7(2) = vk72(i)
        ext_diff_outbeam_tree(num_outbeams_counter)%vk7(3) = vk73(i)
        ext_diff_outbeam_tree(num_outbeams_counter)%verts(1:3,1:3) = transpose(verts(Face1(i,1:3),1:3)) ! x, y, and z components of each vertex
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_out(1) = 0 ! outgoing propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_out(2) = 0 ! outgoing propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_out(3) = -1 ! outgoing propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_in(1) = 0 ! incoming propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_in(2) = 0 ! incoming propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%prop_in(3) = -1 ! incoming propagation vector
        ext_diff_outbeam_tree(num_outbeams_counter)%FOut = i ! face ID from which the beam was emitted

        ! if(num_outbeams_counter .eq. 2341) print*,'num_outbeams_counter:',num_outbeams_counter,' -> face:',i

    end if
end do

ext_diff_outbeam_tree(num_outbeams_counter)%interactionOut = 0 ! interaction counter set to 0 for external diffraction


! print*,'number of visible facets (for external diffraction):',num_outbeams_counter 

ampl_diff = new_in_ampl ! save external diffraction amplitude matrices for later

end subroutine

subroutine applyFresnelMatrices(new_in_ampl, rperp, rpar, tperp, tpar, refl_ampl, trans_ampl)
 
complex(8), dimension(:,:,:), allocatable, intent(in) :: new_in_ampl ! amplitude matrix after rotating into new scattering plane
complex(8), dimension(:), allocatable, intent(in) :: rperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(in) :: rpar ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(in) :: tperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(in) :: tpar ! Fresnel coefficient
complex(8), dimension(:,:,:), allocatable, intent(inout) :: refl_ampl ! reflected amplitude matrices
complex(8), dimension(:,:,:), allocatable, intent(inout) :: trans_ampl ! transmitted amplitude matrices

integer i
complex(8) refl_matrix(1:2,1:2)
complex(8) trans_matrix(1:2,1:2)

! intitialise
refl_matrix = 0
trans_matrix = 0

do i = 1, size(rperp) ! for each facet
    refl_matrix(1,1) = rpar(i)
    refl_matrix(2,2) = rperp(i)
    trans_matrix(1,1) = tpar(i)
    trans_matrix(2,2) = tperp(i)
    refl_ampl(1:2,1:2,i) = matmul(refl_matrix(1:2,1:2),new_in_ampl(1:2,1:2,i))
    trans_ampl(1:2,1:2,i) = matmul(trans_matrix(1:2,1:2),new_in_ampl(1:2,1:2,i))
end do

end subroutine

subroutine rotateAmplitudeMatrices(rot_ampl,ampl_in,new_in_ampl)

real(8), dimension(:,:,:), allocatable :: rot_ampl ! rotation matrix for beams incident on eahc facet
complex(8), dimension(:,:,:), allocatable :: ampl_in ! amplitude matrix over the surface, for a specific recursion
complex(8), dimension(:,:,:), allocatable :: new_in_ampl ! amplitude matrix after rotating into new scattering plane

integer i

do i = 1, size(rot_ampl,3) ! for each facet
    new_in_ampl(1:2,1:2,i) = matmul(rot_ampl(1:2,1:2,i),ampl_in(1:2,1:2,i))
end do

end subroutine

subroutine getRotationMatrices(rot_ampl, vk71, vk72, vk73, ev11, ev12, ev13, ev31, ev32, ev33)

real(8), dimension(:), allocatable, intent(in) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
real(8), dimension(:,:,:), allocatable, intent(inout) :: rot_ampl ! rotation matrix for beams incident on eahc facet
real(8), intent(in) :: ev11, ev12, ev13 ! incident e-perp vector
real(8), intent(in) :: ev31, ev32, ev33 ! incident propagation vector

integer i

do i = 1, size(vk71,1) ! for each facet
    call getRotationMatrix(rot_ampl(1:2,1:2,i),vk71(i),vk72(i),vk73(i),ev11,ev12,ev13,ev31,ev32,ev33) ! get the rotation matrix
end do

end subroutine

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

subroutine getReflectionVectors(Norm, Face2, isShadow, apertureNormals, apertures, &
                                vk71, vk72, vk73, vk91, vk92, vk93)

! modified version of sr reflectionVectors
! returns perp e-field vector and propagation vectors

real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
logical, dimension(:), allocatable, intent(in) :: isShadow ! whether the facet was in shadow (down-facing but within the illuminating beam and part of an illuminated aperture)
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), dimension(:), allocatable, intent(inout) :: vk71, vk72, vk73 ! reflected e-perp vector from each facet
real(8), dimension(:), allocatable, intent(inout) :: vk91, vk92, vk93 ! reflected prop vector from each facet

real(8) ev11,ev12,ev13 ! incident e-perp vector
real(8) ev3(1:3) ! incident propagation vector
real(8) normal(1:3) ! surface normal
real(8) cross_temp(1:3) ! temporary cross product output
real(8) w, help1 ! temporary constants
real(8) nf ! normalisation factor
integer i
real(8), parameter :: thresh = -0.999999

ev11 = 1 ! take care with this in the future
ev12 = 0
ev13 = 0
ev3(1) = 0 ! same with this, although probably less care needed
ev3(2) = 0
ev3(3) = -1

do i = 1, size(Face2,1) ! for each face
    ! if(isShadow(i)) then ! if shadow, use aperture normal because Snell's Law unphysical
        normal(1:3) = apertureNormals(apertures(i),1:3)
        w = -abs(apertureNormals(apertures(i),3)) ! dot product of incidence vector with aperture normal
    ! else
    !     normal(1:3) = Norm(Face2(i),1:3)
    !     w = -abs(Norm(Face2(i),3)) ! dot product of incidence vector with face normal
    ! end if
    if(w .lt. thresh) then
        vk71(i) = 1
        vk72(i) = 0
        vk73(i) = 0
    end if
    help1 = -2*w
    vk91(i) = ev3(1) + help1*normal(1)
    vk92(i) = ev3(2) + help1*normal(2)
    vk93(i) = ev3(3) + help1*normal(3)
    nf = sqrt(vk91(i)**2 + vk92(i)**2 + vk93(i)**2)
    vk91(i) = vk91(i) / nf
    vk92(i) = vk92(i) / nf
    vk93(i) = vk93(i) / nf
    call cross(normal(1:3),ev3(1:3),cross_temp(1:3),.true.) ! call cross product with normalisation
    vk71(i) = cross_temp(1)
    vk72(i) = cross_temp(2)
    vk73(i) = cross_temp(3)
end do

end subroutine

subroutine getPropagationVectors(aperturePropagationVectors,sufficientlyIlluminated,apertureNormals,m)

real(8), dimension(:,:), allocatable, intent(inout) :: aperturePropagationVectors ! propagation vector of beam emitted from each aperture
logical, dimension(:), allocatable, intent(in) :: sufficientlyIlluminated
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
real(8), intent(in) :: m ! real part of the refractive index

integer i
real(8) theta_i ! incident angle
real(8) theta_t ! refracted angle
real(8) A, B ! coefficients to be determined

aperturePropagationVectors = 0 ! initialise

do i = 1, size(sufficientlyIlluminated,1) ! for each aperture
    if(sufficientlyIlluminated(i)) then ! if the aperture was sufficiently illuminated
        theta_i = acos(apertureNormals(i,3)) ! angle of incidence
        theta_t = asin(sin(theta_i)/m) ! angle of refraction
        A = sin(theta_i-theta_t)/sin(pi-theta_i)
        B = sin(theta_t)/sin(pi-theta_i)
        aperturePropagationVectors(i,1) = -A*apertureNormals(i,1)
        aperturePropagationVectors(i,2) = -A*apertureNormals(i,2)
        aperturePropagationVectors(i,3) = -B-A*apertureNormals(i,3)
        ! print'(A,I,A,f10.6,f10.6,f10.6)','aperture',i,' propagation vector: ',aperturePropagationVectors(i,1),aperturePropagationVectors(i,2),aperturePropagationVectors(i,3)
    end if
end do

end subroutine

subroutine getFresnel(Face1, isShadow, Norm, Face2, apertureNormals, apertures, rbi, rperp, rpar, tperp, tpar, ibi)

integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs
logical, dimension(:), allocatable, intent(in) :: isShadow ! whether the facet was in shadow (down-facing but within the illuminating beam and part of an illuminated aperture)
real(8), dimension(:,:), allocatable, intent(in) :: Norm ! face normals
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals ! the normal of each aperture
integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
real(8), intent(in) :: rbi ! real part of the refractive index
real(8), intent(in) :: ibi ! real part of the refractive index
complex(8), dimension(:), allocatable, intent(inout) :: rperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(inout) :: rpar ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(inout) :: tperp ! Fresnel coefficient
complex(8), dimension(:), allocatable, intent(inout) :: tpar ! Fresnel coefficient

integer i
real(8) normal ! z component of face normal
real(8) theta_i ! incident theta in rads
real(8) theta_t ! refracted theta in rads
complex(8) m ! complex refractive index

m = cmplx(rbi,ibi,8)

! for external interactions
do i = 1, size(Face1,1) ! for each face
    if(isShadow(i)) then ! if the facet is in the shadow
        normal = apertureNormals(apertures(i),3) ! just use aperture normal for shadow region (Snells law for down facing doesnt make sense)
    else
        normal = Norm(Face2(i),3)
    end if
    theta_i = acos(normal) ! get incident angle
    theta_t = asin(sin(theta_i)/m) ! get transmitted angle using Snell
    rperp(i) = (cos(theta_i) - m*cos(theta_t))/(cos(theta_i) + m*cos(theta_t))
    tperp(i) = (2*cos(theta_i))/(cos(theta_i) + m*cos(theta_t))
    rpar(i) = (m*cos(theta_i) - cos(theta_t))/(cos(theta_t) + m*cos(theta_i))
    tpar(i) = (2*cos(theta_i))/(cos(theta_t) + m*cos(theta_i))
end do  


end subroutine

subroutine propagateBeam(ampl_in, ampl_in_ps, ampl_beam, isWithinBeam, isWithinBeam_ps, beamIDs, beamIDs_ps, waveno, distances, distances_ps, Face1)

complex(8), dimension(:,:,:), allocatable, intent(inout) :: ampl_in ! amplitude matrix over the surface, for a specific recursion
complex(8), dimension(:,:,:), allocatable, intent(inout) :: ampl_in_ps ! amplitude matrix over the surface, for a specific recursion
complex(8), allocatable, dimension(:,:,:), intent(in) :: ampl_beam ! amplitude matrix of incident beam
logical, dimension(:), allocatable, intent(in) :: isWithinBeam ! whether each visible facet was within the illuminating beam
logical, dimension(:), allocatable, intent(in) :: isWithinBeam_ps ! whether each visible facet was within the illuminating beam, including down-facing facets of illuminated apertures
integer, dimension(:), allocatable, intent(in) :: beamIDs ! the beamF1 ID of the face which illuminated this facet
integer, dimension(:), allocatable, intent(in) :: beamIDs_ps ! the beamF1 ID of the face which illuminated this facet, including down-facing facets of illuminated aperture
real(8), intent(in) :: waveno ! wavenumber in vacuum
real(8), dimension(:), allocatable, intent(in) :: distances ! distances to each illuminated face from the illuminating face ID given by beamIDs
real(8), dimension(:), allocatable, intent(in) :: distances_ps ! distances to each illuminated face from the illuminating face ID given by beamIDs, including down-facing facets of illuminated apertures
integer(8), dimension(:,:), allocatable, intent(in) :: Face1 ! face vertex IDs

integer i, blockingID

do i = 1, size(Face1,1) ! for each face
    if(isWithinBeam(i)) then ! if the face was within the beam
        blockingID = beamIDs(i)
        ! ampl_in(1:2,1:2,i) = ampl_beam(1:2,1:2,blockingID)*exp(1i*waveno*distances(i))
        ampl_in(1:2,1:2,i) = ampl_beam(1:2,1:2,blockingID)*exp2cmplx(waveno*distances(i))
        ampl_in_ps(1:2,1:2,i) = ampl_beam(1:2,1:2,blockingID)*exp2cmplx(waveno*distances(i)) ! fix 12/08/23
        ! ampl_in(1,1,i) = ampl_beam(1,1,blockingID)*exp2cmplx(waveno*distances(i))
        ! if (i .eq. 6282) print*,'distances(i)',distances(i)
        ! if (i .eq. 6282) print*,'exp2cmplx(waveno*distances(i))',exp2cmplx(waveno*distances(i))
    end if
    if(isWithinBeam_ps(i)) then ! if the face was within the beam
        blockingID = beamIDs_ps(i)
        ! ampl_in_ps(1:2,1:2,i) = ampl_beam(1:2,1:2,blockingID)*exp(1i*waveno*distances_ps(i))
        ampl_in_ps(1:2,1:2,i) = ampl_beam(1:2,1:2,blockingID)*exp2cmplx(waveno*distances_ps(i))
    end if
    
end do

end subroutine

subroutine getSufficientlyIlluminated(illuminatedApertureAreas,threshold,sufficientlyIlluminated)

real(8), dimension(:), allocatable, intent(in) :: illuminatedApertureAreas ! the llluminated area of each aperture
logical, dimension(:), allocatable, intent(inout) :: sufficientlyIlluminated
real(8), intent(in) :: threshold ! minimum area of illumination per aperture to create new beam

integer i

sufficientlyIlluminated = .false.
do i = 1, size(illuminatedApertureAreas,1)
    if(illuminatedApertureAreas(i) .gt. threshold) sufficientlyIlluminated(i) = .true.
    ! print*,'aperture:',i,'sufficiently Illuminated? ',sufficientlyIlluminated(i)
end do

end subroutine

subroutine getShadow(isWithinBeam,isWithinBeam_ps,isShadow)

logical, dimension(:), allocatable, intent(in) :: isWithinBeam, isWithinBeam_ps
logical, dimension(:), allocatable, intent(inout) :: isShadow

integer i

isShadow = .false.
do i = 1, size(isWithinBeam,1) ! for each facet
    if(isWithinBeam_ps(i) .and. isWithinBeam(i) .eqv. .false.) isShadow(i) = .true. ! determine if it was visible but in the shadow
end do

end subroutine

subroutine getThreshold(la,threshold,job_params)

real(8), intent(out) :: threshold
real(8), intent(in) :: la
type(job_parameters_type), intent(in) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details

! print*,'========== start sr getThreshold'

threshold = (2*la)**2
threshold = 0
if(job_params%debug >= 1) then
    print*,'threshold: ',threshold
end if

! print*,'========== end sr getThreshold'

end subroutine

subroutine getIlluminatedGeoCrossSection(faceAreas, isWithinBeam, Norm, Face2, illuminatedGeoCrossSection)

real(8), dimension(:), allocatable, intent(in) :: faceAreas
logical, dimension(:), allocatable, intent(in) :: isWithinBeam
real(8), dimension(:,:), allocatable, intent(in) :: Norm
integer, dimension(:), allocatable, intent(in) :: Face2
real(8), intent(out) :: illuminatedGeoCrossSection

integer i

! print*,'========== start sr getIlluminatedGeoCrossSection'

illuminatedGeoCrossSection = 0
do i = 1,size(faceAreas,1) ! for each facet
    if(isWithinBeam(i)) illuminatedGeoCrossSection = illuminatedGeoCrossSection + faceAreas(i)*Norm(Face2(i),3)
end do

! print*,'illuminated geometric cross section: ',illuminatedGeoCrossSection

! print*,'========== end sr getIlluminatedGeoCrossSection'

end subroutine

subroutine initApertures(apertures, Norm, Face2, Midpoints, apertureNormals, apertureMidpoints, apertureAreas, illuminatedApertureAreas, sufficientlyIlluminated, &
                         aperturePropagationVectors)

integer(8), dimension(:), allocatable, intent(in) :: apertures
real(8), dimension(:,:), allocatable, intent(in) :: Norm
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints
real(8), dimension(:,:), allocatable, intent(out) :: aperturePropagationVectors
integer, dimension(:), allocatable, intent(in) :: Face2
real(8), dimension(:,:), allocatable, intent(out) :: apertureNormals
real(8), dimension(:,:), allocatable, intent(out) :: apertureMidpoints
real(8), dimension(:), allocatable, intent(out) :: apertureAreas, illuminatedApertureAreas
logical, dimension(:), allocatable, intent(out) :: sufficientlyIlluminated

integer noFaces, i, j, noApertures, face_counter
real(8) nf

noFaces = size(Face2,1)

!print*,'number of faces: ',noFaces

noApertures = maxval(apertures)

allocate(apertureNormals(1:noApertures,1:3))
allocate(apertureMidpoints(1:noApertures,1:3))
allocate(aperturePropagationVectors(1:noApertures,1:3))
allocate(apertureAreas(1:noApertures))
allocate(illuminatedApertureAreas(1:noApertures))
allocate(sufficientlyIlluminated(1:noApertures))

! initialise
apertureNormals = 0
apertureMidpoints = 0

!print*,'number of apertures: ',noApertures

do i = 1, noApertures ! for each aperture
    face_counter = 0 ! counts the number of faces in this aperture
    do j = 1, noFaces
        if(apertures(j) .eq. i) then ! if the face belongs to this aperture
        face_counter = face_counter + 1 ! update face counter
        apertureNormals(i,1:3) = apertureNormals(i,1:3) + Norm(Face2(j),1:3) ! add the normal of this face to the total for this aperture
        apertureMidpoints(i,1:3) = apertureMidpoints(i,1:3) + Midpoints(j,1:3)
        end if
    end do
    apertureMidpoints(i,1:3) = apertureMidpoints(i,1:3) / face_counter ! average
    nf = sqrt(apertureNormals(i,1)**2 + apertureNormals(i,2)**2 + apertureNormals(i,3)**2)
    apertureNormals(i,1:3) = apertureNormals(i,1:3) / nf ! normalise normals
end do

! do i = 1, noApertures
!     print'(a,i2,a,f10.6,f10.6,f10.6)','aperture ',i,' normal: ',apertureNormals(i,1),apertureNormals(i,2),apertureNormals(i,3)
! end do
! do i = 1, noApertures
!     print'(a,i2,a,f10.6,f10.6,f10.6)','aperture ',i,' midpoint: ',apertureMidpoints(i,1),apertureMidpoints(i,2),apertureMidpoints(i,3)
! end do

end subroutine

subroutine getApertureAreas(apertures, isWithinBeam, apertureAreas, illuminatedApertureAreas, apertureNormals, faceAreas)

integer(8), dimension(:), allocatable, intent(in) :: apertures
logical, dimension(:), allocatable, intent(in) :: isWithinBeam
real(8), dimension(:), allocatable, intent(inout) :: apertureAreas, illuminatedApertureAreas
real(8), dimension(:,:), allocatable, intent(in) :: apertureNormals
real(8), dimension(:), allocatable, intent(in) :: faceAreas

logical, dimension(:), allocatable :: logical_array
integer i, j

allocate(logical_array(1:size(apertures,1))) ! logical array with number of rows equal to number of facets
illuminatedApertureAreas = 0 ! initialise
apertureAreas = 0 ! initialise

do i = 1, size(apertureAreas,1) ! for each aperture
    logical_array = .false. ! initialise
    do j = 1, size(apertures,1)
        if(apertures(j) .eq. i .and. isWithinBeam(j)) logical_array(j) = .true. ! if facet belongs to this aperture and was within the beam, record for later
    end do
    if(apertureNormals(i,3) .gt. 0) then ! if aperture faces upwards
        do j = 1, size(apertures,1) ! for each facet
            if(apertures(j) .eq. i .and. isWithinBeam(j)) illuminatedApertureAreas(i) = illuminatedApertureAreas(i) + faceAreas(j)    
        end do
    end if
    do j = 1, size(apertures,1) ! for each facet
        if(apertures(j) .eq. i) then
            apertureAreas(i) = apertureAreas(i) + faceAreas(j)
        end if
    end do
    !print*,'aperture: ',i,' Illuminated area: ',illuminatedApertureAreas(i)
    !print*,'aperture: ',i,' Total area: ',apertureAreas(i)
end do
    
end subroutine

subroutine findWithinBeam(Face1, Midpoints, isVisible, beamV, beamF1, beamN, beamF2, beamMidpoints, isWithinBeam, distances, beamIDs)

integer(8), dimension(:,:), allocatable, intent(in) :: Face1
integer(8), dimension(:,:), allocatable, intent(in) :: beamF1
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints, beamMidpoints
logical, dimension(:), allocatable, intent(in) :: isVisible
real(8), dimension(:,:), allocatable, intent(in) :: beamV, beamN
integer(8), dimension(:), allocatable, intent(in) :: beamF2
logical, dimension(:), allocatable, intent(inout) :: isWithinBeam
real(8), dimension(:), allocatable, intent(inout) :: distances
integer, dimension(:), allocatable, intent(inout) :: beamIDs

real(8) start, finish
integer, dimension(:), allocatable :: upFacingIndices, indicesToLookAt
integer i, j, k, l, m, numUpFacingIndices, numIndicesToLookAt, noPoints
logical, dimension(:), allocatable :: areBlocking
real(8), dimension(:,:), allocatable :: vecAs, vecBs, edgeNormals, edgeVectors
real(8), dimension(:), allocatable :: AdotNs, NdotKs, vecs, edgeChecks, nfs

! print*,'========== start sr findWithinBeam'
call CPU_TIME(start)

isWithinBeam = .false. ! start by assuming that all facets are not within beam
distances = 0 ! distances from beam to each visible face
beamIDs = 0 ! beam face IDs which blocked each face

! allocations
allocate(upFacingIndices(1:size(Face1,1)))
allocate(indicesToLookAt(1:size(beamF1,1)))
allocate(areBlocking(1:size(beamF1,1)))
allocate(vecAs(1:size(beamF1,1),1:3))
allocate(AdotNs(1:size(beamF1,1)))
allocate(NdotKs(1:size(Face1,1)))
allocate(vecs(1:size(Face1,1)))
allocate(edgeVectors(1:size(beamF1,1),1:2))
allocate(vecBs(1:size(beamF1,1),1:2))
allocate(edgeNormals(1:size(beamF1,1),1:2))
allocate(edgeChecks(1:size(beamF1,1)))
allocate(nfs(1:size(beamF1,1)))

numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
    if(isVisible(i)) then
        numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
        upFacingIndices(numUpFacingIndices) = i
    end if        
end do
numIndicesToLookAt = size(beamF1,1)
do i = 1, numIndicesToLookAt
    indicesToLookAt(i) = i
end do

    do i = 1, numUpFacingIndices ! for each visible, upfacing face within this bounding box
        k = upFacingIndices(i) ! get face ID
        areBlocking(1:numIndicesToLookAt) = .true. ! start by assuming that this beam face blocks the facet
        vecAs(1:numIndicesToLookAt,1) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),1) - Midpoints(k,1) ! vectors from midpoint of face i to midpoints of faces
        vecAs(1:numIndicesToLookAt,2) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),2) - Midpoints(k,2)
        vecAs(1:numIndicesToLookAt,3) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),3) - Midpoints(k,3)
        AdotNs(1:numIndicesToLookAt) =  vecAs(1:numIndicesToLookAt,1)*beamN(beamF2(indicesToLookAt(1:numIndicesToLookAt)),1) + & ! dot products of face normals with vector A (above)
                                        vecAs(1:numIndicesToLookAt,2)*beamN(beamF2(indicesToLookAt(1:numIndicesToLookAt)),2) + &
                                        vecAs(1:numIndicesToLookAt,3)*beamN(beamF2(indicesToLookAt(1:numIndicesToLookAt)),3)
        NdotKs(1:numIndicesToLookAt) = -beamN(beamF2(indicesToLookAt(1:numIndicesToLookAt)),3) ! dot product of face normals with incident beam direction (z-axis)
        vecs(1:numIndicesToLookAt) = -AdotNs(1:numIndicesToLookAt) / NdotKs(1:numIndicesToLookAt)
        do j = 1, numIndicesToLookAt 
            if(vecs(j) .lt. 0) areBlocking(j) = .false. ! if faces are behind, they are not blocking
        end do    
        noPoints = size(beamF1,2)
        do m = 1, noPoints ! for each vertex in blocking face
            ! compute vector from this vertex to the next vertex of the blocking facet
            if(m .eq. noPoints) then
                edgeVectors(1:numIndicesToLookAt,1) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),1),1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),1),2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            else
                edgeVectors(1:numIndicesToLookAt,1) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m+1),1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m+1),2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            end if
            vecBs(1:numIndicesToLookAt,1) = Midpoints(k,1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1) ! get the vector from the vertex of the faces to the midpoint of face i
            vecBs(1:numIndicesToLookAt,2) = Midpoints(k,2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            edgeNormals(1:numIndicesToLookAt,1) = -edgeVectors(1:numIndicesToLookAt,2) ! cross product of edge vector with reverse beam direction
            edgeNormals(1:numIndicesToLookAt,2) = edgeVectors(1:numIndicesToLookAt,1)
            nfs(1:numIndicesToLookAt) = sqrt(edgeNormals(1:numIndicesToLookAt,1)**2 + edgeNormals(1:numIndicesToLookAt,2)**2) ! normalisation factor
            edgeNormals(1:numIndicesToLookAt,1) = edgeNormals(1:numIndicesToLookAt,1) / nfs(1:numIndicesToLookAt) ! normalise
            edgeNormals(1:numIndicesToLookAt,2) = edgeNormals(1:numIndicesToLookAt,2) / nfs(1:numIndicesToLookAt)
            edgeChecks(1:numIndicesToLookAt) = -(-vecBs(1:numIndicesToLookAt,1)*edgeNormals(1:numIndicesToLookAt,1) - vecBs(1:numIndicesToLookAt,2)*edgeNormals(1:numIndicesToLookAt,2))
            do l = 1,numIndicesToLookAt
                if(edgeChecks(l) .gt. 0) areBlocking(l) = .false. ! if its not "inside" the edge vector, face i has failed the check and is therefore not blocked by the face
            end do
            !print*,'areBlocking(1:numIndicesToLookAt)',areBlocking(1:numIndicesToLookAt)
            
        end do
        if(any(areBlocking(1:numIndicesToLookAt))) then
            isWithinBeam(k) = .true. ! visibility of up-facing facets
            distances(k) = vecs(findloc(areBlocking(1:numIndicesToLookAt),.true.,1))
            beamIDs(k) = indicesToLookat(findloc(areBlocking(1:numIndicesToLookAt),.true.,1))
        end if
        
    
    end do

call CPU_TIME(finish)
! print'(A,f12.8,A)',"facets within beam found in: ",finish-start," secs"
! print*,'========== end sr findWithinBeam'

end subroutine
    
subroutine beam_aligned_bounding_boxes(verts, boundingBoxV, boundingBoxF)

real(8), dimension(:,:), allocatable, intent(in) :: verts
real(8), allocatable, dimension(:,:), intent(out) :: boundingBoxV
integer(8), allocatable, dimension(:,:), intent(out) :: boundingBoxF

real(8) min_x, min_y, max_x, max_y, min_z, max_z
integer(8), parameter :: bounding_box_x_dim = 8 ! bounding box x dimension
integer(8), parameter :: bounding_box_y_dim = 8 ! bounding box y dimension
real(8), parameter :: fac = 1.1
real(8), dimension(1:bounding_box_x_dim+1, 1:bounding_box_y_dim+1) :: xvals, yvals ! bounding box x and y vertices
integer(8) i,j, vert_counter, face_counter

! get the max xyz and min xy coordinates
max_x = maxval(verts(1:size(verts,1),1))
min_x = minval(verts(1:size(verts,1),1))
max_y = maxval(verts(1:size(verts,1),2))
min_y = minval(verts(1:size(verts,1),2))
max_z = maxval(verts(1:size(verts,1),3))
min_z = minval(verts(1:size(verts,1),3))

!print*,'making bounding boxes'

! allocate arrays to hold outer bounding box
allocate(boundingBoxV(1:(bounding_box_x_dim+1)*(bounding_box_y_dim+1),1:3)) ! vertices, xyz coords
allocate(boundingBoxF(1:(bounding_box_x_dim*bounding_box_y_dim),1:4)) ! facets, 4 vertices

!print*,'size(boundingBoxV,1)',size(boundingBoxV,1)
!print*,'boundingBoxFSize',boundingBoxFSize

! make the bounding box vertices (meshgrid style)
boundingBoxV = 0 ! initialise
vert_counter = 0
do i = 1, bounding_box_x_dim + 1
    do j = 1, bounding_box_y_dim + 1
        xvals(i,j) = fac*(min_x + (max_x - min_x)/bounding_box_x_dim*(i-1))
        yvals(i,j) = fac*(min_y + (max_y - min_y)/bounding_box_x_dim*(j-1))
        !print*,'xvals(i,j)',xvals(i,j)
        !print*,'yvals(i,j)',yvals(i,j)
        vert_counter = vert_counter + 1
        boundingBoxV(vert_counter,1) = xvals(i,j)
        boundingBoxV(vert_counter,2) = yvals(i,j)
    end do
end do

!print*,'total vertices in vertex array',vert_counter

! make vertex and face arrays using the bottom left vertex as starting point (moving clockwise)
face_counter = 0
do i = 1, bounding_box_x_dim
    do j = 1, bounding_box_y_dim
        face_counter = face_counter + 1
        boundingBoxF(face_counter,1) = (i-1)*(bounding_box_y_dim+1) + j
        boundingBoxF(face_counter,2) = (i-1)*(bounding_box_y_dim+1) + j + 1
        boundingBoxF(face_counter,3) = i*(bounding_box_y_dim+1) + j + 1
        boundingBoxF(face_counter,4) = i*(bounding_box_y_dim+1) + j
    end do
end do

!print*,'bounding box min x',min_x*fac
!print*,'bounding box max x',max_x*fac
!print*,'bounding box min y',min_y*fac
!print*,'bounding box max y',max_y*fac
!
!! print vertex array
!do i = 1,vert_counter
!    do j = 1,3
!        print'(A,I5,f8.4,f8.4,f8.4)','verts',i,boundingBoxV(i,1),boundingBoxV(i,2),boundingBoxV(i,3)
!    end do
!end do
! print face array
!do i = 1,face_counter
!    print'(A,I5,I5,I5,I5,I5)','face',i,boundingBoxF(i,1),boundingBoxF(i,2),boundingBoxF(i,3),boundingBoxF(i,4)
!end do

end subroutine

subroutine findVisibleFacets(verts, Face1, Norm, Face2, midPoints, isVisible, isVisiblePlusShadows, apertures, num_face_vert)

! subroutine findVisibleFacets finds the external visible facets of the particle as viewed in the -z direction
! works for concave particles
! accuracy increases with mesh discretisation
! only difference between matlab version is that some facets were asigned to fuzzy bounding boxes differently
! may need to increase fac to improve consistency

real(8), dimension(:,:), allocatable, intent(in) :: verts, Norm, Midpoints
integer(8), dimension(:,:), allocatable, intent(in) :: Face1
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
logical, dimension(:), allocatable, intent(inout) :: isVisible, isVisiblePlusShadows
integer(8), dimension(:), allocatable, intent(in) :: apertures
integer(8), dimension(:), allocatable, intent(in) :: num_face_vert ! number of vertices in each face

real(8), allocatable, dimension(:,:) :: boundingBoxV
integer(8), allocatable, dimension(:,:) :: boundingBoxF
integer(8) boundingBoxFSize
real(8), dimension(:,:), allocatable :: boundingBoxMidpoints ! unique vertices, face vertex IDs, face normals, face midpoints
real(8), dimension(:), allocatable :: boundingBoxFaceAreas
integer(8), dimension(:), allocatable :: boundingBoxNumFaceVert
integer, dimension(:), allocatable :: F3 ! bounding box IDs
integer, dimension(:,:), allocatable :: F4 ! fuzzy bounding box IDs
real(8), dimension(:), allocatable :: distanceToBB, distanceToFuzzyBB
integer i, j, k, l, m, BB, numUpFacingIndices, numIndicesToLookAt, noPoints
integer, dimension(:), allocatable :: upFacingIndices, indicesToLookAt
logical, dimension(:), allocatable :: areBlocking, areBlocking2
real(8), dimension(:,:), allocatable :: vecAs, edgeVectors, vecBs, edgeNormals
real(8), dimension(:), allocatable :: AdotNs, NdotKs, vecs, edgeChecks, nfs
real(8) start, finish
logical, dimension(:), allocatable :: visibleApertures
integer numApertures
logical isApertureVisible
integer num_verts, num_verts_max

! print*,'========== start sr findVisibleFacets'
call CPU_TIME(start)

!print*,'bounding box x dimension: ',bounding_box_x_dim
!print*,'bounding box y dimension: ',bounding_box_y_dim
!print*,'is isVisible allocated? ',allocated(isVisible)

! ### beam-aligned bounding boxes ###
!print*,'making beam-aligned bounding boxes'
! use the current crystal vertices to create some bounding boxes in x-y plane
call beam_aligned_bounding_boxes(verts, boundingBoxV, boundingBoxF)

allocate(boundingBoxNumFaceVert(1:size(boundingBoxF,1)))
boundingBoxNumFaceVert = 4

! compute bounding box midpoints
call midPointsAndAreas(boundingBoxF, boundingBoxV, boundingBoxMidpoints, boundingBoxFaceAreas, boundingBoxNumFaceVert)

! check midpoints
!do i = 1,size(boundingBoxMidpoints,1)
!    print'(A,f10.6,f10.6)','boundingBoxMidpoints',boundingBoxMidpoints(i,1),boundingBoxMidpoints(i,2)
!end do

num_verts_max = maxval(num_face_vert)
! print*,'max number of vertices: ', num_verts_max

allocate(F3(1:size(Face1,1))) ! array to hold index of bounding box that each face belongs to
allocate(F4(1:size(Face1,1),1:num_verts_max)) ! array to hold index of fuzzy bounding box that each face belongs to
boundingBoxFSize = size(boundingBoxF,1) ! number of bounding box faces
allocate(distanceToBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box
allocate(distanceToFuzzyBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box

!print*,'Number of bounding boxes: ',boundingBoxFSize

! find which bounding box each vertex belongs to
do i = 1, size(Face1,1) ! for each face
    num_verts = num_face_vert(i)
    ! print*,'this face had ',num_verts,' vertices'
    distanceToBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - Midpoints(i,1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - Midpoints(i,2))**2) ! distance to each bb
    !do j = 1,boundingBoxFSize
    !    print*,'dist to box: ',distanceToBB(j)
    !end do
    F3(i) = minloc(distanceToBB,1) ! record which bounding box midpoint this facet was closest to
    !print'(A,i6,A,f10.6)','i',i,'dist to closest bb: ',distanceToBB(F3(i))
    !print'(A,f10.6,f10.6)','midpoints of closest bb: ',boundingBoxMidpoints(F3(i),1),boundingBoxMidpoints(F3(i),2)
    do j = 1, num_verts ! for each vertex
        distanceToFuzzyBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - verts(Face1(i,j),1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - verts(Face1(i,j),2))**2) ! distance to each bb
        F4(i,j) = minloc(distanceToFuzzyBB,1) ! record which bounding box midpoint this facet vertex was closest to
    end do 
end do

!do i = 1, size(Face1,1)
!    print*,'i',i,'F3(i): ',F3(i)
!end do


! next step is the loop for visible facets
isVisible = .true. ! start by assuming that all facets are unblocked and visible to the beam
isVisiblePlusShadows = .true. ! start by assuming that all facets are unblocked and visible to the beam
do i = 1, size(Face1,1)
    if(Norm(Face2(i),3) < 0.01) isVisible(i) = .false. ! down facing-facets are not visible (for now)
end do

! allocations
allocate(upFacingIndices(1:size(Face1,1)))
allocate(indicesToLookAt(1:size(Face1,1)))
allocate(areBlocking(1:size(Face1,1)))
allocate(areBlocking2(1:size(Face1,1)))
allocate(vecAs(1:size(Face1,1),1:3))
allocate(edgeVectors(1:size(Face1,1),1:2))
allocate(vecBs(1:size(Face1,1),1:2))
allocate(edgeNormals(1:size(Face1,1),1:2))
allocate(AdotNs(1:size(Face1,1)))
allocate(NdotKs(1:size(Face1,1)))
allocate(vecs(1:size(Face1,1)))
allocate(edgeChecks(1:size(Face1,1)))
allocate(nfs(1:size(Face1,1)))

do BB = 1, boundingBoxFSize ! for each bounding box
!do BB = 8, 8 ! for each bounding box
    !print*,'care, only using first BB atm'
    numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        !if(F3(i) .eq. BB .and. Norm(Face2(i),3) > -0.01) then
        if(F3(i) .eq. BB) then
            numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
            upFacingIndices(numUpFacingIndices) = i
        end if        
    end do
    numIndicesToLookAt = 0 ! make a counter to count the number of faces within fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        if(any(F4(i,1:3) .eq. BB,1) .and. Norm(Face2(i),3) > -0.01) then
            numIndicesToLookAt = numIndicesToLookAt + 1 ! count number of faces within non-fuzzy bounding box
            indicesToLookAt(numIndicesToLookAt) = i ! record ID of visible facet within fuzzy bounding box
        end if            
    end do
    !print*,'number of visible faces within this BB', numUpFacingIndices
    !print*,'visible face IDs within this BB',upFacingIndices(1:numUpFacingIndices)
    !print*,'number of visible faces within this fuzzy BB', numIndicesToLookAt
    !print*,'visible face IDs within this fuzzy BB',indicesToLookAt(1:numIndicesToLookAt)
    
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),1)
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),2)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),1)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),2)
    !print*,'size(Midpoints(i,1:3),1)',size(Midpoints(i,1:3),1)
    !print*,'size(Midpoints(i,1:3),2)',size(Midpoints(i,1:3),2)
    !stop
    do i = 1, numUpFacingIndices ! for each visible, upfacing face within this bounding box
        k = upFacingIndices(i) ! get face ID
        areBlocking(1:numIndicesToLookAt) = .true. ! logical array containing whether each face blocks face i
        vecAs(1:numIndicesToLookAt,1) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),1) - Midpoints(k,1) ! vectors from midpoint of face i to midpoints of faces
        vecAs(1:numIndicesToLookAt,2) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),2) - Midpoints(k,2)
        vecAs(1:numIndicesToLookAt,3) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),3) - Midpoints(k,3)
        AdotNs(1:numIndicesToLookAt) =  vecAs(1:numIndicesToLookAt,1)*Norm(Face2(indicesToLookAt(1:numIndicesToLookAt)),1) + & ! dot products of face normals with vector A (above)
                                        vecAs(1:numIndicesToLookAt,2)*Norm(Face2(indicesToLookAt(1:numIndicesToLookAt)),2) + &
                                        vecAs(1:numIndicesToLookAt,3)*Norm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3)
        NdotKs(1:numIndicesToLookAt) = -Norm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3) ! dot product of face normals with incident beam direction (z-axis)
        vecs(1:numIndicesToLookAt) = -AdotNs(1:numIndicesToLookAt) / NdotKs(1:numIndicesToLookAt)
        do j = 1, numIndicesToLookAt 
            if(vecs(j) .lt. 0) areBlocking(j) = .false. ! if faces are behind, they are not blocking
            if (indicesToLookAt(j) .eq. k) areBlocking(j) = .false. ! cant block itself
        end do
        
        noPoints = 3 ! number of vertices in blocking face
        do m = 1, noPoints ! for each vertex in blocking face
            ! compute vector from this vertex to the next vertex of the blocking facet
            if(m .eq. noPoints) then
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            else
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            end if
            vecBs(1:numIndicesToLookAt,1) = Midpoints(k,1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1) ! get the vector from the vertex of the faces to the midpoint of face i
            vecBs(1:numIndicesToLookAt,2) = Midpoints(k,2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            edgeNormals(1:numIndicesToLookAt,1) = -edgeVectors(1:numIndicesToLookAt,2) ! cross product of edge vector with reverse beam direction
            edgeNormals(1:numIndicesToLookAt,2) = edgeVectors(1:numIndicesToLookAt,1)
            nfs(1:numIndicesToLookAt) = sqrt(edgeNormals(1:numIndicesToLookAt,1)**2 + edgeNormals(1:numIndicesToLookAt,2)**2) ! normalisation factor
            edgeNormals(1:numIndicesToLookAt,1) = edgeNormals(1:numIndicesToLookAt,1) / nfs(1:numIndicesToLookAt) ! normalise
            edgeNormals(1:numIndicesToLookAt,2) = edgeNormals(1:numIndicesToLookAt,2) / nfs(1:numIndicesToLookAt)
            edgeChecks(1:numIndicesToLookAt) = -vecBs(1:numIndicesToLookAt,1)*edgeNormals(1:numIndicesToLookAt,1) - vecBs(1:numIndicesToLookAt,2)*edgeNormals(1:numIndicesToLookAt,2)
            do l = 1,numIndicesToLookAt
                if(edgeChecks(l) .gt. 0) areBlocking(l) = .false. ! if its not "inside" the edge vector, face i has failed the check and is therefore not blocked by the face
            end do
            !print*,'areBlocking(1:numIndicesToLookAt)',areBlocking(1:numIndicesToLookAt)
            
        end do
        ! shadow region stuff
        areBlocking2(1:numIndicesToLookAt) = areBlocking(1:numIndicesToLookAt) ! shadow region facets will have the same blockers as non-shadow region, but cannot be in the same aperture as the facet
        do m = 1, numIndicesToLookAt ! looping through all blockers
            if(apertures(indicesToLookat(m)) .eq. apertures(k)) areBlocking2(m) = .false. ! set as non-blocking if part of same aperture
        end do
        
        if(any(areBlocking(1:numIndicesToLookAt))) isVisible(k) = .false. ! visibility of up-facing facets
        if(any(areBlocking2(1:numIndicesToLookAt))) isVisiblePlusShadows(k) = .false. ! visibility of all facets of illuminated apertures, including shadow region
        
        !if(k.eq.7573) then ! debugging statement
        !print*,'looking for facets which block facet: ',k
        !print*,'this facet is in BB: ',BB
        !print*,'this BB has midpoints:',boundingBoxMidpoints(BB,1),boundingBoxMidpoints(BB,2)
        !print*,'indicesToLookAt(1:numIndicesToLookAt)',indicesToLookAt(1:numIndicesToLookAt)
        !i2 = 148 ! debugging parameter
        !print*,'indicesToLookAt(i2)',indicesToLookAt(i2)
        !!print*,'vecAs(i2,1:3)',vecAs(i2,1:3)
        !!print*,'vecs(i2)',vecs(i2)
        !!print*,'edgeVectors(i2,1:2)',edgeVectors(i2,1:2)
        !!print*,'vecBs(i2,1:2)',vecBs(i2,1:2)
        !!print*,'edgeChecks(i2)',edgeChecks(i2)
        !print*,'F4(1205,1:3)',F4(1205,1:3)
        !!print*,Norm(Face2(8136),3)
        !print*,'areBlocking2(i2)',areBlocking2(i2)
        !print*,'areBlocking(i2)',areBlocking(i2)
        !print*,'is visible?',isVisible(k)
        !print*,'is visible (plus shadow)?',isVisiblePlusShadows(k)
        !!print*,'blocked by facet ID: ',indicesToLookat(findloc(areBlocking(1:numIndicesToLookAt),.true.))        
        !!print*,'blocked (shadow) by facet ID: ',indicesToLookat(findloc(areBlocking2(1:numIndicesToLookAt),.true.))        
        !!stop
        !end if ! debugging statement
        
    end do
    
    
end do

! bodge to set shadow facets to be not visible if there were no non-shadow facets visible from that aperture
numApertures = maxval(apertures) ! first, get the number of apertures, assuming they start at 1 and dont skip any numbers
!print*,'num apertures:',numApertures
allocate(visibleApertures(1:numApertures)) ! allocate array to hold whether an aperture was visible or not
! determine if any of the nonshadow facets were visible for each aperture 
visibleApertures = .false.
do i = 1, numApertures
    isApertureVisible = .false.
    do j = 1, size(Face1,1)
        if(isVisible(j) .and. apertures(j) .eq. i) visibleApertures(i) = .true.
    end do
    !print*,'i: ',i,'isApertureVisible',visibleApertures(i)   
end do
! if a shdaow facet was visible but part of an aperture that had no non-shadow visible facets, set shadow facets to not visible
do i = 1, size(Face1,1)
    if(visibleApertures(apertures(i)) .eqv. .false.) isVisiblePlusShadows(i) = .false.
end do


call CPU_TIME(finish)
! print'(A,f6.4,A)',"visible facets found in: ",finish-start," secs"
! print*,'========== end sr findVisibleFacets'

end subroutine

subroutine findVisibleFacetsInt(verts, Face1, Norm, Face2, midPoints, isVisible, apertures)

! subroutine findVisibleFacets finds the internal visible facets of the particle as viewed in the -z direction
! works for concave particles
! accuracy increases with mesh discretisation
! only difference between matlab version is that some facets were asigned to fuzzy bounding boxes differently
! may need to increase fac to improve consistency

real(8), dimension(:,:), allocatable, intent(in) :: verts, Norm, Midpoints
integer(8), dimension(:,:), allocatable, intent(in) :: Face1
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
logical, dimension(:), allocatable, intent(inout) :: isVisible
! logical, dimension(:), allocatable, intent(inout) :: isVisiblePlusShadows
integer(8), dimension(:), allocatable, intent(in) :: apertures

real(8), allocatable, dimension(:,:) :: boundingBoxV
integer(8), allocatable, dimension(:,:) :: boundingBoxF
integer(8) boundingBoxFSize
real(8), dimension(:,:), allocatable :: boundingBoxMidpoints ! unique vertices, face vertex IDs, face normals, face midpoints
real(8), dimension(:), allocatable :: boundingBoxFaceAreas
integer(8), dimension(:), allocatable :: boundingBoxNumFaceVert
integer, dimension(:), allocatable :: F3 ! bounding box IDs
integer, dimension(:,:), allocatable :: F4 ! fuzzy bounding box IDs
real(8), dimension(:), allocatable :: distanceToBB, distanceToFuzzyBB
integer i, j, k, l, m, BB, numUpFacingIndices, numIndicesToLookAt, noPoints
integer, dimension(:), allocatable :: upFacingIndices, indicesToLookAt
logical, dimension(:), allocatable :: areBlocking, areBlocking2
real(8), dimension(:,:), allocatable :: vecAs, edgeVectors, vecBs, edgeNormals
real(8), dimension(:), allocatable :: AdotNs, NdotKs, vecs, edgeChecks, nfs
real(8) start, finish
logical, dimension(:), allocatable :: visibleApertures
integer numApertures
logical isApertureVisible
real(8), dimension(:,:), allocatable :: flippedNorm ! variable to hold flipped normals (because this is internal instead of external)

! print*,'========== start sr findVisibleFacetsInt'
call CPU_TIME(start)

allocate(flippedNorm(1:size(Norm,1),1:3))
flippedNorm = -Norm ! flip normals because this is for internal interactions

!print*,'bounding box x dimension: ',bounding_box_x_dim
!print*,'bounding box y dimension: ',bounding_box_y_dim
!print*,'is isVisible allocated? ',allocated(isVisible)

! ### beam-aligned bounding boxes ###
!print*,'making beam-aligned bounding boxes'
! use the current crystal vertices to create some bounding boxes in x-y plane
call beam_aligned_bounding_boxes(verts, boundingBoxV, boundingBoxF)

allocate(boundingBoxNumFaceVert(1:size(boundingBoxF,1)))
boundingBoxNumFaceVert = 4

! compute bounding box midpoints
call midPointsAndAreas(boundingBoxF, boundingBoxV, boundingBoxMidpoints, boundingBoxFaceAreas, boundingBoxNumFaceVert)

! check midpoints
!do i = 1,size(boundingBoxMidpoints,1)
!    print'(A,f10.6,f10.6)','boundingBoxMidpoints',boundingBoxMidpoints(i,1),boundingBoxMidpoints(i,2)
!end do

allocate(F3(1:size(Face1,1))) ! array to hold index of bounding box that each face belongs to
allocate(F4(1:size(Face1,1),1:3)) ! array to hold index of fuzzy bounding box that each face belongs to
boundingBoxFSize = size(boundingBoxF,1) ! number of bounding box faces
allocate(distanceToBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box
allocate(distanceToFuzzyBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box

!print*,'Number of bounding boxes: ',boundingBoxFSize

! find which bounding box each vertex belongs to
do i = 1, size(Face1,1) ! for each face
!do i = 1, 1! for each face
    distanceToBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - Midpoints(i,1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - Midpoints(i,2))**2) ! distance to each bb
    !do j = 1,boundingBoxFSize
    !    print*,'dist to box: ',distanceToBB(j)
    !end do
    F3(i) = minloc(distanceToBB,1) ! record which bounding box midpoint this facet was closest to
    !print'(A,i6,A,f10.6)','i',i,'dist to closest bb: ',distanceToBB(F3(i))
    !print'(A,f10.6,f10.6)','midpoints of closest bb: ',boundingBoxMidpoints(F3(i),1),boundingBoxMidpoints(F3(i),2)
    do j = 1, 3 ! for each vertex
        distanceToFuzzyBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - verts(Face1(i,j),1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - verts(Face1(i,j),2))**2) ! distance to each bb
        F4(i,j) = minloc(distanceToFuzzyBB,1) ! record which bounding box midpoint this facet vertex was closest to
    end do 
end do

!do i = 1, size(Face1,1)
!    print*,'i',i,'F3(i): ',F3(i)
!end do


! next step is the loop for visible facets
isVisible = .true. ! start by assuming that all facets are unblocked and visible to the beam
! isVisiblePlusShadows = .true. ! start by assuming that all facets are unblocked and visible to the beam
do i = 1, size(Face1,1)
    if(flippedNorm(Face2(i),3) < 0.01) isVisible(i) = .false. ! down facing-facets are not visible (for now)
end do

! allocations
allocate(upFacingIndices(1:size(Face1,1)))
allocate(indicesToLookAt(1:size(Face1,1)))
allocate(areBlocking(1:size(Face1,1)))
allocate(areBlocking2(1:size(Face1,1)))
allocate(vecAs(1:size(Face1,1),1:3))
allocate(edgeVectors(1:size(Face1,1),1:2))
allocate(vecBs(1:size(Face1,1),1:2))
allocate(edgeNormals(1:size(Face1,1),1:2))
allocate(AdotNs(1:size(Face1,1)))
allocate(NdotKs(1:size(Face1,1)))
allocate(vecs(1:size(Face1,1)))
allocate(edgeChecks(1:size(Face1,1)))
allocate(nfs(1:size(Face1,1)))

do BB = 1, boundingBoxFSize ! for each bounding box
!do BB = 8, 8 ! for each bounding box
    !print*,'care, only using first BB atm'
    numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        !if(F3(i) .eq. BB .and. flippedNorm(Face2(i),3) > -0.01) then
        if(F3(i) .eq. BB) then
            numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
            upFacingIndices(numUpFacingIndices) = i
        end if        
    end do
    numIndicesToLookAt = 0 ! make a counter to count the number of faces within fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        if(any(F4(i,1:3) .eq. BB,1) .and. flippedNorm(Face2(i),3) > -0.01) then
            numIndicesToLookAt = numIndicesToLookAt + 1 ! count number of faces within non-fuzzy bounding box
            indicesToLookAt(numIndicesToLookAt) = i ! record ID of visible facet within fuzzy bounding box
        end if            
    end do
    !print*,'number of visible faces within this BB', numUpFacingIndices
    !print*,'visible face IDs within this BB',upFacingIndices(1:numUpFacingIndices)
    !print*,'number of visible faces within this fuzzy BB', numIndicesToLookAt
    !print*,'visible face IDs within this fuzzy BB',indicesToLookAt(1:numIndicesToLookAt)
    
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),1)
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),2)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),1)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),2)
    !print*,'size(Midpoints(i,1:3),1)',size(Midpoints(i,1:3),1)
    !print*,'size(Midpoints(i,1:3),2)',size(Midpoints(i,1:3),2)
    !stop
    do i = 1, numUpFacingIndices ! for each visible, upfacing face within this bounding box
        k = upFacingIndices(i) ! get face ID
        areBlocking(1:numIndicesToLookAt) = .true. ! logical array containing whether each face blocks face i
        vecAs(1:numIndicesToLookAt,1) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),1) - Midpoints(k,1) ! vectors from midpoint of face i to midpoints of faces
        vecAs(1:numIndicesToLookAt,2) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),2) - Midpoints(k,2)
        vecAs(1:numIndicesToLookAt,3) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),3) - Midpoints(k,3)
        AdotNs(1:numIndicesToLookAt) =  vecAs(1:numIndicesToLookAt,1)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),1) + & ! dot products of face normals with vector A (above)
                                        vecAs(1:numIndicesToLookAt,2)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),2) + &
                                        vecAs(1:numIndicesToLookAt,3)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3)
        NdotKs(1:numIndicesToLookAt) = -flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3) ! dot product of face normals with incident beam direction (z-axis)
        vecs(1:numIndicesToLookAt) = -AdotNs(1:numIndicesToLookAt) / NdotKs(1:numIndicesToLookAt)
        do j = 1, numIndicesToLookAt 
            if(vecs(j) .lt. 0) areBlocking(j) = .false. ! if faces are behind, they are not blocking
            if (indicesToLookAt(j) .eq. k) areBlocking(j) = .false. ! cant block itself
        end do
        
        noPoints = 3 ! number of vertices in blocking face
        do m = 1, noPoints ! for each vertex in blocking face
            ! compute vector from this vertex to the next vertex of the blocking facet
            if(m .eq. noPoints) then
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            else
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            end if
            vecBs(1:numIndicesToLookAt,1) = Midpoints(k,1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1) ! get the vector from the vertex of the faces to the midpoint of face i
            vecBs(1:numIndicesToLookAt,2) = Midpoints(k,2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            edgeNormals(1:numIndicesToLookAt,1) = edgeVectors(1:numIndicesToLookAt,2) ! cross product of edge vector with reverse beam direction
            edgeNormals(1:numIndicesToLookAt,2) = -edgeVectors(1:numIndicesToLookAt,1) ! flipped signs here relative to external subroutine
            nfs(1:numIndicesToLookAt) = sqrt(edgeNormals(1:numIndicesToLookAt,1)**2 + edgeNormals(1:numIndicesToLookAt,2)**2) ! normalisation factor
            edgeNormals(1:numIndicesToLookAt,1) = edgeNormals(1:numIndicesToLookAt,1) / nfs(1:numIndicesToLookAt) ! normalise
            edgeNormals(1:numIndicesToLookAt,2) = edgeNormals(1:numIndicesToLookAt,2) / nfs(1:numIndicesToLookAt)
            edgeChecks(1:numIndicesToLookAt) = -vecBs(1:numIndicesToLookAt,1)*edgeNormals(1:numIndicesToLookAt,1) - vecBs(1:numIndicesToLookAt,2)*edgeNormals(1:numIndicesToLookAt,2)
            do l = 1,numIndicesToLookAt
                if(edgeChecks(l) .gt. 0) areBlocking(l) = .false. ! if its not "inside" the edge vector, face i has failed the check and is therefore not blocked by the face
            end do
            !print*,'areBlocking(1:numIndicesToLookAt)',areBlocking(1:numIndicesToLookAt)
            
        end do
        ! shadow region stuff
        areBlocking2(1:numIndicesToLookAt) = areBlocking(1:numIndicesToLookAt) ! shadow region facets will have the same blockers as non-shadow region, but cannot be in the same aperture as the facet
        do m = 1, numIndicesToLookAt ! looping through all blockers
            if(apertures(indicesToLookat(m)) .eq. apertures(k)) areBlocking2(m) = .false. ! set as non-blocking if part of same aperture
        end do
        
        if(any(areBlocking(1:numIndicesToLookAt))) isVisible(k) = .false. ! visibility of up-facing facets
        ! if(any(areBlocking2(1:numIndicesToLookAt))) isVisiblePlusShadows(k) = .false. ! visibility of all facets of illuminated apertures, including shadow region
        
        ! if(k.eq.9215) then ! debugging statement
        ! print*,'looking for facets which block facet: ',k
        ! print*,'this facet is in BB: ',BB
        ! print*,'this BB has midpoints:',boundingBoxMidpoints(BB,1),boundingBoxMidpoints(BB,2)
        ! print*,'indicesToLookAt(1:numIndicesToLookAt)',indicesToLookAt(1:numIndicesToLookAt)
        ! i2 = 1 ! debugging parameter (look through indices to look at to pick)
        ! print*,'indicesToLookAt(i2)',indicesToLookAt(i2)
        ! !print*,'vecAs(i2,1:3)',vecAs(i2,1:3)
        ! !print*,'vecs(i2)',vecs(i2)
        ! !print*,'edgeVectors(i2,1:2)',edgeVectors(i2,1:2)
        ! !print*,'vecBs(i2,1:2)',vecBs(i2,1:2)
        ! !print*,'edgeChecks(i2)',edgeChecks(i2)
        ! ! print*,'F4(1205,1:3)',F4(1205,1:3)
        ! print*,'Norm(Face2(3449),1:3)',Norm(Face2(3449),1:3)f
        ! !print*,flippedNorm(Face2(8136),3)
        ! print*,'areBlocking2(i2)',areBlocking2(i2)
        ! print*,'areBlocking(i2)',areBlocking(i2)
        ! print*,'is visible?',isVisible(k)
        ! print*,'is visible (plus shadow)?',isVisiblePlusShadows(k)
        ! !print*,'blocked by facet ID: ',indicesToLookat(findloc(areBlocking(1:numIndicesToLookAt),.true.))        
        ! !print*,'blocked (shadow) by facet ID: ',indicesToLookat(findloc(areBlocking2(1:numIndicesToLookAt),.true.))        
        ! !stop
        ! end if ! debugging statement
        
    end do
    
    
end do

! bodge to set shadow facets to be not visible if there were no non-shadow facets visible from that aperture
numApertures = maxval(apertures) ! first, get the number of apertures, assuming they start at 1 and dont skip any numbers
!print*,'num apertures:',numApertures
allocate(visibleApertures(1:numApertures)) ! allocate array to hold whether an aperture was visible or not
! determine if any of the nonshadow facets were visible for each aperture 
visibleApertures = .false.
do i = 1, numApertures
    isApertureVisible = .false.
    do j = 1, size(Face1,1)
        if(isVisible(j) .and. apertures(j) .eq. i) visibleApertures(i) = .true.
    end do
    !print*,'i: ',i,'isApertureVisible',visibleApertures(i)   
end do
! if a shdaow facet was visible but part of an aperture that had no non-shadow visible facets, set shadow facets to not visible
! do i = 1, size(Face1,1)
!     if(visibleApertures(apertures(i)) .eq. .false.) isVisiblePlusShadows(i) = .false.
! end do


call CPU_TIME(finish)
print'(A,f6.4,A)',"visible facets found in: ",finish-start," secs"
! print*,'========== end sr findVisibleFacetsInt'

end subroutine

subroutine findVisibleFacetsInt_ps(verts, Face1, Norm, Face2, midPoints, isVisiblePlusShadows, apertures, rotatedapertureNormals)

! subroutine findVisibleFacets finds the internal visible facets of the particle as viewed in the -z direction
! works for concave particles
! accuracy increases with mesh discretisation
! only difference between matlab version is that some facets were asigned to fuzzy bounding boxes differently
! may need to increase fac to improve consistency

real(8), dimension(:,:), allocatable, intent(in) :: verts, Norm, Midpoints
integer, dimension(:,:), allocatable, intent(in) :: Face1
integer, dimension(:), allocatable, intent(in) :: Face2 ! face normal ID of each face
! logical, dimension(:), allocatable, intent(inout) :: isVisible
logical, dimension(:), allocatable, intent(inout) :: isVisiblePlusShadows
integer(8), dimension(:), allocatable, intent(in) :: apertures
real(8), dimension(:,:), allocatable, intent(in) :: rotatedapertureNormals

real(8), allocatable, dimension(:,:) :: boundingBoxV
integer(8), allocatable, dimension(:,:) :: boundingBoxF
integer boundingBoxFSize
real(8), dimension(:,:), allocatable :: boundingBoxMidpoints ! unique vertices, face vertex IDs, face normals, face midpoints
real(8), dimension(:), allocatable :: boundingBoxFaceAreas
integer(8), dimension(:), allocatable :: boundingBoxNumFaceVert
integer, dimension(:), allocatable :: F3 ! bounding box IDs
integer, dimension(:,:), allocatable :: F4 ! fuzzy bounding box IDs
real(8), dimension(:), allocatable :: distanceToBB, distanceToFuzzyBB
integer i, j, k, l, m, BB, numUpFacingIndices, numIndicesToLookAt, noPoints
integer, dimension(:), allocatable :: upFacingIndices, indicesToLookAt
logical, dimension(:), allocatable :: areBlocking, areBlocking2
real(8), dimension(:,:), allocatable :: vecAs, edgeVectors, vecBs, edgeNormals
real(8), dimension(:), allocatable :: AdotNs, NdotKs, vecs, edgeChecks, nfs
real(8) start, finish
real(8), dimension(:,:), allocatable :: flippedNorm ! variable to hold flipped normals (because this is internal instead of external)

! print*,'========== start sr findVisibleFacetsInt_ps'
call CPU_TIME(start)

allocate(flippedNorm(1:size(Norm,1),1:3))
flippedNorm = -Norm ! flip normals because this is for internal interactions

!print*,'bounding box x dimension: ',bounding_box_x_dim
!print*,'bounding box y dimension: ',bounding_box_y_dim
!print*,'is isVisible allocated? ',allocated(isVisible)

! ### beam-aligned bounding boxes ###
!print*,'making beam-aligned bounding boxes'
! use the current crystal vertices to create some bounding boxes in x-y plane
call beam_aligned_bounding_boxes(verts, boundingBoxV, boundingBoxF)

allocate(boundingBoxNumFaceVert(1:size(boundingBoxF,1)))
boundingBoxNumFaceVert = 4

! compute bounding box midpoints
call midPointsAndAreas(boundingBoxF, boundingBoxV, boundingBoxMidpoints, boundingBoxFaceAreas, boundingBoxNumFaceVert)

! check midpoints
!do i = 1,size(boundingBoxMidpoints,1)
!    print'(A,f10.6,f10.6)','boundingBoxMidpoints',boundingBoxMidpoints(i,1),boundingBoxMidpoints(i,2)
!end do

allocate(F3(1:size(Face1,1))) ! array to hold index of bounding box that each face belongs to
allocate(F4(1:size(Face1,1),1:3)) ! array to hold index of fuzzy bounding box that each face belongs to
boundingBoxFSize = size(boundingBoxF,1) ! number of bounding box faces
allocate(distanceToBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box
allocate(distanceToFuzzyBB(1:boundingBoxFSize)) ! array to hold the distance of a given vertex to each bounding box

!print*,'Number of bounding boxes: ',boundingBoxFSize

! find which bounding box each vertex belongs to
do i = 1, size(Face1,1) ! for each face
!do i = 1, 1! for each face
    distanceToBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - Midpoints(i,1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - Midpoints(i,2))**2) ! distance to each bb
    !do j = 1,boundingBoxFSize
    !    print*,'dist to box: ',distanceToBB(j)
    !end do
    F3(i) = minloc(distanceToBB,1) ! record which bounding box midpoint this facet was closest to
    !print'(A,i6,A,f10.6)','i',i,'dist to closest bb: ',distanceToBB(F3(i))
    !print'(A,f10.6,f10.6)','midpoints of closest bb: ',boundingBoxMidpoints(F3(i),1),boundingBoxMidpoints(F3(i),2)
    do j = 1, 3 ! for each vertex
        distanceToFuzzyBB(1:boundingBoxFSize) = sqrt((boundingBoxMidpoints(1:boundingBoxFSize,1) - verts(Face1(i,j),1))**2 + (boundingBoxMidpoints(1:boundingBoxFSize,2) - verts(Face1(i,j),2))**2) ! distance to each bb
        F4(i,j) = minloc(distanceToFuzzyBB,1) ! record which bounding box midpoint this facet vertex was closest to
    end do 
end do

!do i = 1, size(Face1,1)
!    print*,'i',i,'F3(i): ',F3(i)
!end do


! next step is the loop for visible facets
! isVisible = .true. ! start by assuming that all facets are unblocked and visible to the beam
isVisiblePlusShadows = .true. ! start by assuming that all facets are unblocked and visible to the beam
do i = 1, size(Face1,1)
    if(rotatedapertureNormals(apertures(i),3) > -0.01) isVisiblePlusShadows(i) = .false. ! set facets that belong to that aperture to be not visible, including shadow facets
end do

! allocations
allocate(upFacingIndices(1:size(Face1,1)))
allocate(indicesToLookAt(1:size(Face1,1)))
allocate(areBlocking(1:size(Face1,1)))
allocate(areBlocking2(1:size(Face1,1)))
allocate(vecAs(1:size(Face1,1),1:3))
allocate(edgeVectors(1:size(Face1,1),1:2))
allocate(vecBs(1:size(Face1,1),1:2))
allocate(edgeNormals(1:size(Face1,1),1:2))
allocate(AdotNs(1:size(Face1,1)))
allocate(NdotKs(1:size(Face1,1)))
allocate(vecs(1:size(Face1,1)))
allocate(edgeChecks(1:size(Face1,1)))
allocate(nfs(1:size(Face1,1)))

do BB = 1, boundingBoxFSize ! for each bounding box
!do BB = 8, 8 ! for each bounding box
    !print*,'care, only using first BB atm'
    numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        !if(F3(i) .eq. BB .and. flippedNorm(Face2(i),3) > -0.01) then
        if(F3(i) .eq. BB) then
            numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
            upFacingIndices(numUpFacingIndices) = i
        end if        
    end do
    numIndicesToLookAt = 0 ! make a counter to count the number of faces within fuzzy bounding box
    do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
        if(any(F4(i,1:3) .eq. BB,1) .and. rotatedapertureNormals(apertures(i),3) < -0.01) then
            numIndicesToLookAt = numIndicesToLookAt + 1 ! count number of faces within non-fuzzy bounding box
            indicesToLookAt(numIndicesToLookAt) = i ! record ID of visible facet within fuzzy bounding box
        end if            
    end do
    !print*,'number of visible faces within this BB', numUpFacingIndices
    !print*,'visible face IDs within this BB',upFacingIndices(1:numUpFacingIndices)
    !print*,'number of visible faces within this fuzzy BB', numIndicesToLookAt
    !print*,'visible face IDs within this fuzzy BB',indicesToLookAt(1:numIndicesToLookAt)
    
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),1)
    !print*,'indicesToLookAt(1:numIndicesToLookAt)',size(Midpoints(indicesToLookAt(1:numIndicesToLookAt),1:3),2)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),1)
    !print*,'size(vecAs(1:numIndicesToLookAt,1:3),1)',size(vecAs(1:numIndicesToLookAt,1:3),2)
    !print*,'size(Midpoints(i,1:3),1)',size(Midpoints(i,1:3),1)
    !print*,'size(Midpoints(i,1:3),2)',size(Midpoints(i,1:3),2)
    !stop
    do i = 1, numUpFacingIndices ! for each visible, upfacing face within this bounding box
        k = upFacingIndices(i) ! get face ID
        areBlocking(1:numIndicesToLookAt) = .true. ! logical array containing whether each face blocks face i
        vecAs(1:numIndicesToLookAt,1) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),1) - Midpoints(k,1) ! vectors from midpoint of face i to midpoints of faces
        vecAs(1:numIndicesToLookAt,2) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),2) - Midpoints(k,2)
        vecAs(1:numIndicesToLookAt,3) = Midpoints(indicesToLookAt(1:numIndicesToLookAt),3) - Midpoints(k,3)
        AdotNs(1:numIndicesToLookAt) =  vecAs(1:numIndicesToLookAt,1)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),1) + & ! dot products of face normals with vector A (above)
                                        vecAs(1:numIndicesToLookAt,2)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),2) + &
                                        vecAs(1:numIndicesToLookAt,3)*flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3)
        NdotKs(1:numIndicesToLookAt) = -flippedNorm(Face2(indicesToLookAt(1:numIndicesToLookAt)),3) ! dot product of face normals with incident beam direction (z-axis)
        vecs(1:numIndicesToLookAt) = -AdotNs(1:numIndicesToLookAt) / NdotKs(1:numIndicesToLookAt)
        do j = 1, numIndicesToLookAt 
            if(vecs(j) .lt. 0) areBlocking(j) = .false. ! if faces are behind, they are not blocking
            if (indicesToLookAt(j) .eq. k) areBlocking(j) = .false. ! cant block itself
        end do
        
        noPoints = 3 ! number of vertices in blocking face
        do m = 1, noPoints ! for each vertex in blocking face
            ! compute vector from this vertex to the next vertex of the blocking facet
            if(m .eq. noPoints) then
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            else
                edgeVectors(1:numIndicesToLookAt,1) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1)
                edgeVectors(1:numIndicesToLookAt,2) = verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m+1),2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            end if
            vecBs(1:numIndicesToLookAt,1) = Midpoints(k,1) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),1) ! get the vector from the vertex of the faces to the midpoint of face i
            vecBs(1:numIndicesToLookAt,2) = Midpoints(k,2) - verts(Face1(indicesToLookAt(1:numIndicesToLookAt),m),2)
            edgeNormals(1:numIndicesToLookAt,1) = edgeVectors(1:numIndicesToLookAt,2) ! cross product of edge vector with reverse beam direction
            edgeNormals(1:numIndicesToLookAt,2) = -edgeVectors(1:numIndicesToLookAt,1) ! flipped signs here relative to external subroutine
            nfs(1:numIndicesToLookAt) = sqrt(edgeNormals(1:numIndicesToLookAt,1)**2 + edgeNormals(1:numIndicesToLookAt,2)**2) ! normalisation factor
            edgeNormals(1:numIndicesToLookAt,1) = edgeNormals(1:numIndicesToLookAt,1) / nfs(1:numIndicesToLookAt) ! normalise
            edgeNormals(1:numIndicesToLookAt,2) = edgeNormals(1:numIndicesToLookAt,2) / nfs(1:numIndicesToLookAt)
            edgeChecks(1:numIndicesToLookAt) = -vecBs(1:numIndicesToLookAt,1)*edgeNormals(1:numIndicesToLookAt,1) - vecBs(1:numIndicesToLookAt,2)*edgeNormals(1:numIndicesToLookAt,2)
            do l = 1,numIndicesToLookAt
                if(edgeChecks(l) .gt. 0) areBlocking(l) = .false. ! if its not "inside" the edge vector, face i has failed the check and is therefore not blocked by the face
            end do
            !print*,'areBlocking(1:numIndicesToLookAt)',areBlocking(1:numIndicesToLookAt)
            
        end do
        ! shadow region stuff
        areBlocking2(1:numIndicesToLookAt) = areBlocking(1:numIndicesToLookAt) ! shadow region facets will have the same blockers as non-shadow region, but cannot be in the same aperture as the facet
        do m = 1, numIndicesToLookAt ! looping through all blockers
            if(apertures(indicesToLookat(m)) .eq. apertures(k)) areBlocking2(m) = .false. ! set as non-blocking if part of same aperture
        end do
        
        ! if(any(areBlocking(1:numIndicesToLookAt))) isVisible(k) = .false. ! visibility of up-facing facets
        if(any(areBlocking2(1:numIndicesToLookAt))) isVisiblePlusShadows(k) = .false. ! visibility of all facets of illuminated apertures, including shadow region
        
        ! if(k.eq.9215) then ! debugging statement
        ! print*,'looking for facets which block facet: ',k
        ! print*,'this facet is in BB: ',BB
        ! print*,'this BB has midpoints:',boundingBoxMidpoints(BB,1),boundingBoxMidpoints(BB,2)
        ! print*,'indicesToLookAt(1:numIndicesToLookAt)',indicesToLookAt(1:numIndicesToLookAt)
        ! i2 = 1 ! debugging parameter (look through indices to look at to pick)
        ! print*,'indicesToLookAt(i2)',indicesToLookAt(i2)
        ! !print*,'vecAs(i2,1:3)',vecAs(i2,1:3)
        ! !print*,'vecs(i2)',vecs(i2)
        ! !print*,'edgeVectors(i2,1:2)',edgeVectors(i2,1:2)
        ! !print*,'vecBs(i2,1:2)',vecBs(i2,1:2)
        ! !print*,'edgeChecks(i2)',edgeChecks(i2)
        ! ! print*,'F4(1205,1:3)',F4(1205,1:3)
        ! print*,'Norm(Face2(3449),1:3)',Norm(Face2(3449),1:3)f
        ! !print*,flippedNorm(Face2(8136),3)
        ! print*,'areBlocking2(i2)',areBlocking2(i2)
        ! print*,'areBlocking(i2)',areBlocking(i2)
        ! print*,'is visible?',isVisible(k)
        ! print*,'is visible (plus shadow)?',isVisiblePlusShadows(k)
        ! !print*,'blocked by facet ID: ',indicesToLookat(findloc(areBlocking(1:numIndicesToLookAt),.true.))        
        ! !print*,'blocked (shadow) by facet ID: ',indicesToLookat(findloc(areBlocking2(1:numIndicesToLookAt),.true.))        
        ! !stop
        ! end if ! debugging statement
        
    end do
    
    
end do

! ! bodge to set shadow facets to be not visible if there were no non-shadow facets visible from that aperture
! numApertures = maxval(apertures) ! first, get the number of apertures, assuming they start at 1 and dont skip any numbers
! !print*,'num apertures:',numApertures
! allocate(visibleApertures(1:numApertures)) ! allocate array to hold whether an aperture was visible or not
! ! determine if any of the nonshadow facets were visible for each aperture 
! visibleApertures = .false.
! do i = 1, numApertures
!     isApertureVisible = .false.
!     do j = 1, size(Face1,1)
!         if(isVisible(j) .and. apertures(j) .eq. i) visibleApertures(i) = .true.
!     end do
!     !print*,'i: ',i,'isApertureVisible',visibleApertures(i)   
! end do
! ! ! if a shdaow facet was visible but part of an aperture that had no non-shadow visible facets, set shadow facets to not visible
! do i = 1, size(Face1,1)
!     if(visibleApertures(apertures(i)) .eq. .false.) isVisiblePlusShadows(i) = .false.
! end do


call CPU_TIME(finish)
print'(A,f6.4,A)',"visible facets found in: ",finish-start," secs"
! print*,'========== end sr findVisibleFacetsInt_ps'

end subroutine

subroutine findWithinBeamInt(Face1, Midpoints, isVisible, beamV, beamF1, beamN, beamF2, beamMidpoints, isWithinBeam, distances, beamIDs)

integer, dimension(:,:), allocatable, intent(in) :: Face1, beamF1
real(8), dimension(:,:), allocatable, intent(in) :: Midpoints, beamMidpoints
logical, dimension(:), allocatable, intent(in) :: isVisible
real(8), dimension(:,:), allocatable, intent(in) :: beamV, beamN
integer, dimension(:), allocatable, intent(in) :: beamF2
logical, dimension(:), allocatable, intent(inout) :: isWithinBeam
real(8), dimension(:), allocatable, intent(inout) :: distances
integer, dimension(:), allocatable, intent(inout) :: beamIDs

real(8) start, finish
integer, dimension(:), allocatable :: upFacingIndices, indicesToLookAt
integer i, j, k, l, m, numUpFacingIndices, numIndicesToLookAt, noPoints
logical, dimension(:), allocatable :: areBlocking
real(8), dimension(:,:), allocatable :: vecAs, vecBs, edgeNormals, edgeVectors
real(8), dimension(:), allocatable :: AdotNs, NdotKs, vecs, edgeChecks, nfs
real(8), dimension(:,:), allocatable :: flippedNorm ! variable to hold flipped normals (because this is internal instead of external)
real(8) beamXmax0, beamXmin0, beamYmax0, beamYmin0
real(8) beamXmax, beamXmin, beamYmax, beamYmin

! print*,'========== start sr findWithinBeam'
call CPU_TIME(start)

allocate(flippedNorm(1:size(beamN,1),1:3))
flippedNorm = -beamN ! flip normals because this is for internal interactions

isWithinBeam = .false. ! start by assuming that all facets are not within beam
distances = 0 ! distances from beam to each visible face
beamIDs = 0 ! beam face IDs which blocked each face

! allocations
allocate(upFacingIndices(1:size(Face1,1)))
allocate(indicesToLookAt(1:size(beamF1,1)))
allocate(areBlocking(1:size(beamF1,1)))
allocate(vecAs(1:size(beamF1,1),1:3))
allocate(AdotNs(1:size(beamF1,1)))
allocate(NdotKs(1:size(Face1,1)))
allocate(vecs(1:size(Face1,1)))
allocate(edgeVectors(1:size(beamF1,1),1:2))
allocate(vecBs(1:size(beamF1,1),1:2))
allocate(edgeNormals(1:size(beamF1,1),1:2))
allocate(edgeChecks(1:size(beamF1,1)))
allocate(nfs(1:size(beamF1,1)))

numIndicesToLookAt = size(beamF1,1)
do i = 1, numIndicesToLookAt
    indicesToLookAt(i) = i
end do

! numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
! do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
!     if(isVisible(i)) then
!         numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
!         upFacingIndices(numUpFacingIndices) = i
!     end if        
! end do

! can speed up by looking at extreme cases
do i = 1, 3 ! bodge
    beamXmax0 = maxval(beamV(beamF1(1:size(beamF1,1),i),1))
    beamXmin0 = minval(beamV(beamF1(1:size(beamF1,1),i),1))
    beamYmax0 = maxval(beamV(beamF1(1:size(beamF1,1),i),2))
    beamYmin0 = minval(beamV(beamF1(1:size(beamF1,1),i),2))
    if(i .eq. 1) then
        beamXmax = beamXmax0
        beamXmin = beamXmin0
        beamYmax = beamYmax0
        beamYmin = beamYmin0
    else
        if(beamXmax0 .gt. beamXmax) beamXmax = beamXmax0
        if(beamXmin0 .lt. beamXmin) beamXmin = beamXmin0
        if(beamYmax0 .gt. beamYmax) beamYmax = beamYmax0
        if(beamYmin0 .lt. beamYmin) beamYmin = beamYmin0        
    end if
end do
numUpFacingIndices = 0 ! make a counter to count the number of faces within non-fuzzy bounding box
do i = 1, size(Face1,1) ! initial loop through to find how many facets are visible and within this bounding box
    if(isVisible(i) .and. Midpoints(i,1) .gt. beamXmin .and. Midpoints(i,1) .lt. beamXmax &
                    .and. Midpoints(i,2) .gt. beamYmin .and. Midpoints(i,2) .lt. beamYmax) then
        numUpFacingIndices = numUpFacingIndices + 1 ! count number of faces within non-fuzzy bounding box
        upFacingIndices(numUpFacingIndices) = i
    end if        
end do


do i = 1, numUpFacingIndices ! for each visible, upfacing face within this bounding box
    k = upFacingIndices(i) ! get face ID
    areBlocking(1:numIndicesToLookAt) = .true. ! start by assuming that this beam face blocks the facet
    vecAs(1:numIndicesToLookAt,1) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),1) - Midpoints(k,1) ! vectors from midpoint of face i to midpoints of faces
    vecAs(1:numIndicesToLookAt,2) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),2) - Midpoints(k,2)
    vecAs(1:numIndicesToLookAt,3) = beamMidpoints(indicesToLookAt(1:numIndicesToLookAt),3) - Midpoints(k,3)
    AdotNs(1:numIndicesToLookAt) =  vecAs(1:numIndicesToLookAt,1)*flippedNorm(beamF2(indicesToLookAt(1:numIndicesToLookAt)),1) + & ! dot products of face normals with vector A (above)
                                    vecAs(1:numIndicesToLookAt,2)*flippedNorm(beamF2(indicesToLookAt(1:numIndicesToLookAt)),2) + &
                                    vecAs(1:numIndicesToLookAt,3)*flippedNorm(beamF2(indicesToLookAt(1:numIndicesToLookAt)),3)
    NdotKs(1:numIndicesToLookAt) = -flippedNorm(beamF2(indicesToLookAt(1:numIndicesToLookAt)),3) ! dot product of face normals with incident beam direction (z-axis)
    vecs(1:numIndicesToLookAt) = -AdotNs(1:numIndicesToLookAt) / NdotKs(1:numIndicesToLookAt)
    do j = 1, numIndicesToLookAt 
        if(vecs(j) .lt. 0) areBlocking(j) = .false. ! if faces are behind, they are not blocking
    end do    
    noPoints = size(beamF1,2)
    do m = 1, noPoints ! for each vertex in blocking face
        ! compute vector from this vertex to the next vertex of the blocking facet
        if(m .eq. noPoints) then
            edgeNormals(1:numIndicesToLookAt,2) = -(beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),1),1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1))
            edgeNormals(1:numIndicesToLookAt,1) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),1),2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)
        else
            edgeNormals(1:numIndicesToLookAt,2) = -(beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m+1),1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1))
            edgeNormals(1:numIndicesToLookAt,1) = beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m+1),2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)
        end if
        nfs(1:numIndicesToLookAt) = sqrt(edgeNormals(1:numIndicesToLookAt,1)**2 + edgeNormals(1:numIndicesToLookAt,2)**2) ! normalisation factor
        edgeNormals(1:numIndicesToLookAt,1) = edgeNormals(1:numIndicesToLookAt,1) / nfs(1:numIndicesToLookAt) ! normalise
        edgeNormals(1:numIndicesToLookAt,2) = edgeNormals(1:numIndicesToLookAt,2) / nfs(1:numIndicesToLookAt)
        vecBs(1:numIndicesToLookAt,1) = Midpoints(k,1) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),1) ! get the vector from the vertex of the faces to the midpoint of face i
        vecBs(1:numIndicesToLookAt,2) = Midpoints(k,2) - beamV(beamF1(indicesToLookAt(1:numIndicesToLookAt),m),2)        
        edgeChecks(1:numIndicesToLookAt) = -(-vecBs(1:numIndicesToLookAt,1)*edgeNormals(1:numIndicesToLookAt,1) - vecBs(1:numIndicesToLookAt,2)*edgeNormals(1:numIndicesToLookAt,2))
        do l = 1,numIndicesToLookAt
            if(edgeChecks(l) .gt. 0) areBlocking(l) = .false. ! if its not "inside" the edge vector, face i has failed the check and is therefore not blocked by the face
        end do
        !print*,'areBlocking(1:numIndicesToLookAt)',areBlocking(1:numIndicesToLookAt)
            
    end do
    if(any(areBlocking(1:numIndicesToLookAt))) then
        isWithinBeam(k) = .true. ! visibility of up-facing facets
        distances(k) = vecAs(findloc(areBlocking(1:numIndicesToLookAt),.true.,1),3)
        beamIDs(k) = indicesToLookat(findloc(areBlocking(1:numIndicesToLookAt),.true.,1))

        ! print*,'blockerID,distances(k)',findloc(areBlocking(1:numIndicesToLookAt),.true.,1),distances(k)
        ! stop
    end if
        
    
end do

call CPU_TIME(finish)
print'(A,f12.8,A)',"facets within beam found in: ",finish-start," secs"
! print*,'========== end sr findWithinBeam'

!stop

end subroutine

subroutine init(face_ids, isVisible, isVisiblePlusShadows, isWithinBeam, distances, beamIDs, &
    isWithinBeam_ps, distances_ps, beamIDs_ps, isShadow, ampl_in, ampl_in_ps, la, waveno, &
    rperp, rpar, tperp, tpar, vk71, vk72, vk73, vk91, vk92, vk93, rot_ampl, new_in_ampl, &
    new_in_ampl_ps, trans_ampl_ps, trans_ampl, refl_ampl, refl_ampl_ps, &
    beam_outbeam_tree, beam_outbeam_tree_counter, interactionCounter)

! subroutine init initialises many variables for the beam tracing loop

logical, dimension(:), allocatable, intent(out) :: isVisible, isVisiblePlusShadows, isWithinBeam, isWithinBeam_ps, isShadow
integer(8), dimension(:,:), allocatable, intent(in) :: face_ids
real(8), dimension(:), allocatable, intent(out) :: distances
real(8), dimension(:), allocatable, intent(out) :: distances_ps
integer, dimension(:), allocatable, intent(out) :: beamIDs, beamIDs_ps
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
type(outbeamtype), dimension(:), allocatable, intent(out) :: beam_outbeam_tree ! outgoing beams from the beam tracing
integer, intent(out) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
integer, intent(out) :: interactionCounter ! counts the current number of interactions

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
allocate(beam_outbeam_tree(1:1000000)) ! set to 100000 as guess for max outbeams

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

ext_cross_section = 0 ! test
! print*,'========== end sr init'

end subroutine

subroutine get_beam_params( F1Mapping, &
                            isWithinBeam, &
                            ampl, &
                            vk71, &
                            vk72, &
                            vk73, &
                            aperturePropagationVectors, &
                            sufficientlyIlluminated, &
                            propagationVectors, &
                            i, &
                            illuminatedFaceIDs, &
                            refl_ampl_out11Int, &
                            refl_ampl_out12Int, &
                            refl_ampl_out21Int, &
                            refl_ampl_out22Int, &
                            vk71Int, vk72Int, vk73Int)

real(8), dimension(:,:,:), allocatable, intent(in) :: propagationVectors ! propagation vector of beam emitted from each aperture
integer, dimension(:,:), allocatable, intent(in) :: F1Mapping ! the face indices of each row of variables to go into main loop
integer, intent(in) :: i
complex(8), dimension(:,:), allocatable, intent(in) :: refl_ampl_out11Int ! the amplitude matrix that goes into the main loop
complex(8), dimension(:,:), allocatable, intent(in) :: refl_ampl_out12Int ! the amplitude matrix that goes into the main loop
complex(8), dimension(:,:), allocatable, intent(in) :: refl_ampl_out21Int ! the amplitude matrix that goes into the main loop
complex(8), dimension(:,:), allocatable, intent(in) :: refl_ampl_out22Int ! the amplitude matrix that goes into the main loop
real(8), dimension(:,:), allocatable, intent(in) :: vk71Int, vk72Int, vk73Int ! reflected e-perp vector

logical, dimension(:), allocatable, intent(inout) :: isWithinBeam
complex(8), dimension(:,:,:), allocatable, intent(inout) :: ampl
real(8), dimension(:), allocatable, intent(inout) :: vk71, vk72, vk73
real(8), dimension(:,:), allocatable, intent(inout) :: aperturePropagationVectors
logical, dimension(:), allocatable, intent(inout) :: sufficientlyIlluminated
integer, dimension(:), allocatable, intent(inout) :: illuminatedFaceIDs

integer j

    ! step 1: retrieve the parameters needed to propagate the next set of beams
    ! ampl, vk71, vk72, vk73, isWithinBeam, aperturePropagationVectors, sufficientlyIlluminated

    illuminatedFaceIDs(1:size(F1Mapping,1)) = F1Mapping(1:size(F1Mapping,1),i) ! extract the column
    isWithinBeam = .false. ! initialise
    ampl = 0 ! initialise
    vk71 = 0 ! initialise
    vk72 = 0 ! initialise
    vk73 = 0 ! initialise
    do j = 1, size(illuminatedFaceIDs,1) ! loop through face IDs
        if(illuminatedFaceIDs(j) .gt. 0) then ! if face was illuminated (had to change to 0 from nan vs. Matlab code)
            isWithinBeam(illuminatedFaceIDs(j)) = .true. ! set visible
            ampl(1,1,illuminatedFaceIDs(j)) = refl_ampl_out11Int(j,i)
            ampl(1,2,illuminatedFaceIDs(j)) = refl_ampl_out12Int(j,i)
            ampl(2,1,illuminatedFaceIDs(j)) = refl_ampl_out21Int(j,i)
            ampl(2,2,illuminatedFaceIDs(j)) = refl_ampl_out22Int(j,i)
            vk71(illuminatedFaceIDs(j)) = vk71Int(j,i)
            vk72(illuminatedFaceIDs(j)) = vk72Int(j,i)
            vk73(illuminatedFaceIDs(j)) = vk73Int(j,i)
        end if
    end do
    aperturePropagationVectors = 0 ! initialise
    sufficientlyIlluminated = .false. ! initialise
    do j = 1, size(propagationVectors,1) ! for each aperture
        aperturePropagationVectors(j,1) = propagationVectors(j,1,i) ! retrieve propagation vector
        aperturePropagationVectors(j,2) = propagationVectors(j,2,i)
        aperturePropagationVectors(j,3) = propagationVectors(j,3,i)
        ! if all propagation components arent set to 0, aperture was sufficiently illuminated
        if(sum(abs(aperturePropagationVectors(j,1:3))) .gt. 0.00001) sufficientlyIlluminated(j) = .true. ! bit of a bodge since cant do nans in fortran
    end do

end subroutine

subroutine get_interaction( InteractionInt2, &
                            FInt2, &
                            apertures, &
                            interactionCounter)

    integer, dimension(:,:), allocatable, intent(out) :: InteractionInt2
    integer, dimension(:,:), allocatable, intent(in) :: FInt2
    integer(8), dimension(:), allocatable, intent(in) :: apertures ! the aperture which each facet belongs to
    integer, intent(inout) :: interactionCounter ! counts the current number of interactions

    integer FTemp, counter
    integer, dimension(:), allocatable :: illuminated_apertures_temp
    integer, dimension(:), allocatable :: unique_illuminated_apertures
    integer j, k, m

    ! sort interactionOut array - bit messy but seems to do the job
    ! allocate an array to hold the interaction counter
    if(allocated(InteractionInt2)) deallocate(InteractionInt2)
    allocate(InteractionInt2(1:size(FInt2,1),1:size(FInt2,2)))
    InteractionInt2 = -1 ! init

    do j = 1, size(FInt2,2) ! for each illuminating aperture
        counter = 0 ! init counter
        do k = 1, size(FInt2,1) ! for each of the illuminated facets
            FTemp = FInt2(k,j) ! get the illuminated face ID
            if(FTemp .ne. 0) then ! if there was an illuminated face (0 is assumed as padding)
                counter = counter + 1 ! count the number of illluminated facets for this illuminating aperture
            end if
        end do
        ! print*,'column:',j,'had:',counter,' illuminated facets'
        if(counter .gt. 0) then ! skip if no illuminated facets
            if(allocated(illuminated_apertures_temp)) deallocate(illuminated_apertures_temp)
            allocate(illuminated_apertures_temp(1:counter))
            do k = 1, size(FInt2,1) ! for each of the illuminated facets
                FTemp = FInt2(k,j) ! get the illuminated face ID
                if(FTemp .ne. 0) then ! if there was an illuminated face (0 is assumed as padding)
                    illuminated_apertures_temp(k) = apertures(FTemp)
                end if
            end do
            call unique_int(illuminated_apertures_temp,unique_illuminated_apertures) ! get an array containing unique illuminated apertures
            ! print*,'unique_illuminated_apertures',unique_illuminated_apertures
            do k = 1, size(unique_illuminated_apertures,1) ! for each unique ill. aperture
                interactionCounter = interactionCounter + 1 ! update interaction counter
                do m = 1, size(FInt2,1) ! for each illuminated facet
                    FTemp = FInt2(m,j) ! get the illuminated face ID
                    if(FTemp .ne. 0) then
                        if(apertures(FTemp) .eq. unique_illuminated_apertures(k)) InteractionInt2(m,j) = interactionCounter
                    end if
                end do
            end do
        end if
    end do

end subroutine

subroutine add_to_outbeam_tree( FInt2, &
                                beam_outbeam_tree_counter, &
                                beam_outbeam_tree, &
                                trans_ampl_out11_2, &
                                trans_ampl_out12_2, &
                                trans_ampl_out21_2, &
                                trans_ampl_out22_2, &
                                vk71Int2, vk72Int2, vk73Int2, &
                                vk121Int2, vk122Int2, vk123Int2, &
                                vk91Int2, vk92Int2, vk93Int2, &
                                InteractionInt2)

integer, dimension(:,:), allocatable, intent(in) :: FInt2
integer, intent(inout) :: beam_outbeam_tree_counter ! counts the current number of beam outbeams
type(outbeamtype), dimension(:), allocatable, intent(inout) :: beam_outbeam_tree ! outgoing beams from the beam tracing
real(8), dimension(:,:), allocatable, intent(in) :: vk71Int2, vk72Int2, vk73Int2
real(8), dimension(:,:), allocatable, intent(in) :: vk121Int2 ,vk122Int2, vk123Int2
real(8), dimension(:,:), allocatable, intent(in) :: vk91Int2, vk92Int2, vk93Int2
complex(8), dimension(:,:), allocatable, intent(in) :: trans_ampl_out11_2, trans_ampl_out12_2, trans_ampl_out21_2, trans_ampl_out22_2
integer, dimension(:,:), allocatable, intent(in) :: InteractionInt2

integer j, k

do j = 1, size(FInt2,2) ! looping over illuminated apertures
    do k = 1, size(FInt2,1) ! looping over faces of the illuminated apertures
        if(FInt2(k,j) .eq. 0) then ! if we have reached some padding
            ! do nothing
        else ! else, add to the outbeam tree
            beam_outbeam_tree_counter = beam_outbeam_tree_counter + 1
            if(beam_outbeam_tree_counter .gt. size(beam_outbeam_tree,1)) then
                print*,'error: need more space in outbeam_tree. please increase in sr init'
            end if
                beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(1,1) = trans_ampl_out11_2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(1,2) = trans_ampl_out12_2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(2,1) = trans_ampl_out21_2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%ampl(2,2) = trans_ampl_out22_2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(1) = vk71Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(2) = vk72Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%vk7(3) = vk73Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(1) = vk121Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(2) = vk122Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_out(3) = vk123Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(1) = vk91Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(2) = vk92Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%prop_in(3) = vk93Int2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%FOut = FInt2(k,j)
                beam_outbeam_tree(beam_outbeam_tree_counter)%interactionOut = InteractionInt2(k,j)
        end if
    end do
end do

! print*,'finished adding to outbeam tree, ending at:',beam_outbeam_tree_counter,' outbeams'

end subroutine

subroutine get_beam(beam)

type(beam_type), intent(out) :: beam



end subroutine

end module beam_loop_mod
