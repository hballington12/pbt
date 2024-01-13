! types_mod.f90
! contains data structure types shared across the project
    
module types_mod

implicit none

type outbeamtype ! type for outbeam trees (used for diffraction)
    complex(8) ampl(1:2,1:2) ! amplitude matrix
    real(8) vk7(1:3) ! perpendicular e-field direction
    real(8) verts(1:3,1:3) ! x, y, and z components of each vertex
    real(8) prop_out(1:3) ! outgoing propagation vector
    real(8) prop_in(1:3) ! incoming propagation vector
    integer(8) FOut ! face ID from which the beam was emitted
    integer interactionOut ! interaction counter
end type outbeamtype

type cc_hex_params_type ! type for holding parameters needed for C. Collier Gaussian Random hexagonal columns/plates
    real(8) l ! L from Muinonen & Saarinen - needs to be large compared to correlation length(s)
    real(8) hr ! hexagon radius hr
    integer nfhr ! number of subfacets along each hexagon edge - should be large enough to correctly plot the correlation length     nfhr
    real(8) pfl ! length of prism surfaces   pfl
    integer nfpl ! number of subfacets along prism facet length - should be even and large enough to correctly plot the correlation length nfpl
    integer pher ! number of rotations to perform at prism facet-basal facet edges (10% of #subfacets along prism edge) pher
    integer pper ! number of rotations to perform at prism facet-prism facet edges (10% of #subfacets along hexagon edge) pper
    integer nscales ! number of roughness scales sc required - supply correlation length and standard deviation for each below in pairs
    real(8), dimension(:), allocatable :: cls ! correlation lengths
    real(8), dimension(:), allocatable :: sds ! standard deviations
end type cc_hex_params_type

type output_parameters_type ! type for hold various output paramters
    real(8) abs ! absorption cross section
    real(8) scatt ! scattering cross section
    real(8) ext ! extinction cross section
    real(8) albedo ! single scattering albdeo
    real(8) asymmetry ! asymmetry parameter
    real(8) abs_eff ! absorption efficiency
    real(8) scatt_eff ! scattering efficiency
    real(8) ext_eff ! extinction efficiency
    real(8) geo_cross_sec ! illuminated geometric cross section
    real(8) back_scatt ! back-scattering cross section
    real(8) beam_energy_out ! beam energy into far-field
    real(8) ext_energy_out ! external diffraction energy into far-field
end type output_parameters_type 

type job_parameters_type
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
    integer  offs(1:2)
    real(8)  eulers(1:3)
    type(cc_hex_params_type) cc_hex_params ! parameters for C. Collier Gaussian Random hexagonal columns/plates
    real(8), dimension(:), allocatable :: theta_vals
    real(8), dimension(:), allocatable :: phi_vals
    logical suppress_2d ! whether or not to suppress 2d output of the mueller matrix
    logical tri ! enable auto triangulation
    real(8) tri_edge_length ! auto triangulation max edge length
    real(8) tri_roughness ! auto triangulation roughness magnitude
    real(8) time_limit ! job time limit in hours
    logical resume ! enable resume of cached data
    integer cache_id ! cached data to resume from
    logical scaling ! enable energy scaling for diffraction
    real(8) beta_lims(1:2) ! min and max beta values for orientation averaging
    real(8) gamma_lims(1:2) ! min and max gamma values for orientation averaging
    logical output_eulers ! enables an output of the euler angles to a file
    integer debug ! level of debugging output (0-3)
    logical timing ! code timing output
    real(8) threshold ! area threshold for new beams
end type job_parameters_type

type field_in_type
    ! field: contains the information about the electric field at a facet of the illuminating surface
    ! this is used as an allocatable entry in the beam type structure
    integer(8) fi ! the facet id
    complex(8) ampl(1:2,1:2) ! the amplitude matrix
    real(8) e_perp(1:3) ! electric field perpendicular vector (previously known as vk7)
end type field_in_type

type field_out_type
    ! field: contains the information about the electric field at a facet of the illuminated surface
    ! this is used as an allocatable entry in the beam type structure
    integer(8) fi ! the facet id
    complex(8) ampl_int(1:2,1:2) ! the internal amplitude matrix
    complex(8) ampl_ext(1:2,1:2) ! the external amplitude matrix
    real(8) e_perp(1:3) ! electric field perpendicular vector (previously known as vk7)
    integer(8) ap ! the aperture id
    logical is_tir ! whether or not this was total internal reflection
    real(8) prop_int(1:3) ! internal propagation direction (computed from aperture normal)
    real(8) prop_ext(1:3) ! external propagation direction (computed from facet normal)
    real(8) scatt_int ! the internal scattering cross section contribution of this facet
    real(8) scatt_ext ! the external scattering cross section contribution of this facet
    real(8) proj_area ! the area of this facet projected along the incoming beam propagation direction
end type field_out_type

type beam_type
    ! beam type: contains all information of a single beam
    ! a beam is a collection of facets with a single propagation direction
    ! the amplitude matrix is stored at the centroid of each facet
    ! generally, the information here is passed to the beam recursion subroutine to be propagated
    type(field_in_type), dimension(:), allocatable :: field_in ! information about the e-field at each facet in the illuminating surface (see above)
    type(field_out_type), dimension(:), allocatable :: field_out ! information about the e-field at each facet in the illuminated surface (see above)
    real(8) prop(1:3) ! the propagation direction of the beam
    integer(8) nf_in ! total number of facets that belong to this beam
    integer(8) nf_out ! total number of facets illuminated by this beam
    integer(8) ap ! the aperture from which this beam is propagating
    real(8) abs ! absorption cross section (energy absorbed)
    real(8) scatt_in ! input scattering cross section of the beam
    real(8) scatt_out ! output scattering cross section of the beam
    logical is_int ! whether or not the beam is propagating inside the particle
    integer(8) id ! a label for the position of a beam in a beam tree
    real(8) proj_area_in ! the total area of all facets in this beam when projected along the beam propagation direction
    real(8) proj_area_out ! the total area of all illuminated facets when projected along the beam propagation direction
end type beam_type

type facet_type
    ! facet type: information about each facet of the particle geometry
    integer(8), dimension(:), allocatable :: vi ! vertex index - points to a position in the vertex array
    integer(8) ni ! normal index - points to a position in the normal array
    real(8) mid(1:3) ! midpoint
    real(8) area ! area
    integer(8) nv ! number of vertices on this face
    integer(8) ap ! the aperture to which this face belongs
end type facet_type

type aperture_type
    ! aperture type: information about an aperture of the particle geometry
    real(8) n(1:3) ! normal
    real(8) mid(1:3) ! midpoint
    real(8) area
    integer(8) nf ! number of faces in this aperture
end type aperture_type

type geometry_type
    ! geometry type: contains information about a particle geometry
    ! a particle geometry is defined by its vertices and facets
    ! vertices in each facet should be ordered in an anti-clockwise fashion as viewed from outside the particle
    real(8), dimension(:,:), allocatable :: v ! vertices in the geometry, dimension N x 3
    real(8), dimension(:,:), allocatable :: n ! normals in the geometry, dimension N x 3
    type(facet_type), dimension(:), allocatable :: f ! data structure with information about each facet in the geometry
    integer(8) nv ! total number of unique vertices in the geometry
    integer(8) nf ! number of faces in the geometry
    integer(8) nn ! number of unique normals
    integer(8) na ! number of apertures
    type(aperture_type), dimension(:), allocatable :: ap ! data structure with information about each aperture in the geometry
end type geometry_type

! format specifiers

character(len=111), parameter :: fmt_mueller_2d = '(f12.4,f12.4,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8)'
character(len=105), parameter :: fmt_mueller_1d = '(f12.4,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8)'


contains

end module types_mod