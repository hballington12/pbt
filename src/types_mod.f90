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
	integer(8) interactionOut ! interaction counter
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
end type output_parameters_type 

contains

end module types_mod