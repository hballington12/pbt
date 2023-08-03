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

type beamtype ! type for beam trees (used for beam tracing)
! needs to contain:
! facet IDs, and at each facet ID, must contain:
	! perp-field direction
	! propagation direction
	! amplitude matrix
	! thats it???
	! maybe an interaction counter
	! complex(8) ampl(1:2,1:2) ! amplitude matrix
	! integer(8), dimension(:), allocatable :: FOut ! Illuminated face IDs
	! integer(8) interactionOut
	! real(8) prop
	! real(8) perp
end type beamtype

contains

end module types_mod