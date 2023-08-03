! outputs_mod.f90
! contains subroutines for outputting
    
module outputs_mod

   use misc_submod

   implicit none
   
   
   contains
   
      subroutine finalise( ampl_far_beam11,     & ! amplitude matrix (1,1) due to beam diffraction
                           ampl_far_beam12,     & ! amplitude matrix (1,2) due to beam diffraction
                           ampl_far_beam21,     & ! amplitude matrix (2,1) due to beam diffraction
                           ampl_far_beam22,     & ! amplitude matrix (2,2) due to beam diffraction
                           ampl_far_ext_diff11, & ! amplitude matrix (1,1) due to external diffraction
                           ampl_far_ext_diff12, & ! amplitude matrix (1,2) due to external diffraction
                           ampl_far_ext_diff21, & ! amplitude matrix (2,1) due to external diffraction
                           ampl_far_ext_diff22, & ! amplitude matrix (2,2) due to external diffraction
                           theta_vals, &
                           phi_vals, &
                           energy_out_beam, &
                           energy_out_ext_diff, &
                           mueller, &
                           mueller_1d,&
                           la)   

   complex(8), dimension(:,:), allocatable, intent(in) :: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! beam
   complex(8), dimension(:,:), allocatable, intent(in)  :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! ext diff
   real(8), dimension(:), allocatable, intent(in) :: theta_vals, phi_vals
   real(8), intent(in) :: energy_out_beam
   real(8), intent(in) :: energy_out_ext_diff
   real(8), dimension(:,:,:), allocatable, intent(out) :: mueller ! mueller matrices
   real(8), dimension(:,:), allocatable, intent(out) :: mueller_1d ! phi-integrated mueller matrices
   real(8), intent(in) :: la ! wavelength

   complex(8), dimension(:,:), allocatable :: ampl_far11, ampl_far12, ampl_far21, ampl_far22 ! total
   real(8), dimension(:,:,:), allocatable :: mueller_beam, mueller_ext_diff ! mueller matrices
   real(8), dimension(:,:), allocatable :: mueller_beam_1d, mueller_ext_diff_1d ! phi-integrated mueller matrices
   integer i, j
   real(8) scatt, scatt_beam, scatt_ext_diff, asymmetry, asymmetry_beam, asymmetry_ext_diff
   real(8) waveno

   waveno = 2*pi/la

   ! allocate total amplitude matrix (beam - ext diff)
   allocate(ampl_far11(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
   allocate(ampl_far12(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
   allocate(ampl_far21(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
   allocate(ampl_far22(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))

   ! allocate mueller matrix
   allocate(mueller(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element
   allocate(mueller_beam(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element
   allocate(mueller_ext_diff(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element

   ! allocate 1d mueller matrix
   allocate(mueller_1d(1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element
   allocate(mueller_beam_1d(1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element
   allocate(mueller_ext_diff_1d(1:size(ampl_far_beam11,2),1:16)) ! 1:1 is for each element

   ! far field = beam diffraction - external diffraction
   ampl_far11 = ampl_far_beam11 - ampl_far_ext_diff11
   ampl_far12 = ampl_far_beam12 - ampl_far_ext_diff12
   ampl_far21 = ampl_far_beam21 - ampl_far_ext_diff21
   ampl_far22 = ampl_far_beam22 - ampl_far_ext_diff22

   print*,'making 2d mueller matrices...'

   ! add to mueller matrices
   do i = 1, size(ampl_far_beam11,2) ! loop over theta
      do j = 1, size(ampl_far_beam11,1) ! loop over phi
         ! Bohren & Huffman
         ! ##########################--S11--##########################
         ! s11 beam
         mueller_beam(j,i,1) =     real(0.5*(ampl_far_beam11(j,i)*conjg(ampl_far_beam11(j,i)) + &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam12(j,i)) + &
                                             ampl_far_beam21(j,i)*conjg(ampl_far_beam21(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam22(j,i))))    
         ! s11 ext diff
         mueller_ext_diff(j,i,1) = real(0.5*(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff11(j,i)) + &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff12(j,i)) + &
                                             ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff21(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff22(j,i))))   
         ! s11 total
         mueller(j,i,1) =          real(0.5*(ampl_far11(j,i)*conjg(ampl_far11(j,i)) + &
                                             ampl_far12(j,i)*conjg(ampl_far12(j,i)) + &
                                             ampl_far21(j,i)*conjg(ampl_far21(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far22(j,i))))   
         ! ##########################--S12--##########################                                                
         ! s12 beam
         mueller_beam(j,i,2) =     real(0.5*(ampl_far_beam11(j,i)*conjg(ampl_far_beam11(j,i)) - &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam12(j,i)) + &
                                             ampl_far_beam21(j,i)*conjg(ampl_far_beam21(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam22(j,i))))    
         ! s12 ext diff
         mueller_ext_diff(j,i,2) = real(0.5*(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff11(j,i)) - &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff12(j,i)) + &
                                             ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff21(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff22(j,i))))   
         ! s12 total
         mueller(j,i,2) =          real(0.5*(ampl_far11(j,i)*conjg(ampl_far11(j,i)) - &
                                             ampl_far12(j,i)*conjg(ampl_far12(j,i)) + &
                                             ampl_far21(j,i)*conjg(ampl_far21(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far22(j,i)))) 
         ! ##########################--S13--##########################
         ! s13 beam
         mueller_beam(j,i,3) =          real(ampl_far_beam11(j,i)*conjg(ampl_far_beam12(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s13 ext diff
         mueller_ext_diff(j,i,3) =      real(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff12(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s13 total
         mueller(j,i,3) =               real(ampl_far11(j,i)*conjg(ampl_far12(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far21(j,i)))                                                
         ! ##########################--S14--##########################
         ! s14 beam
         mueller_beam(j,i,4) =          imag(ampl_far_beam11(j,i)*conjg(ampl_far_beam12(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s14 ext diff
         mueller_ext_diff(j,i,4) =      imag(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff12(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s14 total
         mueller(j,i,4) =               imag(ampl_far11(j,i)*conjg(ampl_far12(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far21(j,i)))
         ! ##########################--S21--##########################
         ! s21 beam
         mueller_beam(j,i,5) =     real(0.5*(ampl_far_beam11(j,i)*conjg(ampl_far_beam11(j,i)) + &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam12(j,i)) - &
                                             ampl_far_beam21(j,i)*conjg(ampl_far_beam21(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam22(j,i))))    
         ! s21 ext diff
         mueller_ext_diff(j,i,5) = real(0.5*(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff11(j,i)) + &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff12(j,i)) - &
                                             ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff21(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff22(j,i))))   
         ! s21 total
         mueller(j,i,5) =          real(0.5*(ampl_far11(j,i)*conjg(ampl_far11(j,i)) + &
                                             ampl_far12(j,i)*conjg(ampl_far12(j,i)) - &
                                             ampl_far21(j,i)*conjg(ampl_far21(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far22(j,i))))
         ! ##########################--S22--##########################                                                             
         ! s22 beam
         mueller_beam(j,i,6) =     real(0.5*(ampl_far_beam11(j,i)*conjg(ampl_far_beam11(j,i)) - &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam12(j,i)) - &
                                             ampl_far_beam21(j,i)*conjg(ampl_far_beam21(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam22(j,i))))    
         ! s22 ext diff
         mueller_ext_diff(j,i,6) = real(0.5*(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff11(j,i)) - &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff12(j,i)) - &
                                             ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff21(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff22(j,i))))   
         ! s22 total
         mueller(j,i,6) =          real(0.5*(ampl_far11(j,i)*conjg(ampl_far11(j,i)) - &
                                             ampl_far12(j,i)*conjg(ampl_far12(j,i)) - &
                                             ampl_far21(j,i)*conjg(ampl_far21(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far22(j,i))))
         ! ##########################--S23--##########################   
         ! s23 beam
         mueller_beam(j,i,7) =          real(ampl_far_beam11(j,i)*conjg(ampl_far_beam12(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s23 ext diff
         mueller_ext_diff(j,i,7) =      real(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff12(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s23 total
         mueller(j,i,7) =               real(ampl_far11(j,i)*conjg(ampl_far12(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far21(j,i)))
         ! ##########################--S24--##########################   
         ! s24 beam
         mueller_beam(j,i,8) =          imag(ampl_far_beam11(j,i)*conjg(ampl_far_beam12(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s24 ext diff
         mueller_ext_diff(j,i,8) =      imag(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff12(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s24 total
         mueller(j,i,8) =               imag(ampl_far11(j,i)*conjg(ampl_far12(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far21(j,i))) 
         ! ##########################--S31--##########################   
         ! s31 beam
         mueller_beam(j,i,9) =          real(ampl_far_beam11(j,i)*conjg(ampl_far_beam21(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam12(j,i))) 
         ! s31 ext diff
         mueller_ext_diff(j,i,9) =      real(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff21(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff12(j,i)))   
         ! s31 total
         mueller(j,i,9) =               real(ampl_far11(j,i)*conjg(ampl_far21(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far12(j,i)))
         ! ##########################--S32--##########################  
         ! s32 beam
         mueller_beam(j,i,10) =         real(ampl_far_beam11(j,i)*conjg(ampl_far_beam21(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam12(j,i))) 
         ! s32 ext diff
         mueller_ext_diff(j,i,10) =     real(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff21(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff12(j,i)))   
         ! s32 total
         mueller(j,i,10) =              real(ampl_far11(j,i)*conjg(ampl_far21(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far12(j,i))) 
         ! ##########################--S33--##########################  
         ! s33 beam
         mueller_beam(j,i,11) =         real(ampl_far_beam11(j,i)*conjg(ampl_far_beam22(j,i)) + &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s33 ext diff
         mueller_ext_diff(j,i,11) =     real(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff22(j,i)) + &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s33 total
         mueller(j,i,11) =              real(ampl_far11(j,i)*conjg(ampl_far22(j,i)) + &
                                             ampl_far12(j,i)*conjg(ampl_far21(j,i)))  
         ! ##########################--S34--##########################   
         ! s34 beam
         mueller_beam(j,i,12) =         imag(ampl_far_beam11(j,i)*conjg(ampl_far_beam22(j,i)) + &
                                             ampl_far_beam21(j,i)*conjg(ampl_far_beam12(j,i))) 
         ! s34 ext diff
         mueller_ext_diff(j,i,12) =     imag(ampl_far_ext_diff11(j,i)*conjg(ampl_far_ext_diff22(j,i)) + &
                                             ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff12(j,i)))   
         ! s34 total
         mueller(j,i,12) =              imag(ampl_far11(j,i)*conjg(ampl_far22(j,i)) + &
                                             ampl_far21(j,i)*conjg(ampl_far12(j,i))) 
         ! ##########################--S41--########################## 
         ! s41 beam
         mueller_beam(j,i,13) =         imag(ampl_far_beam21(j,i)*conjg(ampl_far_beam11(j,i)) + &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam12(j,i))) 
         ! s41 ext diff
         mueller_ext_diff(j,i,13) =     imag(ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff11(j,i)) + &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff12(j,i)))   
         ! s41 total
         mueller(j,i,13) =              imag(ampl_far21(j,i)*conjg(ampl_far11(j,i)) + &
                                             ampl_far22(j,i)*conjg(ampl_far12(j,i)))  
         ! ##########################--S42--##########################  
         ! s42 beam
         mueller_beam(j,i,14) =         imag(ampl_far_beam21(j,i)*conjg(ampl_far_beam11(j,i)) - &
                                             ampl_far_beam22(j,i)*conjg(ampl_far_beam12(j,i))) 
         ! s42 ext diff
         mueller_ext_diff(j,i,14) =     imag(ampl_far_ext_diff21(j,i)*conjg(ampl_far_ext_diff11(j,i)) - &
                                             ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff12(j,i)))   
         ! s42 total
         mueller(j,i,14) =              imag(ampl_far21(j,i)*conjg(ampl_far11(j,i)) - &
                                             ampl_far22(j,i)*conjg(ampl_far12(j,i)))  
         ! ##########################--S43--##########################  
         ! s43 beam
         mueller_beam(j,i,15) =         imag(ampl_far_beam22(j,i)*conjg(ampl_far_beam11(j,i)) - &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s43 ext diff
         mueller_ext_diff(j,i,15) =     imag(ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff11(j,i)) - &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s43 total
         mueller(j,i,15) =              imag(ampl_far22(j,i)*conjg(ampl_far11(j,i)) - &
                                             ampl_far12(j,i)*conjg(ampl_far21(j,i)))  
         ! ##########################--S44--##########################   
         ! s44 beam
         mueller_beam(j,i,16) =         real(ampl_far_beam22(j,i)*conjg(ampl_far_beam11(j,i)) - &
                                             ampl_far_beam12(j,i)*conjg(ampl_far_beam21(j,i))) 
         ! s44 ext diff
         mueller_ext_diff(j,i,16) =     real(ampl_far_ext_diff22(j,i)*conjg(ampl_far_ext_diff11(j,i)) - &
                                             ampl_far_ext_diff12(j,i)*conjg(ampl_far_ext_diff21(j,i)))   
         ! s44 total
         mueller(j,i,16) =              real(ampl_far22(j,i)*conjg(ampl_far11(j,i)) - &
                                             ampl_far12(j,i)*conjg(ampl_far21(j,i)))     


      end do
   end do

   print*,'making 1d mueller matrices...'

   ! phi integrations...
   do i = 1, size(theta_vals,1) ! for each theta bin...
      
      ! p11
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,1),mueller_1d(i,1)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,1),mueller_beam_1d(i,1)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,1),mueller_ext_diff_1d(i,1)) ! integrate ext diff

      ! p12
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,2),mueller_1d(i,2)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,2),mueller_beam_1d(i,2)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,2),mueller_ext_diff_1d(i,2)) ! integrate ext diff      

      ! p13
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,3),mueller_1d(i,3)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,3),mueller_beam_1d(i,3)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,3),mueller_ext_diff_1d(i,3)) ! integrate ext diff          
   
      ! p14
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,4),mueller_1d(i,4)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,4),mueller_beam_1d(i,4)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,4),mueller_ext_diff_1d(i,4)) ! integrate ext diff   

      ! p21
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,5),mueller_1d(i,5)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,5),mueller_beam_1d(i,5)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,5),mueller_ext_diff_1d(i,5)) ! integrate ext diff  

      ! p22
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,6),mueller_1d(i,6)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,6),mueller_beam_1d(i,6)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,6),mueller_ext_diff_1d(i,6)) ! integrate ext diff  
      
      ! p23
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,7),mueller_1d(i,7)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,7),mueller_beam_1d(i,7)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,7),mueller_ext_diff_1d(i,7)) ! integrate ext diff  
      
      ! p24
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,8),mueller_1d(i,8)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,8),mueller_beam_1d(i,8)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,8),mueller_ext_diff_1d(i,8)) ! integrate ext diff    
      
      ! p31
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,9),mueller_1d(i,9)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,9),mueller_beam_1d(i,9)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,9),mueller_ext_diff_1d(i,9)) ! integrate ext diff  

      ! p32
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,10),mueller_1d(i,10)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,10),mueller_beam_1d(i,10)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,10),mueller_ext_diff_1d(i,10)) ! integrate ext diff  
      
      ! p33
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,11),mueller_1d(i,11)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,11),mueller_beam_1d(i,11)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,11),mueller_ext_diff_1d(i,11)) ! integrate ext diff  
      
      ! p34
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,12),mueller_1d(i,12)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,12),mueller_beam_1d(i,12)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,12),mueller_ext_diff_1d(i,12)) ! integrate ext diff   
      
      ! p41
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,13),mueller_1d(i,13)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,13),mueller_beam_1d(i,13)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,13),mueller_ext_diff_1d(i,13)) ! integrate ext diff  

      ! p42
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,14),mueller_1d(i,14)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,14),mueller_beam_1d(i,14)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,14),mueller_ext_diff_1d(i,14)) ! integrate ext diff  
      
      ! p43
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,15),mueller_1d(i,15)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,15),mueller_beam_1d(i,15)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,15),mueller_ext_diff_1d(i,15)) ! integrate ext diff  
      
      ! p44
      call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,16),mueller_1d(i,16)) ! integrate total
      call simpne(size(phi_vals,1),phi_vals,mueller_beam(1:size(mueller,1),i,16),mueller_beam_1d(i,16)) ! integrate beam
      call simpne(size(phi_vals,1),phi_vals,mueller_ext_diff(1:size(mueller,1),i,16),mueller_ext_diff_1d(i,16)) ! integrate ext diff   

   end do

   print*,'calculating integrated parameters...'

   ! theta integrations...
   call simpne(size(theta_vals,1),theta_vals,mueller_1d(1:size(theta_vals,1),1)*sin(theta_vals),scatt) ! p11*sin(theta)
   print'(A40,f16.8)','scattering cross section (total):',scatt

   call simpne(size(theta_vals,1),theta_vals,mueller_beam_1d(1:size(theta_vals,1),1)*sin(theta_vals),scatt_beam) ! p11*sin(theta)
   print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (beam):',scatt_beam," (",scatt_beam/energy_out_beam*100," %)"

   call simpne(size(theta_vals,1),theta_vals,mueller_ext_diff_1d(1:size(theta_vals,1),1)*sin(theta_vals),scatt_ext_diff) ! p11*sin(theta)
   print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (ext diff):',scatt_ext_diff," (",scatt_ext_diff/energy_out_ext_diff*100," %)"

   call simpne(size(theta_vals,1),theta_vals,mueller_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt,asymmetry) ! p11*sin(theta)
   print'(A40,f16.8)','asymmetry parameter (total):',asymmetry

   call simpne(size(theta_vals,1),theta_vals,mueller_beam_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt_beam,asymmetry_beam) ! p11*sin(theta)
   print'(A40,f16.8)','asymmetry parameter (beam):',asymmetry_beam

   call simpne(size(theta_vals,1),theta_vals,mueller_ext_diff_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt_ext_diff,asymmetry_ext_diff) ! p11*sin(theta)
   print'(A40,f16.8)','asymmetry parameter (ext diff):',asymmetry_ext_diff  

   print*,'scattering cross section via optical theorem:',abs(2*pi/waveno*imag(ampl_far11(1,1) + ampl_far22(1,1))) ! Jackson 10.137

   ! open(10,file="mueller_scatgrid")
   ! do i = 1, size(ampl_far_beam11,2)
   !    do j = 1, size(ampl_far_beam11,1)
   !       write(10,'(f12.4,f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
   !       theta_vals(i)*180/pi, phi_vals(j)*180/pi, &
   !       mueller_beam(j,i,1), mueller_beam(j,i,2), mueller_beam(j,i,3), mueller_beam(j,i,4), &
   !       mueller_beam(j,i,5), mueller_beam(j,i,6), mueller_beam(j,i,7), mueller_beam(j,i,8), &
   !       mueller_beam(j,i,9), mueller_beam(j,i,10), mueller_beam(j,i,11), mueller_beam(j,i,12), &
   !       mueller_beam(j,i,13), mueller_beam(j,i,14), mueller_beam(j,i,15), mueller_beam(j,i,16)
                                                                                             
   !    end do
   ! end do
   ! close(10)   

   ! open(10,file="mueller_scatgrid")
   ! do i = 1, size(ampl_far_beam11,2)
   !    do j = 1, size(ampl_far_beam11,1)
   !       write(10,'(f12.4,f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
   !       theta_vals(i)*180/pi, phi_vals(j)*180/pi, &
   !       mueller_ext_diff(j,i,1), mueller_ext_diff(j,i,2), mueller_ext_diff(j,i,3), mueller_ext_diff(j,i,4), &
   !       mueller_ext_diff(j,i,5), mueller_ext_diff(j,i,6), mueller_ext_diff(j,i,7), mueller_ext_diff(j,i,8), &
   !       mueller_ext_diff(j,i,9), mueller_ext_diff(j,i,10), mueller_ext_diff(j,i,11), mueller_ext_diff(j,i,12), &
   !       mueller_ext_diff(j,i,13), mueller_ext_diff(j,i,14), mueller_ext_diff(j,i,15), mueller_ext_diff(j,i,16)
                                                                                             
   !    end do
   ! end do
   ! close(10)    

   end subroutine

   subroutine writeup(mueller, mueller_1d, theta_vals, phi_vals)

   real(8), dimension(:,:,:), allocatable, intent(in) :: mueller ! mueller matrices
   real(8), dimension(:,:), allocatable, intent(in) :: mueller_1d ! phi-integrated mueller matrices
   real(8), dimension(:), allocatable, intent(in) :: theta_vals, phi_vals

   integer i, j

   print*,'writing mueller to file...'

   open(10,file="mueller_scatgrid")
   do i = 1, size(theta_vals,1)
      do j = 1, size(phi_vals,1)
         write(10,'(f12.4,f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
         theta_vals(i)*180/pi, phi_vals(j)*180/pi, &
         mueller(j,i,1), mueller(j,i,2), mueller(j,i,3), mueller(j,i,4), &
         mueller(j,i,5), mueller(j,i,6), mueller(j,i,7), mueller(j,i,8), &
         mueller(j,i,9), mueller(j,i,10), mueller(j,i,11), mueller(j,i,12), &
         mueller(j,i,13), mueller(j,i,14), mueller(j,i,15), mueller(j,i,16)                                                                             
      end do
   end do
   close(10)

   open(10,file="mueller_scatgrid_1d")
   do j = 1, size(theta_vals,1)
      write(10,'(f12.4,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') &
      theta_vals(j)*180/pi, &
      mueller_1d(j,1), mueller_1d(j,2), mueller_1d(j,3), mueller_1d(j,4), &
      mueller_1d(j,5), mueller_1d(j,6), mueller_1d(j,7), mueller_1d(j,8), &
      mueller_1d(j,9), mueller_1d(j,10), mueller_1d(j,11), mueller_1d(j,12), &
      mueller_1d(j,13), mueller_1d(j,14), mueller_1d(j,15), mueller_1d(j,16)   
   end do
   close(10)

   end subroutine

   subroutine summation(mueller, mueller_total, mueller_1d, mueller_1d_total)

      real(8), dimension(:,:,:), allocatable, intent(in) :: mueller ! mueller matrices
      real(8), dimension(:,:,:), allocatable, intent(inout) :: mueller_total ! mueller matrices
      real(8), dimension(:,:), allocatable, intent(in) :: mueller_1d ! phi-integrated mueller matrices
      real(8), dimension(:,:), allocatable, intent(inout) :: mueller_1d_total ! phi-integrated mueller matrices

      ! if its the first call to summation, allocate the total mueller 1d and 2d arrays
      if(.not. allocated(mueller_total)) then
         allocate(mueller_total(1:size(mueller,1),1:size(mueller,2),1:size(mueller,3)))
         mueller_total = 0 ! init
      end if
      if(.not. allocated(mueller_1d_total)) then
         allocate(mueller_1d_total(1:size(mueller_1d,1),1:size(mueller_1d,2)))
         mueller_1d_total = 0 ! init
      end if

      ! sum
      mueller_total = mueller_total + mueller
      mueller_1d_total = mueller_1d_total + mueller_1d

   end subroutine

    end module outputs_mod