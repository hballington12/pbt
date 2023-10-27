! outputs_mod.f90
! contains subroutines for outputting
    
   module outputs_mod
   
   use misc_submod
   
   implicit none
   
   
   contains
   
   subroutine cache_remaining_orients_seq(cache_dir,i,num_remaining_orients,remaining_orients)

      ! saves the remaining orientation numbers to cached file

      integer, intent(in) :: i
      character(len=255), intent(in) :: cache_dir ! cached files directory (if job stops early)
      integer j
      integer, dimension(:), allocatable, intent(in) :: remaining_orients
      integer, intent(in) ::  num_remaining_orients

      open(10,file=trim(cache_dir)//"/orient_remaining.dat") ! open file
         do j = i + 1, num_remaining_orients
            write(10,*) remaining_orients(j)
         end do
      close(10)      

   end subroutine
   
   subroutine finalise( ampl_far_beam11,     & ! amplitude matrix (1,1) due to beam diffraction
                        ampl_far_beam12,     & ! amplitude matrix (1,2) due to beam diffraction
                        ampl_far_beam21,     & ! amplitude matrix (2,1) due to beam diffraction
                        ampl_far_beam22,     & ! amplitude matrix (2,2) due to beam diffraction
                        ampl_far_ext_diff11, & ! amplitude matrix (1,1) due to external diffraction
                        ampl_far_ext_diff12, & ! amplitude matrix (1,2) due to external diffraction
                        ampl_far_ext_diff21, & ! amplitude matrix (2,1) due to external diffraction
                        ampl_far_ext_diff22, & ! amplitude matrix (2,2) due to external diffraction
                        energy_out_beam,     & ! beam energy remaining before diffraction
                        energy_out_ext_diff, & ! external diffraction energy remaining before diffraction
                        mueller,             & ! 2d mueller matrix
                        mueller_1d,          & ! 1d mueller matrix
                        energy_abs_beam,     & ! energy absorbed inside the particle
                        output_parameters,   & ! output parameters
                        job_params)
      
      ! sr finalise is called at the end of each rotation
      ! it combines the beam and external diffraction far-field amplitude matrices to yield the total far-field
      ! far-field amplitude matrix is then used to compute the 2d mueller matrix
      ! the 2d mueller matrix is integrated using sr simpne, which yields the 1d mueller matrices
      !  simpne interpolates a 3-point Lagrangian polynomial to the data and integrates that exactly
      ! the 1d mueller matrices are integrated to compute various integrated scattering parameters
      
      complex(8), dimension(:,:), allocatable, intent(inout) :: ampl_far_beam11, ampl_far_beam12, ampl_far_beam21, ampl_far_beam22 ! beam
      complex(8), dimension(:,:), allocatable, intent(inout)  :: ampl_far_ext_diff11, ampl_far_ext_diff12, ampl_far_ext_diff21, ampl_far_ext_diff22 ! ext diff
      real(8), dimension(:), allocatable :: theta_vals, phi_vals
      real(8), intent(in) :: energy_out_beam
      real(8), intent(in) :: energy_out_ext_diff
      real(8), dimension(:,:,:), allocatable, intent(out) :: mueller ! mueller matrices
      real(8), dimension(:,:), allocatable, intent(out) :: mueller_1d ! phi-integrated mueller matrices
      real(8) la ! wavelength
      real(8), intent(in) :: energy_abs_beam
      type(output_parameters_type), intent(out) :: output_parameters 
      type(job_parameters_type), intent(in) :: job_params
      
      complex(8), dimension(:,:), allocatable :: ampl_far11, ampl_far12, ampl_far21, ampl_far22 ! total
      real(8), dimension(:,:,:), allocatable :: mueller_beam, mueller_ext_diff ! mueller matrices
      real(8), dimension(:,:), allocatable :: mueller_beam_1d, mueller_ext_diff_1d ! phi-integrated mueller matrices

      real(8) scatt, scatt_beam, scatt_ext_diff, asymmetry, asymmetry_beam, asymmetry_ext_diff, ext, absorption, albedo, back_scatt
      real(8) waveno
      integer(8) i
      
      theta_vals = job_params%theta_vals
      phi_vals = job_params%phi_vals
      la = job_params%la
      
      waveno = 2d0*pi/la
      
      ! allocate total amplitude matrix (beam - ext diff)
      allocate(ampl_far11(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
      allocate(ampl_far12(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
      allocate(ampl_far21(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))
      allocate(ampl_far22(1:size(ampl_far_beam11,1),1:size(ampl_far_beam11,2)))

      ! account for 1/waveno factor in Bohren & Huffman eq 3.12
      ampl_far_beam11 = ampl_far_beam11 * waveno
      ampl_far_beam12 = ampl_far_beam12 * waveno
      ampl_far_beam21 = ampl_far_beam21 * waveno
      ampl_far_beam22 = ampl_far_beam22 * waveno
      ampl_far_ext_diff11 = ampl_far_ext_diff11 * waveno
      ampl_far_ext_diff12 = ampl_far_ext_diff12 * waveno
      ampl_far_ext_diff21 = ampl_far_ext_diff21 * waveno
      ampl_far_ext_diff22 = ampl_far_ext_diff22 * waveno      

      ! far field = beam diffraction - external diffraction
      ampl_far11 = ampl_far_beam11 + ampl_far_ext_diff11
      ampl_far12 = ampl_far_beam12 + ampl_far_ext_diff12
      ampl_far21 = ampl_far_beam21 + ampl_far_ext_diff21
      ampl_far22 = ampl_far_beam22 + ampl_far_ext_diff22
      
      print*,'making 2d mueller matrices...'
 
      call ampl_to_mueller(ampl_far11,ampl_far12,ampl_far21,ampl_far22,mueller)
      call ampl_to_mueller(ampl_far_beam11,ampl_far_beam12,ampl_far_beam21,ampl_far_beam22,mueller_beam)
      call ampl_to_mueller(ampl_far_ext_diff11,ampl_far_ext_diff12,ampl_far_ext_diff21,ampl_far_ext_diff22,mueller_ext_diff)

      print*,'making 1d mueller matrices...'
      write(101,*)'------------------------------------------------------'
      
      call get_1d_mueller(mueller, mueller_1d, theta_vals, phi_vals)
      call get_1d_mueller(mueller_beam, mueller_beam_1d, theta_vals, phi_vals)
      call get_1d_mueller(mueller_ext_diff, mueller_ext_diff_1d, theta_vals, phi_vals)

      print*,'calculating integrated parameters...'
      
      ! theta integrations...
      call simpne(size(theta_vals,1),theta_vals,mueller_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt) ! p11*sin(theta)
      write(101,'(A40,f16.8)')'scatt. cross (total):',scatt
      print'(A40,f16.8)','scattering cross section (total):',scatt
      
      call simpne(size(theta_vals,1),theta_vals,mueller_beam_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt_beam) ! p11*sin(theta)
      write(101,'(A40,f16.8,A2,f10.6,A3)')'scatt. cross (beam):',scatt_beam," (",scatt_beam/energy_out_beam*100," %)"
      print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (beam):',scatt_beam," (",scatt_beam/energy_out_beam*100," %)"
      
      call simpne(size(theta_vals,1),theta_vals,mueller_ext_diff_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt_ext_diff) ! p11*sin(theta)
      write(101,'(A40,f16.8,A2,f10.6,A3)')'scatt. cross (ext diff):',scatt_ext_diff," (",scatt_ext_diff/energy_out_ext_diff*100," %)"
      print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (ext diff):',scatt_ext_diff," (",scatt_ext_diff/energy_out_ext_diff*100," %)"
      
      if (job_params%scaling) then
         print'(A40)','applying diffraction energy scaling...'

         ! scale the beam diffraction energy to match near field
         ampl_far_beam11 = ampl_far_beam11/sqrt(scatt_beam/energy_out_beam)
         ampl_far_beam12 = ampl_far_beam12/sqrt(scatt_beam/energy_out_beam)
         ampl_far_beam21 = ampl_far_beam21/sqrt(scatt_beam/energy_out_beam)
         ampl_far_beam22 = ampl_far_beam22/sqrt(scatt_beam/energy_out_beam)

         ! scale the external diffraction energy to match near field
         ampl_far_ext_diff11 = ampl_far_ext_diff11/sqrt(scatt_ext_diff/energy_out_ext_diff)
         ampl_far_ext_diff12 = ampl_far_ext_diff12/sqrt(scatt_ext_diff/energy_out_ext_diff)
         ampl_far_ext_diff21 = ampl_far_ext_diff21/sqrt(scatt_ext_diff/energy_out_ext_diff)
         ampl_far_ext_diff22 = ampl_far_ext_diff22/sqrt(scatt_ext_diff/energy_out_ext_diff)

         ampl_far11 = ampl_far_beam11 + ampl_far_ext_diff11
         ampl_far12 = ampl_far_beam12 + ampl_far_ext_diff12
         ampl_far21 = ampl_far_beam21 + ampl_far_ext_diff21
         ampl_far22 = ampl_far_beam22 + ampl_far_ext_diff22
         
         print*,'remaking 2d mueller matrices...'
    
         call ampl_to_mueller(ampl_far11,ampl_far12,ampl_far21,ampl_far22,mueller)
         call ampl_to_mueller(ampl_far_beam11,ampl_far_beam12,ampl_far_beam21,ampl_far_beam22,mueller_beam)
         call ampl_to_mueller(ampl_far_ext_diff11,ampl_far_ext_diff12,ampl_far_ext_diff21,ampl_far_ext_diff22,mueller_ext_diff)
   
         print*,'remaking 1d mueller matrices...'
         write(101,*)'------------------------------------------------------'
         
         call get_1d_mueller(mueller, mueller_1d, theta_vals, phi_vals)
         call get_1d_mueller(mueller_beam, mueller_beam_1d, theta_vals, phi_vals)
         call get_1d_mueller(mueller_ext_diff, mueller_ext_diff_1d, theta_vals, phi_vals)
   
         print*,'recalculating integrated parameters...'
         
         ! theta integrations...
         call simpne(size(theta_vals,1),theta_vals,mueller_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt) ! p11*sin(theta)
         write(101,'(A40,f16.8)')'scatt. cross (total):',scatt
         print'(A40,f16.8)','scattering cross section (total):',scatt
         
         call simpne(size(theta_vals,1),theta_vals,mueller_beam_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt_beam) ! p11*sin(theta)
         write(101,'(A40,f16.8,A2,f10.6,A3)')'scatt. cross (beam):',scatt_beam," (",scatt_beam/energy_out_beam*100," %)"
         print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (beam):',scatt_beam," (",scatt_beam/energy_out_beam*100," %)"
         
         call simpne(size(theta_vals,1),theta_vals,mueller_ext_diff_1d(1:size(theta_vals,1),1)*sin(theta_vals)/(waveno**2),scatt_ext_diff) ! p11*sin(theta)
         write(101,'(A40,f16.8,A2,f10.6,A3)')'scatt. cross (ext diff):',scatt_ext_diff," (",scatt_ext_diff/energy_out_ext_diff*100," %)"
         print'(A40,f16.8,A2,f10.6,A3)','scattering cross section (ext diff):',scatt_ext_diff," (",scatt_ext_diff/energy_out_ext_diff*100," %)"








      end if

      ! ext = absorption(2*pi/waveno*real(ampl_far11(1,1) + ampl_far12(1,1) + ampl_far21(1,1) + ampl_far22(1,1))) ! extinction cross section
      ! write(101,'(A40,f16.8)'),'ext. cross section via opt. theorem (real):',ext 
      ! print'(A40,f16.8)','ext. cross (opt. theorem) real:',ext 
      
      ! ext = absorption(2*pi/waveno*imag(ampl_far11(1,1) + ampl_far12(1,1) + ampl_far21(1,1) + ampl_far22(1,1))) ! extinction cross section
      ! write(101,'(A40,f16.8)'),'ext. cross section via opt. theorem (imag):',ext 
      ! print'(A40,f16.8)','ext. cross (opt. theorem) imag:',ext 
      
      absorption = energy_abs_beam
      write(101,'(A40,f16.8)')'abs. cross section:',absorption
      print'(A40,f16.8)','abs. cross:',absorption 
      
      ext = absorption + scatt
      write(101,'(A40,f16.8)')'ext. cross section:',ext
      print'(A40,f16.8)','ext. cross:',ext 
      
      albedo = 1-(ext-scatt)/ext 
      write(101,'(A40,f16.8)')'single-scattering albedo:',albedo
      print'(A40,f16.8)','single-scatt. albedo:',albedo
      
      call simpne(size(theta_vals,1),theta_vals,mueller_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt/(waveno**2),asymmetry) ! p11*sin(theta)
      write(101,'(A40,f16.8)')'asymmetry parameter (total):',asymmetry
      print'(A40,f16.8)','asymmetry parameter (total):',asymmetry
      
      call simpne(size(theta_vals,1),theta_vals,mueller_beam_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt_beam/(waveno**2),asymmetry_beam) ! p11*sin(theta)
      write(101,'(A40,f16.8)')'asymmetry parameter (beam):',asymmetry_beam
      print'(A40,f16.8)','asymmetry parameter (beam):',asymmetry_beam
      
      call simpne(size(theta_vals,1),theta_vals,mueller_ext_diff_1d(1:size(theta_vals,1),1)*sin(theta_vals)*cos(theta_vals)/scatt_ext_diff/(waveno**2),asymmetry_ext_diff) ! p11*sin(theta)
      write(101,'(A40,f16.8)')'asymmetry parameter (ext diff):',asymmetry_ext_diff  
      print'(A40,f16.8)','asymmetry parameter (ext diff):',asymmetry_ext_diff  
      
      ! calculate back-scattering cross section
      ! find closest theta value to direct back-scattering
      i = minloc(abs(theta_vals - 2*pi),1)
      ! print*,'i = ',i,'(',theta_vals(i)*180D0/pi,')' ! print the value of theta at direct back-scattering
      back_scatt = mueller(1,i,1)*4*pi/waveno**2 ! calculate back-scattering cross section (B&H definition)
      ! back_scatt = mueller(1,i,1)/waveno**2 ! calculate back-scattering cross section (intuitive definition)
      write(101,'(A40,f16.8)')'back scatt. cross section:',back_scatt
      print'(A40,f16.8)','back scatt. cross:',back_scatt

      output_parameters%scatt = scatt
      output_parameters%abs = absorption
      output_parameters%ext = ext
      output_parameters%albedo = albedo
      output_parameters%asymmetry = asymmetry
      output_parameters%abs_eff = absorption / output_parameters%geo_cross_sec
      output_parameters%scatt_eff = scatt / output_parameters%geo_cross_sec
      output_parameters%ext_eff = ext / output_parameters%geo_cross_sec
      output_parameters%back_scatt = back_scatt



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
   
   subroutine get_1d_mueller(mueller, mueller_1d, theta_vals, phi_vals)

      ! integrates over phi at each theta value to get the 1d mueller matrix from a 2d mueller matrix

      real(8), dimension(:,:,:), allocatable, intent(in) :: mueller
      real(8), dimension(:,:), allocatable, intent(out) :: mueller_1d
      real(8), dimension(:), allocatable, intent(in) :: theta_vals, phi_vals

      integer(8) i

      ! allocate 1d mueller matrix
      allocate(mueller_1d(1:size(theta_vals,1),1:16)) ! 1:1 is for each element

! phi integrations...
      do i = 1, size(theta_vals,1) ! for each theta bin...
         ! p11
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,1),mueller_1d(i,1)) ! integrate total
         ! p12
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,2),mueller_1d(i,2)) ! integrate total
         ! p13
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,3),mueller_1d(i,3)) ! integrate total
         ! p14
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,4),mueller_1d(i,4)) ! integrate total
         ! p21
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,5),mueller_1d(i,5)) ! integrate total
         ! p22
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,6),mueller_1d(i,6)) ! integrate total
         ! p23
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,7),mueller_1d(i,7)) ! integrate total
         ! p24
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,8),mueller_1d(i,8)) ! integrate total
         ! p31
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,9),mueller_1d(i,9)) ! integrate total
         ! p32
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,10),mueller_1d(i,10)) ! integrate total
         ! p33
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,11),mueller_1d(i,11)) ! integrate total
         ! p34
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,12),mueller_1d(i,12)) ! integrate total
         ! p41
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,13),mueller_1d(i,13)) ! integrate total
         ! p42
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,14),mueller_1d(i,14)) ! integrate total
         ! p43
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,15),mueller_1d(i,15)) ! integrate total
         ! p44
         call simpne(size(phi_vals,1),phi_vals,mueller(1:size(mueller,1),i,16),mueller_1d(i,16)) ! integrate total
      end do

   end subroutine

   subroutine writeup(  mueller,    &
      mueller_1d, &
      output_dir, &
      output_parameters_total, &
      job_params)
      
      ! sr writeup writes the 1d and 2d mueller matrices to the job directory
      
      real(8), dimension(:,:,:), allocatable, intent(in) :: mueller ! mueller matrices
      real(8), dimension(:,:), allocatable, intent(in) :: mueller_1d ! phi-integrated mueller matrices
      real(8), dimension(:), allocatable :: theta_vals, phi_vals
      character(len=*), intent(in) :: output_dir
      type(output_parameters_type), intent(inout) :: output_parameters_total
      type(job_parameters_type), intent(in) :: job_params
      
      integer i, j
      
      theta_vals = job_params%theta_vals
      phi_vals = job_params%phi_vals
      
      print*,'writing mueller to file...'
      
      if(.not. job_params%suppress_2d) then
         open(10,file=trim(output_dir)//"/"//"mueller_scatgrid")
         do i = 1, size(theta_vals,1)
            do j = 1, size(phi_vals,1)
               write(10,'(f12.4,f12.4,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8)') &
               theta_vals(i)*180/pi, phi_vals(j)*180/pi, &
               mueller(j,i,1), mueller(j,i,2), mueller(j,i,3), mueller(j,i,4), &
               mueller(j,i,5), mueller(j,i,6), mueller(j,i,7), mueller(j,i,8), &
               mueller(j,i,9), mueller(j,i,10), mueller(j,i,11), mueller(j,i,12), &
               mueller(j,i,13), mueller(j,i,14), mueller(j,i,15), mueller(j,i,16)                                                                             
            end do
         end do
         close(10)
      end if
      open(10,file=trim(output_dir)//"/"//"mueller_scatgrid_1d")
      do j = 1, size(theta_vals,1)
         write(10,'(f12.4,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8,f22.8)') &
         theta_vals(j)*180/pi, &
         mueller_1d(j,1), mueller_1d(j,2), mueller_1d(j,3), mueller_1d(j,4), &
         mueller_1d(j,5), mueller_1d(j,6), mueller_1d(j,7), mueller_1d(j,8), &
         mueller_1d(j,9), mueller_1d(j,10), mueller_1d(j,11), mueller_1d(j,12), &
         mueller_1d(j,13), mueller_1d(j,14), mueller_1d(j,15), mueller_1d(j,16)   
      end do
      close(10)
      
      open(10,file=trim(output_dir)//"/"//"params")
      write(10,*) 'scattering parameters (orientation averaged)...'
      write(10,'(A30,f16.8)') 'geo. cross section: ',output_parameters_total%geo_cross_sec
      write(10,'(A30,f16.8)') 'abs. cross section: ',output_parameters_total%abs
      write(10,'(A30,f16.8)') 'scatt. cross section: ',output_parameters_total%scatt
      write(10,'(A30,f16.8)') 'ext. cross section: ',output_parameters_total%ext
      write(10,'(A30,f16.8)') 'back scatt. cross section: ',output_parameters_total%back_scatt
      write(10,'(A30,f16.8)') 'single scattering albedo: ',output_parameters_total%albedo
      write(10,'(A30,f16.8)') 'asymmetry parameter: ',output_parameters_total%asymmetry
      write(10,'(A30,f16.8)') 'abs. efficiency: ',output_parameters_total%abs / output_parameters_total%geo_cross_sec
      write(10,'(A30,f16.8)') 'scatt. efficiency: ',output_parameters_total%scatt / output_parameters_total%geo_cross_sec
      write(10,'(A30,f16.8)') 'ext. efficiency: ',output_parameters_total%ext / output_parameters_total%geo_cross_sec
      ! below is commented because: different orientations have dif. geo. cross sections, so the efficiencies dont add linearly
      ! write(10,'(A30,f16.8)') 'abs. efficiency: ',output_parameters_total%abs_eff
      ! write(10,'(A30,f16.8)') 'scatt. efficiency: ',output_parameters_total%scatt_eff
      ! write(10,'(A30,f16.8)') 'ext. efficiency: ',output_parameters_total%ext_eff
      close(10)
      
      ! print*,'finished mueller to file.'
      
   end subroutine
   
   subroutine summation(mueller,                & ! current 2d mueller
      mueller_total,          & ! total 2d mueller
      mueller_1d,             & ! current 1d mueller
      mueller_1d_total,       & ! total 1d mueller
      output_parameters,      & 
      output_parameters_total)
      
      ! sr summation adds the current mueller matrices to the total
      
      real(8), dimension(:,:,:), allocatable, intent(in) :: mueller ! mueller matrices
      real(8), dimension(:,:,:), allocatable, intent(inout) :: mueller_total ! mueller matrices
      real(8), dimension(:,:), allocatable, intent(in) :: mueller_1d ! phi-integrated mueller matrices
      real(8), dimension(:,:), allocatable, intent(inout) :: mueller_1d_total ! phi-integrated mueller matrices
      type(output_parameters_type), intent(in) :: output_parameters 
      type(output_parameters_type), intent(inout) :: output_parameters_total
      
      ! if its the first call to summation, allocate the total mueller 1d and 2d arrays
      if(.not. allocated(mueller_total)) then
         allocate(mueller_total(1:size(mueller,1),1:size(mueller,2),1:size(mueller,3)))
         mueller_total = 0 ! init
         output_parameters_total%abs = 0 ! init
         output_parameters_total%scatt = 0 ! init
         output_parameters_total%ext = 0 ! init
         output_parameters_total%albedo = 0 ! init
         output_parameters_total%asymmetry = 0 ! init
         output_parameters_total%abs_eff = 0 ! init
         output_parameters_total%scatt_eff = 0 ! init
         output_parameters_total%ext_eff = 0 ! init
         output_parameters_total%geo_cross_sec = 0 ! init
         output_parameters_total%back_scatt = 0 ! init
      end if
      if(.not. allocated(mueller_1d_total)) then
         allocate(mueller_1d_total(1:size(mueller_1d,1),1:size(mueller_1d,2)))
         mueller_1d_total = 0 ! init
      end if
      
      ! sum
      mueller_total = mueller_total + mueller
      mueller_1d_total = mueller_1d_total + mueller_1d
      output_parameters_total%abs = output_parameters_total%abs + output_parameters%abs
      output_parameters_total%scatt = output_parameters_total%scatt + output_parameters%scatt
      output_parameters_total%back_scatt = output_parameters_total%back_scatt + output_parameters%back_scatt
print*,'output_parameters_total%back_scatt =',output_parameters_total%back_scatt
      output_parameters_total%ext = output_parameters_total%ext + output_parameters%ext
      output_parameters_total%albedo = output_parameters_total%albedo + output_parameters%albedo
      output_parameters_total%asymmetry = output_parameters_total%asymmetry + output_parameters%asymmetry
      output_parameters_total%abs_eff = output_parameters_total%abs_eff + output_parameters%abs_eff
      output_parameters_total%scatt_eff = output_parameters_total%scatt_eff + output_parameters%scatt_eff
      output_parameters_total%ext_eff = output_parameters_total%ext_eff + output_parameters%ext_eff
      output_parameters_total%geo_cross_sec = output_parameters_total%geo_cross_sec + output_parameters%geo_cross_sec
      
   end subroutine
   
   subroutine ampl_to_mueller(ampl11,ampl12,ampl21,ampl22,mueller)

      ! returns a mueller matrix from an amplitude matrix

      complex(8), dimension(:,:), allocatable, intent(in) :: ampl11, ampl12, ampl21, ampl22 
      real(8), dimension(:,:,:), allocatable, intent(out) :: mueller

      integer(8) i, j

      ! allocate mueller matrix
      allocate(mueller(1:size(ampl11,1),1:size(ampl11,2),1:16)) ! 1:16 is for each element

      do i = 1, size(ampl11,2) ! loop over theta
         do j = 1, size(ampl11,1) ! loop over phi
            ! Bohren & Huffman
            ! ##########################--S11--##########################
            mueller(j,i,1) = real(0.5*(ampl11(j,i)*conjg(ampl11(j,i)) + &
                                       ampl12(j,i)*conjg(ampl12(j,i)) + &
                                       ampl21(j,i)*conjg(ampl21(j,i)) + &
                                       ampl22(j,i)*conjg(ampl22(j,i))))   
            ! ##########################--S12--##########################                                                
            mueller(j,i,2) = real(0.5*(ampl11(j,i)*conjg(ampl11(j,i)) - &
                                       ampl12(j,i)*conjg(ampl12(j,i)) + &
                                       ampl21(j,i)*conjg(ampl21(j,i)) - &
                                       ampl22(j,i)*conjg(ampl22(j,i)))) 
            ! ##########################--S13--##########################
            mueller(j,i,3) =      real(ampl11(j,i)*conjg(ampl12(j,i)) + &
                                       ampl22(j,i)*conjg(ampl21(j,i))) 
            ! ##########################--S14--##########################
            mueller(j,i,4) =      imag(ampl11(j,i)*conjg(ampl12(j,i)) - &
                                       ampl22(j,i)*conjg(ampl21(j,i)))
            ! ##########################--S21--##########################
            mueller(j,i,5) = real(0.5*(ampl11(j,i)*conjg(ampl11(j,i)) + &
                                       ampl12(j,i)*conjg(ampl12(j,i)) - &
                                       ampl21(j,i)*conjg(ampl21(j,i)) - &
                                       ampl22(j,i)*conjg(ampl22(j,i))))
            ! ##########################--S22--##########################                                                             
            mueller(j,i,6) = real(0.5*(ampl11(j,i)*conjg(ampl11(j,i)) - &
                                       ampl12(j,i)*conjg(ampl12(j,i)) - &
                                       ampl21(j,i)*conjg(ampl21(j,i)) + &
                                       ampl22(j,i)*conjg(ampl22(j,i))))
            ! ##########################--S23--##########################   
            mueller(j,i,7) =      real(ampl11(j,i)*conjg(ampl12(j,i)) - &
                                       ampl22(j,i)*conjg(ampl21(j,i)))
            ! ##########################--S24--##########################   
            mueller(j,i,8) =      imag(ampl11(j,i)*conjg(ampl12(j,i)) + &
                                       ampl22(j,i)*conjg(ampl21(j,i))) 
            ! ##########################--S31--##########################   
            mueller(j,i,9) =      real(ampl11(j,i)*conjg(ampl21(j,i)) + &
                                       ampl22(j,i)*conjg(ampl12(j,i)))
            ! ##########################--S32--##########################  
            mueller(j,i,10) =     real(ampl11(j,i)*conjg(ampl21(j,i)) - &
                                       ampl22(j,i)*conjg(ampl12(j,i))) 
            ! ##########################--S33--##########################  
            mueller(j,i,11) =     real(ampl11(j,i)*conjg(ampl22(j,i)) + &
                                       ampl12(j,i)*conjg(ampl21(j,i)))  
            ! ##########################--S34--##########################   
            mueller(j,i,12) =     imag(ampl11(j,i)*conjg(ampl22(j,i)) + &
                                       ampl21(j,i)*conjg(ampl12(j,i))) 
            ! ##########################--S41--########################## 
            mueller(j,i,13) =     imag(ampl21(j,i)*conjg(ampl11(j,i)) + &
                                       ampl22(j,i)*conjg(ampl12(j,i)))  
            ! ##########################--S42--##########################  
            mueller(j,i,14) =     imag(ampl21(j,i)*conjg(ampl11(j,i)) - &
                                       ampl22(j,i)*conjg(ampl12(j,i)))  
            ! ##########################--S43--##########################  
            mueller(j,i,15) =     imag(ampl22(j,i)*conjg(ampl11(j,i)) - &
                                       ampl12(j,i)*conjg(ampl21(j,i)))  
            ! ##########################--S44--##########################   
            mueller(j,i,16) =     real(ampl22(j,i)*conjg(ampl11(j,i)) - &
                                       ampl12(j,i)*conjg(ampl21(j,i)))     
            
            
         end do
      end do

   end subroutine

   subroutine cache_job(vert_in,                    & ! unrotated vertices
                        face_ids,                   & ! face ids
                        num_face_vert,              & ! num vertices per face
                        apertures,                  & ! apertures
                        job_params,                 & ! job parameters
                        i_loop,                     & ! current loop index
                        output_parameters_total,    & ! total output parameters
                        mueller_total,              & ! total 2d mueller
                        mueller_1d_total,           & ! total 1d mueller
                        cache_dir)
      
      ! saves the job to the cache directory, possibly to be resumed later
      
      real(8), dimension(:,:), allocatable, intent(in) :: vert_in ! unique vertices (unrotated)
      integer, dimension(:,:), allocatable, intent(in) :: face_ids ! face vertex IDs
      integer, dimension(:), allocatable, intent(in) :: num_face_vert ! number of vertices in each face
      integer, dimension(:), allocatable, intent(in) :: apertures ! apertures asignments for each facet
      type(job_parameters_type), intent(in) :: job_params ! job parameters, contains wavelength, rbi, etc., see types mod for more details
      integer, intent(in) :: i_loop
      type(output_parameters_type), intent(inout) :: output_parameters_total
      real(8), dimension(:,:,:), allocatable , intent(in):: mueller_total ! mueller matrices
      real(8), dimension(:,:), allocatable, intent(in) :: mueller_1d_total ! phi-integrated mueller matrices
      character(len=255), intent(in) :: cache_dir ! cached files directory (if job stops early)
      
      print*,'job exit point detected. saving job files to cache...'
      ! call make_cache_dir("cache/",cache_dir)
      print*,'cache directory is "',trim(cache_dir),'"'
      
      call PDAS(  vert_in,        & ! <-  rotated vertices
      face_ids,       & ! <-  face vertex IDs
      cache_dir,      & ! <-  output directory
      num_face_vert,  & ! <-  number of verices in each face
      "unrotated")      ! <-  filename
      
      call save_apertures(apertures,  & ! <-  apertures
      cache_dir)    ! <-  cache directory
      
      call save_params(job_params,i_loop,cache_dir,output_parameters_total)
      
      call writeup(mueller_total, mueller_1d_total, cache_dir, output_parameters_total, job_params) ! write to file
      
      print*,'saved job files to cache.'
      print*,'to resume this job, include the "-resume '//trim(cache_dir(7:len(cache_dir)))//'" flag when you call abt.'
            
   end subroutine
   
end module outputs_mod