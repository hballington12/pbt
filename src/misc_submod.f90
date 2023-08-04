! misc_submod.f90
! submodule for misc subroutines and functions used across all modules and the main program
! also contains fundamental constants
    
module misc_submod

use ifport

implicit none

real(8), parameter, public :: pi = 3.14159265358979 

contains

subroutine make_dir(dir_path_in,cwd_out)

    ! makes job directory - quick and dirty method

    character(len=*), intent(in) :: dir_path_in
    character(len=255) dir_path
    character(len=255), intent(out) :: cwd_out
    logical exists ! for checking if directories or files exist
    integer result ! true if subdirectory was made, false if subdirectory was not made
    integer job_num ! job number
    character(len=255) job_num_string ! job number

    job_num = 1
    result = .false.

    do while (result .eq. .false.) ! while a new directory has not been made

        ! append job_num to directory name
        write(job_num_string,"(I)") job_num
        call StripSpaces(job_num_string)
        ! print*,'job_num_string:',trim(job_num_string)
        ! print*,'dir_path: ',trim(dir_path_in)//trim(job_num_string)
        write(dir_path,*) trim(dir_path_in)//trim(job_num_string)
        call StripSpaces(dir_path)
        ! print*,'attempting to make directory with name "',trim(dir_path),'"'

        inquire(directory = trim(dir_path), exist = exists)
        if(exists .eq. .false.) then 
            ! print*,'Creating job directory at "',trim(dir_path),'"'
            write(cwd_out,*) trim(dir_path)
            result = makedirqq(trim(dir_path))
        else
            ! if already exists, add 1 to job number and try again
            job_num = job_num + 1
            ! print*,'error: job directory already exists'
            ! stop
        end if

    end do

end subroutine

subroutine PDAS(verts,face_ids,output_dir)

    ! writes rotated particle to file
    ! to do: add Macke-style output

real(8), dimension(:,:), allocatable, intent(in) :: verts ! unique vertices
integer(8), dimension(:,:), allocatable, intent(in) :: face_ids ! face vertex IDs
character(len=*), intent(in) :: output_dir

integer num_verts, num_faces, i, j
character(100) my_string, my_string2

print*,'========== end sr PDAS'

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

open(10,file=trim(output_dir)//"/"//"rotated_particle.obj")
write(10,*) '# rotated particle'
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

    write(10,'(A1,f12.8,f12.8,f12.8)') "v ", verts(i,1), verts(i,2), verts(i,3)
end do
do i = 1, num_faces
    write(10,'(A1,I6,I6,I6)') 'f', face_ids(i,1), face_ids(i,2), face_ids(i,3)
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

real(8), dimension(:), allocatable, intent(in) :: array_1, array_2
real(8), dimension(:,:), allocatable, intent(out) :: mesh_1, mesh_2

integer(8) array_1_dim, array_2_dim
integer(8) i, j

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
    ! wrote this in a giffy so its messy but will clean up later

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

    ! var1 is what is to be concatenated
    ! var2 is what is to be concatenated on to
    ! version for integer arrays

integer(8), dimension(:,:), allocatable, intent(inout) :: var1
integer(8), dimension(:,:), allocatable, intent(inout) :: var2

integer(8), dimension(:,:), allocatable :: var1_temp
integer(8), dimension(:,:), allocatable :: var2_temp
integer(8) var1_size1, var1_size2 ! array dimension sizes
integer(8) var2_size1, var2_size2
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

    ! var1 is what is to be concatenated
    ! var2 is what is to be concatenated on to
    ! version for reall arrays

real(8), dimension(:,:), allocatable, intent(inout) :: var1
real(8), dimension(:,:), allocatable, intent(inout) :: var2

real(8), dimension(:,:), allocatable :: var1_temp
real(8), dimension(:,:), allocatable :: var2_temp
integer(8) var1_size1, var1_size2 ! array dimension sizes
integer(8) var2_size1, var2_size2
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

    ! var1 is what is to be concatenated
    ! var2 is what is to be concatenated on to
    ! version for complex arrays

complex(8), dimension(:,:), allocatable, intent(inout) :: var1
complex(8), dimension(:,:), allocatable, intent(inout) :: var2

complex(8), dimension(:,:), allocatable :: var1_temp
complex(8), dimension(:,:), allocatable :: var2_temp
integer(8) var1_size1, var1_size2 ! array dimension sizes
integer(8) var2_size1, var2_size2
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

! trims an integer(8) 1d array according to a logical mask

logical, dimension(:), allocatable, intent(in) :: mask
integer(8), dimension(:), allocatable, intent(inout) :: array

integer(8), dimension(:), allocatable :: array_temp
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
! input and output should be of type "integer(8), dimension(:), allocatable"
! inspired by: https://stackoverflow.com/questions/44198212/a-fortran-equivalent-to-unique

integer(8), dimension(:), allocatable, intent(in) :: array_in
integer(8), dimension(:), allocatable :: unique
integer(8), dimension(:), allocatable, intent(out) :: array_unique
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

exp2cmplx = cmplx(cos(arg),sin(arg))

end function

subroutine StripSpaces(string) ! taken from: https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
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

subroutine midPointsAndAreas(face_ids, verts, midpoints, face_areas)

! subroutine midPointsAndAreas computes midpoints and areas of each facet (assumes triangle facets)

integer(8), dimension(:,:) ,allocatable, intent(in) :: face_ids
real(8), dimension(:,:), allocatable, intent(in) :: verts
real(8), dimension(:,:), allocatable, intent(out) :: midpoints
real(8), dimension(:), allocatable, intent(out) :: face_areas
real(8), dimension(1:3) :: vec_a, vec_b, a_cross_b
real(8) temp_area

integer i,j, num_verts

!print*,'========== start sr midPointsAndAreas'

allocate(midpoints(size(face_ids,1),3))

num_verts = size(face_ids,2) ! assumes all faces have the same number of vertices
!print*,'num verts per face: ',num_verts
    
do i = 1, size(face_ids,1) ! for each face
    midpoints(i,1:3) = sum(verts(face_ids(i,1:num_verts),1:3),1)/num_verts
end do

! allocate face_areas array
allocate(face_areas(1:size(face_ids,1)))
face_areas = 0 ! initialise

do i = 1, size(face_ids,1) ! for each face
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
    if(normalise .eq. .false.) nf = 1
end if

c = c/nf


end subroutine

subroutine surfSubdivide(VV, FF1)

! subroutine surfSubdivide divides each facet of surface mesh into 4 smaller ones
! each new facet is made from the midpoint of the facet,  2 edge midpoints, and an existing vertex
! input mesh is stored in a temporary array, deallocated, then rellocated with new array dimensions using Euler's characteristic formula

real(8), dimension(:,:), allocatable, intent(inout) :: VV ! input and output vertex array
integer(8), dimension(:,:), allocatable, intent(inout) :: FF1 ! input and output face array

real(8), dimension(:,:), allocatable :: V ! temporary array to hold input vertex array
integer(8), dimension(:,:), allocatable :: F1 ! temporary array to hold input face array
integer(8), dimension(:,:,:), allocatable :: E ! temporary array to hold edge array
integer(8) i,j,k, vertOutCount, faceOutCount ! counters
integer(8) num_in_verts, num_out_verts, num_in_faces, num_out_faces, numberOfVerticesIn
integer(8) midPointID, fwdEdgeID, fwdMidPointId, bwdEdgeID, bwdMidPointId
integer(8), dimension(:), allocatable :: vertexIDs
real(8), dimension(:,:), allocatable :: faceVerticesIn ! temporary array to hold input vertex array
real(8) midPointFace(1:3), fwdMidPoint(1:3), bwdMidPoint(1:3)

num_in_verts = size(VV,1)
num_in_faces = size(FF1,1)

print*,'number of input vertices: ',num_in_verts
print*,'number of input faces: ',num_in_faces
print*,'max verts per face: ',size(FF1,2)

! allocate temporary arrays to hold the input data
allocate(V(num_in_verts,1:3))
allocate(F1(num_in_faces,size(FF1,2)))

! hold the input data in the temporarry arrays
do i = 1,num_in_verts ! for each input vertex
    do j = 1,3 ! for each xyz coord
        V(i,j) = VV(i,j) ! save vertex array
    end do
end do
do i = 1,num_in_faces ! for each face
    do j = 1,size(FF1,2) ! for each vertex in face
        F1(i,j) = FF1(i,j)
    end do
end do

! deallocate input/output arrays
deallocate(VV)
deallocate(FF1)

! calculate number of output vertices
! doesnt account for shared vertices at the moment but doesnt really matter
num_out_faces = num_in_faces*4
num_out_verts = 2*(num_in_verts+num_in_faces+1)
num_out_verts = num_in_faces*(3*num_in_verts+1)

print*,'number of output vertices: ',num_out_verts
print*,'number of output faces: ',num_out_faces

! reallocate input/output arrays
allocate(VV(num_out_verts,1:3))
allocate(FF1(num_out_faces,1:4))

! allocate array to hold edge information
allocate(E(1:2,size(F1,2),num_in_faces))

! get edge information
allocate(vertexIDs(size(F1,2)))
do i = 1,num_in_faces ! for each face
    numberOfVerticesIn = size(F1,2) ! get the number of vertices in the face
    vertexIDs(1:numberOfVerticesIn) = F1(i,1:numberOfVerticesIn) ! save the vertex IDs
    ! for each vertex in a polygon, there is also edge
    ! for each vertex, save the ID of the vertex and the one in front of it
    do j = 1, numberOfVerticesIn
        if(j .eq. numberOfVerticesIn) then
            E(1,j,i) = vertexIDs(j)
            E(2,j,i) = vertexIDs(1)            
        else
            E(1,j,i) = vertexIDs(j)
            E(2,j,i) = vertexIDs(j+1)           
        end if
    end do
end do

VV(1:num_in_verts,1:3) = V(1:num_in_verts,1:3) ! store old vertices in the new array
vertOutCount = num_in_verts ! vertex out counter
faceOutCount = 0 ! face out counter

! do the subdivision
allocate(faceVerticesIn(size(F1,2),3))
do i = 1, num_in_faces ! for each input face, subdivide into 4 new faces
    numberOfVerticesIn = size(F1,2) ! get the number of vertices in the face
    vertexIDs(1:numberOfVerticesIn) = F1(i,1:numberOfVerticesIn) ! save the vertex IDs
    do j = 1, numberOfVerticesIn ! for each vertex in the face
        faceVerticesIn(j,1:3) = V(vertexIDs(j),1:3) ! retrieve the vertex coordinates
    end do
    midPointFace(1:3) = sum(faceVerticesIn,1)/numberOfVerticesIn ! compute the midpoint of the face
    vertOutCount = vertOutCount + 1 ! update counter
    midPointID = vertOutCount ! keep track of the midpoint ID for this face
    VV(vertOutCount,1:3) = midPointFace(1:3) ! add the face midpoint to the out-vertices array
    
    do j = 1, numberOfVerticesIn ! for each vertex in the face
        fwdEdgeID = E(2,j,i); ! get the vertex ID of the vertex in front
        fwdMidPoint(1:3) = (V(fwdEdgeID,1:3) + V(vertexIDs(j),1:3)) / 2;
        vertOutCount = vertOutCount + 1; ! update counter
        fwdMidPointId = vertOutCount; ! keep track of the fwd midpoint ID for this face
        VV(vertOutCount,1:3) = fwdMidPoint(1:3); ! add the edge midpoint to the out-vertices array
        if(j .eq. 1) then ! f its the first vertex, need to look at the last as the one behind
            bwdEdgeID = E(1,numberOfVerticesIn,i); ! get the vertex ID of the vertex behind
        else
            bwdEdgeID = E(1,j-1,i); ! get the vertex ID of the vertex behind
        end if
        bwdMidPoint(1:3) = (V(bwdEdgeID,1:3) + V(vertexIDs(j),1:3)) / 2;
        vertOutCount = vertOutCount + 1; ! update counter
        bwdMidPointId = vertOutCount; ! keep track of the fwd midpoint ID for this face
        VV(vertOutCount,1:3) = bwdMidPoint(1:3); ! add the edge midpoint to the out-vertices array
        
        ! now construct a new sub-face from these vertices
        faceOutCount = faceOutCount + 1; ! update counter
        FF1(faceOutCount,1) = vertexIDs(j); ! first vertexID is the original
        FF1(faceOutCount,2) = bwdMidPointId; ! second vertexID is the backward midpoint
        FF1(faceOutCount,3) = midPointID; ! third vertexID is the face midpoint
        FF1(faceOutCount,4) = fwdMidPointId; ! first vertexID is the forward midpoint        
    end do
    
end do

end subroutine

subroutine simpne ( ntab, x, y, result )

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

    SUBROUTINE PROUST(T)	!Remembrance of time passed.

        ! taken from: https://rosettacode.org/wiki/Convert_seconds_to_compound_duration#Fortran

        INTEGER T		!The time, in seconds. Positive only, please.
        INTEGER NTYPES		!How many types of time?
        PARAMETER (NTYPES = 5)	!This should do.
        INTEGER USIZE(NTYPES)	!Size of the time unit.
        CHARACTER*3 UNAME(NTYPES)!Name of the time unit.
        PARAMETER (USIZE = (/7*24*60*60, 24*60*60, 60*60,   60,    1/))	!The compiler does some arithmetic.
        PARAMETER (UNAME = (/      "wk",      "d",  "hr","min","sec"/))	!Approved names, with trailing spaces.
        CHARACTER*28 TEXT	!A scratchpad.
        INTEGER I,L,N,S		!Assistants.
         S = T			!A copy I can mess with.
         L = 0			!No text has been generated.
         DO I = 1,NTYPES		!Step through the types to do so.
           N = S/USIZE(I)	!Largest first.
           IF (N.GT.0) THEN	!Above the waterline?
             S = S - N*USIZE(I)		!Yes! Remove its contribution.
             IF (L.GT.0) THEN		!Is this the first text to be rolled?
               L = L + 2				!No.
               TEXT(L - 1:L) = ", "		!Cough forth some punctuation.
             END IF			!Now ready for this count.
             WRITE (TEXT(L + 1:),1) N,UNAME(I)	!Place, with the unit name.
     1       FORMAT (I0,1X,A)		!I0 means I only: variable-length, no leading spaces.
             L = LEN_TRIM(TEXT)		!Find the last non-blank resulting.
           END IF			!Since I'm not keeping track.
         END DO			!On to the next unit.
!  Cast forth the result.
         print*, ">> ",TEXT(1:L)," <<"	!With annotation.
        END			!Simple enough with integers.
 

end module misc_submod