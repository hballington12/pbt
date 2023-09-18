# Fortran Console Application : ABT Project Overview

This file contains a summary of the aperture beam tracer (ABT), a physical optics hybrid code developed by Harry Ballington and Evelyn Hesse.

Compile instructions (gcc-13.0.1 & openmpi-4.0.5)

- `cd src; make` - compiles sequential code
- `cd src; make abt_mpi` - compiles mpi code

Example shell scripts for submitting jobs may be found in `template/`

## `./src/seq/abt.f90`

Main source file. It contains the program entry point. The ABT reads input parameters from the command line. There are a few required input arguments, as well as a few optional ones. The code will stop if a mandatory argument is missing from the command line. Each input parameter is defined by a keyphrase, with arguments following a space delimiter. The ordering of most arguments in the command line is not important and the ABT will search through lines from left to right in an attempt to find each argument. A summary of command line arguments is given below:

- `-lambda` `<value>` - Defines the wavelength of incident light. Is required.

- `-rbi` `<value>` - Defines the real component of the particle refractive index. Is required.

- `-ibi` `<value>` - Defines the imaginary component of the particle refractive index. Is required.

- `-cmethod` `<string>` - Defines the method of particle input. Is required. Current supported methods are:

    - `read` - attempts to read the particle from the current directory. If "read" is specified, the following arguments must also be specified:
    
        - `cft` `<string>` - Defines the particle input filetype. The supported particle file input types are:
        
            - `obj` - wavefront style geometry file
        
            - `mrt` - Macke ray-tracing style
        
        - `cfn` `<string>` - Defines the particle filename. Is required. Currently, input file must be a wavefront geometry file. The particle should be triangulated.
    
        - `afn` `<string>` - Defines the apertures filename. Is required. The apertures file containes a single column defining which aperture each face belongs to. The number of lines in the apertures file must match the total number of faces in the particle file.
    
    - `cc_hex` - attempts to make gaussian rough hexagonal column/plate. Uses method developed by C. Collier, based on Muinonen & Saarinen 2000. If this flag flag is used, several other flags must also be specified:
        - `-cc_hex_l` `<value>` - L from Muinonen & Saarinen 2000. Should be large compared to the correlation length (see below).
        - `-cc_hex_hr` - hexagonal edge length.
        - `-cc_hex_nfhr` - number of subdivisions along each hexagonal edge.
        - `-cc_hex_pfl` - prism edge length.
        - `-cc_hex_nfpl` - number of subdivisions along each prism edge.
        - `-cc_hex_pher` - number of rotations to perform at prism facet-basal facet edges (10% of no. of subfacets along prism edge).
        - `-cc_hex_pper` - number of rotations to perform at prism facet-prism facet edges (10% of subfacets along hexagon edge).
        - `-cc_hex_nscales` - number of roughness scales.
        - `-cc_hex_cls` - correlation lengths for each roughness scale, separated by spaces.
        - `-cc_hex_sds` - standard deviations for each roughness scale, separated by spaces.

- `-rec` `<value>` - Defines the total number of beam recursions per orientation. Is required.

- `-rot` `<args>` - Defines the orientation of the particle. Is optional. If omitted, the ABT will not rotate the input particle. The <args> parameter is used to define the method of rotation, or to define random orientation. Current supported methods are:

    - `euler` <alpha> <beta> <gamma> - Choose to rotate the particle according to the 3 Euler angles, given in degrees. There are several ways to rotate via Euler angles; the ABT follows the method of Mishchenko.

        - example: `-rot euler 11 25 32`

    - `off` `30` `0,10,20,30` - Choose to rotate the particle according to the "off" convention. Only 4 different values are currently supported for this method. Input particle should be oriented lengthways, with prism axis lying in the xy plane.

        - example: `-rot off 30 20`

    - `none` - Choose to not rotate the particle. This behaviour is identical to not including the "rot" argument.

        - example: `rot none`

    - `multi` `<value>` - Choose to randomly orient the particle. <value> defines the number of orientations. For reproducability, the random_seed subroutine can be used to set the random seed to a specified value, which will cause the random orientations to be reproducable.

        - example: `-rot multi 1000`

- `-mt` - Defines whether the code should attempt to use multithreading, where appropriate. Is optional. If omitted, multithreading is disabled by default. If enabled, the user should ensure that the relevant omp environment variables are set up for their system. ie. OMP_STACKSIZE, OMP_NUM_THREADS, etc. If the code throws a segmentation fault at the diffraction subroutine, the user will likely need to increase their OMP_STACKSIZE.

- `-jobname` `<string>` - Specifies the name of the directory within which the output files should be place. Is optional. If omitted, "my_job#" is used, where # is an integer. If the directory already exists, an integer is appended to the directory name so that no files are overwritten.

- `-theta` `values` - Specifies the polar angles at which the far-field should be evaluated. See below for example usage:
    - `-theta 0 1 180` - Evaluate the far-field from 0 in 1 degree steps to 180.
    - `-theta 6 0.1 25 150 175 0.25 180` - Evaluate the far-field from 6 in 0.1 degree steps to 25, then step to 175, then in 0.25 degree steps to 180.

- `-phi` `values` - Specifies the azimuthal angles at which the far-field should be evaluated. See below for example usage:
    - `-phi 0 2 360` - Evaluate the far-field from 0 in 2 degree steps to 360.

- `-no2d` - Suppresses the output of the 2D mueller matrix, which can be a large file if many far-field evaluation points are specified.

- `-tri` - Enables automatic triangulation. Note that use of this flag requires compiling the triangle code in `./src/tri/`.

- `-tri_edge` - Sets the maximum edge length for triangulation.

- `-tri_rough` - Sets the standard deviation for roughness derived from the triangulation.

 ## Examples

 `abt -lambda 0.532 -rbi 1.3117 -ibi 0 -rec 10 -rot euler 20 35 2 -cmethod read -cft obj -cfn my_particle.obj -afn my_apertures.dat -mt -theta 0 1 180 -phi 0 2 360` - run `abt` executable with wavelength 0.532, refractive index 1.3117 + 0i, 10 beam recursions, particle rotated with euler angles 20, 35, read the particle with particle file type `.obj`, filename `my_particle.obj`, apertures specified in the file `my_apertures.dat`, with multithreading enabled. Evaluate the far-field at polar angle 0 in 1 degree steps to 180 and at azithmual angle 0 in 2 degree steps to 360.
