========================================================================
    Fortran Console Application : "abt" Project Overview
========================================================================

This file contains a summary of the aperture beam tracer (ABT), a physical optics hybrid code developed by Harry Ballington and Evelyn Hesse

## abt.f90
Main source file. It contains the program entry point.

## input.txt
Input file. The ABT reads input parameters from this file. There are a few required input parameters, as well as a few optional ones. The code will stop if a mandatory argument is missing from the input file. Each input parameter is defined by a keyphrase, with arguments following a space delimiter. The ordering of lines in the input file is not important and the ABT will search through lines from top to bottom in an attempt to find each argument. If the same keyphrase is specified on multiple lines, the ABT will take the first one that appears. A summary of input file arguments:

- `lambda` `<value>` - Defines the wavelength of incident light. Is required.

- `rbi` `<value>` - Defines the real component of the particle refractive index. Is required.

- `ibi` `<value>` - Defines the imaginary component of the particle refractive index. Is required.

- `cmethod` `<string>` - Defines the method of particle input. Is required. Current supported methods are:

    - `read` - attempts to read the particle from the current directory. If "read" is specified, the following arguments must also be specified:
    
        - `cft` `<string>` - Defines the particle input filetype. The supported particle file input types are:
        
            - `obj` - wavefront style geometry file
        
            - `mrt` - Macke ray-tracing style
        
        - `cfn` `<string>` - Defines the particle filename. Is required. Currently, input file must be a wavefront geometry file. The particle should be triangulated.
    
        - `afn` `<string>` - Defines the apertures filename. Is required. The apertures file containes a single column defining which aperture each face belongs to. The number of lines in the apertures file must match the total number of faces in the particle file.
    
    - `cc_hex` - attempts to make gaussian rough hexagonal column/plate. Input parameters should be included in "vals.in" file placed in the current directory. Uses method developed by C. Collier, based on Muinonen & Saarinen 2000

- `rec` `<value>` - Defines the total number of beam recursions per orientation. Is required.

- `rot` `<args>` - Defines the orientation of the particle. Is optional. If omitted, the ABT will not rotate the input particle. The <args> parameter is used to define the method of rotation, or to define random orientation. Current supported methods are:

    - `euler` <alpha> <beta> <gamma> - Choose to rotate the particle according to the 3 Euler angles, given in degrees. There are several ways to rotate via Euler angles; the ABT follows the method of Mishchenko.

        - example: `rot euler 11 25 32`

    - `off` `30` `0,10,20,30` - Choose to rotate the particle according to the "off" convention. Only 4 different values are currently supported for this method. Input particle should be oriented lengthways, with prism axis lying in the xy plane.

        - example: `rot off 30 20`

    - `none` - Choose to not rotate the particle. This behaviour is identical to not including the "rot" argument.

        - example: `rot none`

    - `multi` `<value>` - Choose to randomly orient the particle. <value> defines the number of orientations. For reproducability, the random_seed subroutine can be used to set the random seed to a specified value, which will cause the random orientations to be reproducable.

        - example: `rot multi 1000`

- `mt` - Defines whether the code should attempt to use multithreading, where appropriate. Is optional. If omitted, multithreading is disabled by default. If enabled, the user should ensure that the relevant omp environment variables are set up for their system. ie. OMP_STACKSIZE, OMP_NUM_THREADS, etc. If the code throws a segmentation fault at the diffraction subroutine, the user will likely need to increase their OMP_STACKSIZE.

- `jobname` `<string>` - Specifies the name of the directory within which the output files should be place. Is optional. If omitted, "my_job#" is used, where # is an integer. If the directory already exists, an integer is appended to the directory name so that no files are overwritten.

/////////////////////////////////////////////////////////////////////////////
Other notes:

/////////////////////////////////////////////////////////////////////////////
