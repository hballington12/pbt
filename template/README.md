# template

This directory contains various template files you may wish to use for running the code.

- `abt.sh`: example shell file for sequential code.

- `abt_mpi.sh`: example shell file for mpi code.

### summary of input flags ###

```
-lambda 0.2 # wavelength
-rbi 1.3118 # real part refractive index
-ibi 0 # imaginary part refractive index
-rec 10 # max number of recursions
-refl 15 # max number of total internal reflections
-rot euler 0 30 30 # euler rotation
-rot off 30 20 # off rotation
-rot none # no rotation
-rot multi 1000 # 1000 random orientations
-jobname a_job # job name
-cmethod read # particle input by read from file
-cft obj # wavefront style input file
-cft mrt # macke style input file
-cfn hollow_column.obj # particle filename
-afn my_apertures.dat # apertures filename
-cmethod cc_hex # particle input by read from file
-cc_hex_l 20 # L from Muinonen & Saarinen 2000
-cc_hex_hr 30 # hex edge length
-cc_hex_nfhr 8 # number of subdivisions along each hexagonal edge
-cc_hex_pfl 3.4 # prism edge length
-cc_hex_nfpl 4 # number of subdivisions along each prism edge (divisible by 4)
-cc_hex_pher 1 # number of rotations to perform at prism facet-basal facet edges
-cc_hex_pper 1 # number of rotations to perform at prism facet-prism facet edges 
-cc_hex_nscales 1 # number of roughness scales
-cc_hex_cls 1 # correlation lengths
-cc_hex_sds 0 # standard deviations
-mt 1 # multithreading (0 to disable)
-theta 0 0.01 0.1 0.1 1 0.5 179 0.1 180 # theta values
-phi 0 4 360 # phi values
-no2d # disable 2d mueller output
-tri # enable automatic triangulation (triangle must be compiled in ./src/tri)
-tri_edge 0.25 # set maximum triangle edge length
-tri_rough 0.05 # set standard devitation for roughness
-time_limit 22 # set time limit to 22 hours
-resume 10 # resume job from directory "10" in cache directory
-scaling # scale beam and ext diffraction to have equal contribution
-timing # enable timing
-debug 3 # set debugging level (0-3)
-export_beam rec 2 # export beams in json format up to and including the 2nd recursion
-export_beam num 50 # export beams in json format up to the 50th beam
-fast_diff # enable faster (but less accurate) diffraction
```