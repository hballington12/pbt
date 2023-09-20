BIN = ./
SRC = ./
SEQ_DIR =./seq/
MPI_DIR =./mpi/

# FSWITCHES = -traceback -debug:full -openmp # debug - windows
# FSWITCHES = -traceback -g -O -qopenmp # debug - linux
# FSWITCHES = -O3 -fopenmp -ffree-line-length-512 -I 'C:\MinGW\lib\gcc\mingw32\6.3.0\include\' # release - windows
# FSWITCHES = -O3 -fopenmp -ffree-line-length-512 -frecursive -fallow-argument-mismatch # release - linux (cluster)
FSWITCHES = -O3 -fopenmp -march=alderlake -ffree-line-length-512 -frecursive -fallow-argument-mismatch # release - linux (Harry's PC)
# FSWITCHES = -O -Wall -Wtabs -g -fopenmp -ffree-line-length-512 -fallow-argument-mismatch # debug - linux

CSWITCHES = -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib

FC = gfortran # fortran compiler
CC = cc # c compiler

# O_EXT = .obj # windows
O_EXT = .o # linux
MOD_EXT = .mod

# objects for sequential code
objects = 	types_mod$(O_EXT) \
			misc_submod$(O_EXT) \
			input_mod$(O_EXT) \
			beam_loop_mod$(O_EXT) \
			diff_mod$(O_EXT) \
			outputs_mod$(O_EXT) \
			cc_hex_mod$(O_EXT)

# all modules
modules = 	types_mod$(MOD_EXT) \
			misc_submod$(MOD_EXT) \
			input_mod$(MOD_EXT) \
			beam_loop_mod$(MOD_EXT) \
			diff_mod$(MOD_EXT) \
			outputs_mod$(MOD_EXT) \
			cc_hex_mod$(MOD_EXT) \
			mpi_mod$(MOD_EXT)

.PHONY: seq # sequential target
.PHONY: mpi # mpi target
.PHONY: triangle # triangle target

seq: $(SEQ_DIR)main.f90 $(objects) $(TRI_BIN)triangle
	$(FC) $(FSWITCHES) $(objects) $(SEQ_DIR)main.f90 -o $(SEQ_DIR)abt

mpi: $(MPI_DIR)main.f90 $(objects) mpi_mod$(O_EXT)
	mpifort $(FSWITCHES) $(objects) mpi_mod$(O_EXT) $(MPI_DIR)main.f90 -o $(MPI_DIR)abt

outputs_mod$(O_EXT) : outputs_mod.f90 misc_submod$(O_EXT)
	$(FC) $(FSWITCHES) -c outputs_mod.f90

diff_mod$(O_EXT) : diff_mod.f90 misc_submod$(O_EXT) types_mod$(O_EXT)
	$(FC) $(FSWITCHES) -c diff_mod.f90

beam_loop_mod$(O_EXT) : beam_loop_mod.f90 input_mod$(O_EXT) misc_submod$(O_EXT) types_mod$(O_EXT)
	$(FC) $(FSWITCHES) -c beam_loop_mod.f90

input_mod$(O_EXT) : input_mod.f90 misc_submod$(O_EXT) types_mod$(O_EXT) cc_hex_mod$(O_EXT)
	$(FC) $(FSWITCHES) -c input_mod.f90

misc_submod$(O_EXT) : misc_submod.f90
	$(FC) $(FSWITCHES) -c misc_submod.f90

types_mod$(O_EXT) : types_mod.f90
	$(FC) $(FSWITCHES) -c types_mod.f90

cc_hex_mod$(O_EXT) : cc_hex_mod.f90
	$(FC) $(FSWITCHES) -c cc_hex_mod.f90

mpi_mod$(O_EXT) : mpi_mod.f90 types_mod$(O_EXT)
	mpifort $(FSWITCHES) -c mpi_mod.f90

./tri/triangle: ./tri/triangle.c
	$(CC) $(CSWITCHES) ./tri/triangle.c -lm -o ./tri/triangle

clean :
	rm $(objects) $(modules) $(SEQ_DIR)abt $(MPI_DIR)abt
