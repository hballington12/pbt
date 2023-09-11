BIN = ./
SRC = ./
SEQ_DIR =./seq/
MPI_DIR =./mpi/

# FSWITCHES = -traceback -debug:full -openmp # debug - windows
# FSWITCHES = -traceback -g -O -qopenmp # debug - linux
# FSWITCHES = -O3 -openmp # release - windows
FSWITCHES = -O3 -qopenmp -xAVX # release - linux

# O_EXT = .obj # windows
O_EXT = .o # linux
MOD_EXT = .mod

objects = 	types_mod$(O_EXT) \
			misc_submod$(O_EXT) \
			input_mod$(O_EXT) \
			beam_loop_mod$(O_EXT) \
			diff_mod$(O_EXT) \
			outputs_mod$(O_EXT) \
			cc_hex_mod$(O_EXT)

modules = 	types_mod$(MOD_EXT) \
			misc_submod$(MOD_EXT) \
			input_mod$(MOD_EXT) \
			beam_loop_mod$(MOD_EXT) \
			diff_mod$(MOD_EXT) \
			outputs_mod$(MOD_EXT) \
			cc_hex_mod$(MOD_EXT)

.PHONY: abt # sequential target
.PHONY: abt_mpi # mpi target

abt: $(SEQ_DIR)main.f90 $(objects)
	ifort $(FSWITCHES) $(objects) $(SEQ_DIR)main.f90 -o $(SEQ_DIR)abt

abt_mpi: $(MPI_DIR)main.f90 $(objects)
	mpif90 -f90=ifort $(FSWITCHES) $(objects) $(MPI_DIR)main.f90 -o $(MPI_DIR)abt

outputs_mod$(O_EXT) : outputs_mod.f90 misc_submod$(O_EXT)
	ifort $(FSWITCHES) -c outputs_mod.f90

diff_mod$(O_EXT) : diff_mod.f90 misc_submod$(O_EXT) types_mod$(O_EXT)
	ifort $(FSWITCHES) -c diff_mod.f90

beam_loop_mod$(O_EXT) : beam_loop_mod.f90 input_mod$(O_EXT) misc_submod$(O_EXT) types_mod$(O_EXT)
	ifort $(FSWITCHES) -c beam_loop_mod.f90

input_mod$(O_EXT) : input_mod.f90 misc_submod$(O_EXT) types_mod$(O_EXT) cc_hex_mod$(O_EXT)
	ifort $(FSWITCHES) -c input_mod.f90

misc_submod$(O_EXT) : misc_submod.f90
	ifort $(FSWITCHES) -c misc_submod.f90

types_mod$(O_EXT) : types_mod.f90
	ifort $(FSWITCHES) -c types_mod.f90

cc_hex_mod$(O_EXT) : cc_hex_mod.f90
	ifort $(FSWITCHES) -c cc_hex_mod.f90

clean :
	rm $(objects) $(modules) $(SEQ_DIR)abt

clean_mpi :
	rm $(objects) $(modules) $(MPI_DIR)abt