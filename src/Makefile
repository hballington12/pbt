BIN = ./
SRC = ./
SEQ_DIR =./seq/
MPI_DIR =./seq/

debug: $(SEQ_DIR)abt.f90
	ifort types_mod.f90 misc_submod.f90 input_mod.f90 beam_loop_submod.f90 diff_mod.f90 outputs_mod.f90 $(SEQ_DIR)abt.f90 -g -traceback -O -qopenmp -o $(SEQ_DIR)abt.o

debug_mpi: $(MPI_DIR)abt_mpi.f90
	mpif90 -f90=ifort types_mod.f90 misc_submod.f90 input_mod.f90 beam_loop_submod.f90 diff_mod.f90 outputs_mod.f90 $(MPI_DIR)abt_mpi.f90 -g -traceback -O -qopenmp -o $(MPI_DIR)abt.o

release: $(SEQ_DIR)abt.f90
	ifort types_mod.f90 misc_submod.f90 input_mod.f90 beam_loop_submod.f90 diff_mod.f90 outputs_mod.f90 $(SEQ_DIR)abt.f90 -O3 -qopenmp -xAVX -o $(SEQ_DIR)abt.o

release_mpi: $(MPI_DIR)abt_mpi.f90
	mpif90 -f90=ifort types_mod.f90 misc_submod.f90 input_mod.f90 beam_loop_submod.f90 diff_mod.f90 outputs_mod.f90 $(MPI_DIR)abt_mpi.f90 -O3 -qopenmp -xAVX -o $(MPI_DIR)abt.o

release_debug: abt.exe
	ifort types_mod.f90 misc_submod.f90 input_mod.f90 beam_loop_submod.f90 diff_mod.f90 outputs_mod.f90 abt.f90 -O3 -debug:full -Qopt-report=5 -openmp -o abt.exe