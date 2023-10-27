#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l nodes=32 -l pmem=1gb
#PBS -k oe
#PBS -o /adda/clusterTest/output.txt
#PBS -e /adda/clusterTest/error.txt
#PBS -N adda
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

# usage: qsub -t 1-N%M adda_loop.sh
# where:
# N is the total number of orientations to be read from your file (see below)
# M is the maximum number of jobs running at once
# note that EACH job will request to use the number of nodes specified above, so use M
#   to ensure that you don't ask for too many resources at once 

# copy-paste from chatgpt: https://chat.openai.com/c/e10ce1bb-6f98-4a5e-87e1-384d7ea70644

# This Bash script defines a function `read_file_into_arrays` that reads values from a file with three columns and stores these values in three arrays. Here's a step-by-step explanation of what's happening in the code:

# 1. `#!/bin/bash`: This is called a shebang, and it specifies that this script should be interpreted using the Bash shell.

# 2. `read_file_into_arrays() { ... }`: This defines a Bash function named `read_file_into_arrays`. The function is responsible for reading data from a file and storing it in three arrays: `column1`, `column2`, and `column3`.

# 3. Inside the `read_file_into_arrays` function:
#    - `local file="$1"`: This line stores the first argument passed to the function (the filename) in a local variable named `file`.
#    - `local column1=()`, `local column2=()`, and `local column3=()`: These lines declare three arrays, `column1`, `column2`, and `column3`, to store the values from the three columns.

#    - `while read -r col1 col2 col3; do ... done`: This `while` loop reads lines from the file and separates them into three variables, `col1`, `col2`, and `col3`, using white spaces as delimiters. This loop continues until there are no more lines to read.
#      - `column1+=("$col1")`, `column2+=("$col2")`, and `column3+=("$col3")`: These lines append the values of `col1`, `col2`, and `col3` to their respective arrays (`column1`, `column2`, and `column3`).

#    - `done < <(cat "$file" ; echo)`: This construct uses process substitution to read the contents of the file specified by the `file` variable. The `; echo` is used to ensure that the last line (which might not have a newline character) is properly read.

#    - The script concludes the function by echoing the values in each of the three arrays.

# 4. Following the function definition, there's an if statement:
#    - `if [ $# -eq 0 ]; then`: This checks if the script was called without any command-line arguments (i.e., the count of arguments, `$#`, is zero).
#      - If there are no arguments, it prints a usage message and exits the script with an error code.

# 5. Finally, the script calls the `read_file_into_arrays` function with the first argument passed to the script (the filename). This will read the file and store the data in the arrays, and the values in the arrays are echoed to the terminal.

# So, this script can be executed with a filename as an argument, and it will read the file, store the data in arrays, and display the values in those arrays.

# FUNCTIONS ======================

# Function to read values from a file into arrays
read_file_into_arrays() {
  local file="$1"
  declare -g alphas
  declare -g betas
  declare -g gammas

  while read -r col1 col2 col3; do
    alphas+=("$col1")
    betas+=("$col2")
    gammas+=("$col3")
  done < <(cat "$file" ; echo)

  # Output the arrays
  # echo "Read alpha values: ${alphas[*]}"
  # echo "Read beta values: ${betas[*]}"
  # echo "Read gamma values: ${gammas[*]}"
}

# SCRIPT ======================

# set/load environment
ulimit -s unlimited
export OMP_NUM_THREADS=32
eval `/usr/bin/modulecmd bash load openmpi-4.0.5`
cd /home/hballington/adda/clusterTest

# filename for reading angles from
my_file="96eulers.dat"

# Call the function with the provided file
read_file_into_arrays "$my_file"

# use the PBS_ARRAYID variable to give each array job a different set of angles
my_alpha="${alphas[$(($PBS_ARRAYID - 1))]}"
my_beta="${betas[$(($PBS_ARRAYID - 1))]}"
my_gamma="${gammas[$(($PBS_ARRAYID - 1))]}"

echo "my alpha: $my_alpha"
echo "my beta: $my_beta"
echo "my gamma: $my_gamma"

# set the crystal filename directory
cfn="hr5_nfhr16_pfl10_nfpl32_cl0.5_sd0.1_size66.55"

# execute adda
mpiexec /home/hballington/adda/src/mpi/adda_mpi \
-orient $my_alpha $my_beta $my_gamma \
-dir /beegfs/cair/hballington/adda/reference_data/${cfn}/a${my_alpha}_b${my_beta}_g${my_gamma} \
-store_scat_grid \
-shape read hr5_nfhr16_pfl10_nfpl32_cl0.5_sd0.1_size66.55.dat \
-lambda 0.532 \
-dpl 13.103448 \
-m 1.31 0 \
-eps 5 \
-iter bcgs2 \
-prop 0 0 -1
