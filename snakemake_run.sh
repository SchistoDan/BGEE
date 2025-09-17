#!/bin/bash
## Conda environment
source /mnt/apps/users/dparsons/conda/etc/profile.d/conda.sh
conda activate bgee_env

# Thread and resource limits
export OMP_NUM_THREADS=1		# Controls how many threads OpenMP (e.g. in NumPy, etc.) can spawn. 1 = single-threaded.
export MKL_NUM_THREADS=1		# Forces Intel’s Math Kernel Library (MKL) to use a single thread.
ulimit -n 8192					# Sets the maximum number of file descriptors (open files, sockets, pipes, etc.) a process can have open simultaneously to handle more concurrent I/O ops.
 
# Setup logging
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_ID="Snakemake BGEE workflow"				
LOG_FILE="snakemake_${TIMESTAMP}.log"

# Function for timestamped logging
log_with_timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Start logging
log_with_timestamp "=== Snakemake Workflow Start ==="
log_with_timestamp "Run ID: $RUN_ID"
log_with_timestamp "Log file: $LOG_FILE"
log_with_timestamp "Screen session: $STY"

# Show system info (ADD RESOURCE INFO)
log_with_timestamp "Running on: $(hostname)"
log_with_timestamp "Working directory: $(pwd)"
log_with_timestamp "Conda environment: $CONDA_DEFAULT_ENV"
log_with_timestamp "File descriptor limit: $(ulimit -n)"
log_with_timestamp "OMP_NUM_THREADS: $OMP_NUM_THREADS"

##Snakemake
# Unlock directory
log_with_timestamp "Unlocking Snakemake directory..."
snakemake --profile ./profiles/slurm/ \
    --snakefile ./workflow/Snakefile-196beta3-5 \
	--configfile ./config/config.yaml \
	--unlock

# Run snakemake workflow with profile
log_with_timestamp "Starting workflow execution..."
snakemake --profile ./profiles/slurm/ \
    --snakefile ./workflow/Snakefile-196beta3-5 \
    --configfile ./config/config.yaml \
    --rerun-incomplete \
    2>&1 | tee -a "$LOG_FILE"

log_with_timestamp "Log file saved as: $LOG_FILE"
