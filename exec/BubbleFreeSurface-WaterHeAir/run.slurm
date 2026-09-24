#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Example Slurm job script. Replace the <...> placeholders for your machine.
#----------------------------------------------------

#SBATCH -J BFS_6EQ_NP                   # Job name
#SBATCH -o out.%j                       # Name of stdout output file
#SBATCH -e err.%j                       # Name of stderr error file
#SBATCH -p <partition>                  # Queue (partition) name
#SBATCH --nodes=2                       # Total number of nodes
#SBATCH --ntasks=256                    # Total number of MPI ranks
#SBATCH -t 6:00:00                      # Run time (hh:mm:ss)
#SBATCH -A <allocation>                 # Allocation name
##SBATCH --mail-user=<your-email>       # Uncomment to receive job emails
##SBATCH --mail-type=all
###### End Slurm Preamble

EXEC=$(ls ./main*.ex | head -n 1)
INPUTS=./prob/inputs

srun -n "${SLURM_NTASKS}" "${EXEC}" "${INPUTS}"
