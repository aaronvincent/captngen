#!/bin/bash
#SBATCH -c 64		# Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=128G	# Memory requested in megabytes. If omitted, the default is 1024 MB. [K|M|G|T]
#SBATCH -t 0-1:0:0	# d-h:m:s  How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J capOMP	# Job name
#SBATCH -p reserved     # Partition (reserved or standard), assumes standard
#SBATCH --mail-user=neal.aviskozar@queensu.ca
#SBATCH --mail-type=ALL	#BEGIN, END, FAIL, ALL
date # Prints the current date and time
hostname # Prints the current node

echo 'starting test job...'

make clean
make
time gentest.x

echo 'our job worked!'

date # Prints the current date and time
