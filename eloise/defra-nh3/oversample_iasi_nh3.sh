#!/bin/bash
#
#PBS -N iasi_nh3_oversample
#PBS -l nodes=1:ppn=1
#PBS -l walltime=05:00:00
#PBS -l vmem=4gb

# Load environment modules
export OMP_NUM_THREADS=1

# Seasons:
#seas=("jja" "son" "djf" "mam")

#-----------------------------------------------------------------
# Python processing of IASI NH3:
#-----------------------------------------------------------------
# Move to relevant working directory that has the Python script:
cd /data/uptrop/Projects/DEFRA-NH3/python/

# Activate your virtual environment:
source /home/e/em440/miniconda/envs/gc_uptrop/bin/activate

for SS in `seq 1 12`
    do

    MM=$SS

    if [ $MM -le 9 ]; then month="0"$MM; fi
    if [ $MM -gt 9 ]; then month=$MM;    fi

    # Define log file:
    log_py='log_prepare_iasi_nh3_'$month
    # Run the code:
    python process_iasi_nh3.py --month=$month > $log_py   

done

# Deactivate the virtual environment:
conda deactivate

# Resolution options:
resval=("0.5" "0.25" "0.2" "0.15" "0.1" "0.05")

for R in `seq 1 6`
    do

    echo ${resval[R-1]}

    for SS in `seq 1 12`
        do 

	MM=$SS

	if [ $MM -le 9 ]; then month="0"$MM; fi
	if [ $MM -gt 9 ]; then month=$MM;    fi

	#-----------------------------------------------------------------
	# Run Fortran oversampling code:
	#-----------------------------------------------------------------
	# Move to relevant working directory that has the Python script:
	cd /data/uptrop/Projects/DEFRA-NH3/Fortran/RegridPixels/
    
	# Compile the fortran code to make sure it's up to date:
	#gfortran -o RegridPixels.x cakecut_m.f90 tools_m.f90 RegridPixels.f90

	# Set paths:
	Input_Dir='/data/uptrop/Projects/DEFRA-NH3/Data/IASI_Oversampled/'
	Output_Dir='/data/uptrop/Projects/DEFRA-NH3/Data/IASI_Oversampled/'

	# Set output resolution:
	Res=${resval[R-1]}

	# Set file names
	Input_Filename='iasi_nh3_'$month'_2008-2018'
	Output_Filename='iasi_nh3_oversampled_'$month'_2008-2018_'$Res

	# Define fortran log file:
	log_fort='log_oversample_iasi_nh3_'$month'_'$Res

	# run fortran code:
	./RegridPixels.x<<EOF > $log_fort    
$Input_Dir
$Output_Dir 
$Input_Filename
$Output_Filename
$Res
EOF

done

done

exit
