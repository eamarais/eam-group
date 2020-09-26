#!/bin/bash

#################### Begin SLURM header ##########################
# Settings for the SLURM scheduler
#SBATCH -t 3:00:00 
#SBATCH --job-name TROPOMI_NO2_oversampling_1_preparation_job_1
#SBATCH --qos bblargemem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=300G
####################  End SLURM header  ##########################

# load Python packages needed
module purge; module load bluebear
module load IPython/7.9.0-foss-2019b-Python-3.7.4
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load netcdf4-python/1.5.3-foss-2019b-Python-3.7.4

######################################################################################################################
# remove the "#" to choose the job, assign a unique job name and submit it, then repeat the process

#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20190801 --end_date 20190831
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20190901 --end_date 20190930
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20191001 --end_date 20191031
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20191101 --end_date 20191130
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20191201 --end_date 20191231
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200101 --end_date 20200131
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200201 --end_date 20200229
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200301 --end_date 20200331
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200401 --end_date 20200430
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200501 --end_date 20200531
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200601 --end_date 20200630
#python TROPOMI_NO2_oversampling_1_preparation.py --domain AF --qa_flag 0.75 --start_date 20200701 --end_date 20200731

#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20190801 --end_date 20190831
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20190901 --end_date 20190930
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20191001 --end_date 20191031
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20191101 --end_date 20191130
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20191201 --end_date 20191231
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200101 --end_date 20200131
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200201 --end_date 20200229
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200301 --end_date 20200331
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200401 --end_date 20200430
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200501 --end_date 20200531
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200601 --end_date 20200630
#python TROPOMI_NO2_oversampling_1_preparation.py --domain EU --qa_flag 0.75 --start_date 20200701 --end_date 20200731

######################################################################################################################

# Exit normally
exit 0
#EOC
