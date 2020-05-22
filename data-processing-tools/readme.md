This folder includes some useful code to process data relevant to the group.

This includes:
(1) bootstrap.py: function to calculate reduced major axis regression slopes and intercepts and the error on these.

(2) tropomi_ut_no2.py: research-grade (as in, still undergoing updates and improvements) to process S5P/TROPOMI NO2 columns to
                       obtain NO2 vertical profiles using the cloud-slicing technique. This is provided to the group as an 
                       example of steps to follow to read and extract data from the S5P/TROPOMI NO2_OFFL NetCDF file.
                       
(3) submit_python: sample code to submit python scripts to a PBS job scheduler and a conda virtual environment. The commands
                   should be the same for a SLURM, except for the queue submission specification statements. 
                   For details of how to output relevant information to a log file (such as tracking the progress of the code),
                   see information included in the tropomi_ut_no2.py file (detailed above). This includes statements to open,
                   close, name, and print to the log file (including flush=True to ensure that the information is immediately
                   written to the log file).                   
