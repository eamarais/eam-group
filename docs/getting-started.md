Resources for students and postdocs using the Marais Research Group resources on the University of Leicester (UoL) servers. 

General Marais Group Resources:
================================
Group website: http://maraisresearchgroup.co.uk/

Group meeting schedule: http://maraisresearchgroup.co.uk/meetings.html

The GitHub /eam-group/doc/ directory inclues .bashrc and .my_personal_settings files for custom settings and aliases for the UoL servers.  

GEOS-Chem Model Resources:
============================
The GEOS-Chem 3D atmospheric chemistry transport model is used extensively in the research group for a range of purposes. Helpful links and resources for beginner, intermediate and advanced used of the model are detailed below.

GEOS-Chem model website: http://acmg.seas.harvard.edu/geos/index.html

GEOS-Chem Youtube channel with helpful tutorials: https://www.youtube.com/channel/UCyh8HWxxxiBCy30xXU8_UDQ

The GitHub /eam-group/doc/ directory includes samples of GEOS-Chem submit and run scripts for the Leicester HPC.
These are specific to version 12.1.0 of the model. Consult the GEOS-Chem user manual to update compile switches to the
version you're using. 

To compile the module, enter "bash script_name" at the command line, where script_name is self-explanatory (name of the compile script). First initiate an interactie nodes (options to do this with an alias in the .my_personal_settings file). 

The sample submit script is for a nested grid simulation over Europe 
(http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids#0.25_x_0.3125_EU_nested_grid). 
This predates the new flexible nested grid feature in GEOS-Chem.

The version of the model you're using may also require more resources than are specified in the submit
script. To determine this, conduct test simulations on the development queue before submitting the job.
(More on queue options here: https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/alice/job-queues/job-queues) 

For simulations other than the nested European grid request the following resources for a one-month simulation as a starting point and run the model in the development queue (see above) to determine whether these are sufficient:<br/>
&nbsp;&nbsp;&nbsp;Global 4x5: 8 cpus, 20 GB, 6 hour walltime<br/>
&nbsp;&nbsp;&nbsp;Global 2x2.5: 8 cpus, 60 GB, 40 hour walltime<br/>
&nbsp;&nbsp;&nbsp;Nested Europe 0.25x0.3125: 16 cpus, 40 GB, 70 hour walltime<br/>
&nbsp;&nbsp;&nbsp;Nested China 0.25x0.3125: 16 cpus, 70 GB, 110 hour walltime<br/>
&nbsp;&nbsp;&nbsp;Nested North America 0.25x0.3125: 16 cpus, 70 GB, 145 hour walltime<br/>

The GEOS-Chem input files are located on the server at /data/uptrop/nobackup/legcfs01/gcgrid/gcdata/ExtData/. These are up-to-date for v12.1.0, so for the version you use may be missing input files (in particular for the HEMCO package). GEOS-Chem now has a neat feature to run the model and print out the files that are needed (YouTube video link: https://www.youtube.com/watch?v=L7T5QtWehLs)

GEOS-Chem input files are centrally located on the Compute Canada servers (http://geoschemdata.computecanada.ca/ExtData/). Send me a list of additional files that need to be added to the server.

University of Leicester IT-Related Resources:
==============================================
Leicester HPC general introduction page: https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc

HPC FAQs: https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/faq

IT isssues (internal only): Submit ticket via https://ithelp.le.ac.uk/<br/>
HPC-related issues: email Research Computing with your query (rcs [dot] support [at] le [dot] ac [dot] uk).

More on how to get IT help (although, some information may not be valid during the lockdown):<br/>
https://uniofleicester.sharepoint.com/sites/staff/get-it-help/SitePages/how-get-it-help.aspx

Online IT, HPC, and programming training slides (replacing in-person training options during covid-19 lockdown):<br/>
https://www2.le.ac.uk/offices/staff-development/events/courses/it/copy_of_it-training-temporary-help

Recommended training:<br/> 
&nbsp;&nbsp;&nbsp;High Performance Computing at Leicester (go through this before starting to use the system)<br/>
&nbsp;&nbsp;&nbsp;Introduction to Version Control<br/>
&nbsp;&nbsp;&nbsp;Linux Introduction (if new to Unix/Linux systems)<br/>
&nbsp;&nbsp;&nbsp;Python Programming at various levels<br/>

HPC service days (mark your calendar with these so you're not caught unawares):<br/> https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/service-days

Weekly software surgery is on Thursdays at 4-5pm (via Microsoft Teams during lockdown). Delivered by a UoL python expert.<br/> 
Either join to get help with coding or to help others learning to code. Email reminder and Teams link sent weekly to the HPC email list by a Research Software Engineering Specialist. You should be added to the HPC email list once you register to access the UoL servers.

NoMachine virtual desktop to access Leicester services:<br/>
https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/spectre/access/nx5<br/>
https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/spectre/access/nx-basic-usage<br/>

ZendTo file share tool for sharing large files with people at and outside UoL: https://filedrop.le.ac.uk/

Other UoL Resources:
=====================
College external and inaugural seminar series are posted here: 
https://www2.le.ac.uk/colleges/scieng/internal/college-research-seminar-series-2014-15/college-research-seminar-series-2016-17 (link is correct; ignore the dates)

You can also keep pace with college news in the college newsletter, E-ZINE: 
https://www2.le.ac.uk/colleges/scieng/internal/e-zine-folder

EOS Group Resources and Information:
======================================
There is an EOS email list. Check with Bev Robinson that you have been added to this to get updates on group meetings and
any details on group activities.

EOS wiki (internal) with EOS group seminar series schedule: https://wiki.lamp.le.ac.uk/eos/index.php/Seminars<br/>
If you can't access the wiki with your UoL login, ask Robert Parker to give you access.

The wiki also includes IT and Research Computing guidelines: <br/>
https://wiki.lamp.le.ac.uk/eos/index.php/IT<br/>
https://wiki.lamp.le.ac.uk/eos/index.php/Computing<br/>

The EOS seminar organizer is Neil Humpage. If you'd like to invite someone to present at the seminar series (virtually
during the lockdown), contact Neil with details of the person you'd like to invite.

EOS group meetings occur monthly or biweekly depending on the time of year. These can also occur with the Centre 
for Landscape and Climate Research (CLCR; https://www2.le.ac.uk/colleges/scieng/research/centres/clcr). EOS only and 
EOS-CLCR seminars are on Fridays at 10-11:30am and typically include 2 speakers. Kamil Mroz, postdoc in EOS, organizes 
these and sends calendar invites of upcoming talks to the EOS mailing list. 

Duncan Ross is Health and Safetly Officer. After the lockdown, when we're permitted to be back on campus, schedule a 
health and safety training session with him. This includes rules to follow when accessing the building after hours. 
