*****************************
*         README            *
* Edelweiss data analysis   *
*  author : Hugues Lattaud  *
*  date   : 02/03/2021      *
*****************************

This framework aims to analyse data comming out of 
NEPAL processing chain. 
It provides :
    - Instance to select events passing desired 
    quality cut and produce skimmed files and histogram.  
    - Instance to Plot various distribution for comparison
    sake.
    -Smooting instance to fit a continuous model onto a 
    discrete distribution. 
    -Smearing instance to apply electron hole pairs statistics,
    selection efficiency and detector resolution smearing.
In the following all line begining by $ is a command to execute.

****************************STEP 1 : Skimming ************************************************************************
To run this step, you need to  create a Txt file listing the name of the run lists you want to run onto.
exemple is provided in List_Run_physics.txt the temperature is also needed (by convention write 20mk this is only needed for labelling) .   
To run the skimmer use the following command line :
Depending if you are running onto CC 
$ source CC_ENV.sh
or Lyoserv cluster
$ source Lyoserv_ENV.sh
$ mkdir List
$ make clean
$ make
$ ./Skimmer --Parity 0 --Temp 20 -o TEST_TO_DELETE -d NbSi209 --inputList Dummy.txt --IonCut 0 --prod prodk 

this command should produce a repository called TEST_TO_DELETE_0.00eV_NbSi209 , this repository should contain one rootfile per run you've skimmed.
These rootfiles contain one tree with all the selected event and the relevant variables (Eh , Ei etc...) and a bunch a histograms tracking the rejected and accepted 
events for the different cuts. The mains histograms (the one containing energy spectrum) are called Ephonon_pos66_* and Ephonon_pos66_keVee* (for spectrum in keVee)

Then  you got the possibility to merge the root file to add up statistic.
The script MERGE_Physicsrun_per_Month.py is here to automatise this procedure, open it and modify the relevant lines.

****************************STEP 2 : Plotting ************************************************************************



