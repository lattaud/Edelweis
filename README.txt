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

****************************STEP 1 : Skimming ******************************
To run this step, you need to  create a Txt file listing the name of the run lists you want to run onto.
exemple is provided in List_Run_physics.txt the temperature is also needed (by convention write 20mk this is only needed for labelling) .   
To run the skimmer use the following command line :
$ mkdir List



