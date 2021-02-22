#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime


DirList = ["0.40"]#"-1000.00","-1000.00","0.00","0.20"
Month_list =["Mars2019","Avril2019","Octobre2019","Novembre2019","Avril2020","Mai2020","Juin2020","ALL_physics"]
Month_Acro_list =["tc","td","tj","tk","ud","ue","uf","*"]


for iDir in DirList:
    iterator = 0 
    for month in Month_list:
        cmd="hadd -f PHYSIC_RUN_Pair_Partition_Compton_study_"+iDir+"eV_NbSi209/"+month+".root PHYSIC_RUN_Pair_Partition_Compton_study_"+str(iDir)+"eV_NbSi209/Eh_perrun_22.000000mk_"+Month_Acro_list[iterator]+"*" 
        print(cmd)
        os.system(cmd)
        iterator = iterator + 1
