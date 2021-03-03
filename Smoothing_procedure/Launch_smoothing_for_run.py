#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen


List_Folder = ["PHYSIC_RUN_Pair_Partition_Compton_study_0.00eV_NbSi209","PHYSIC_RUN_Pair_Partition_Compton_study_0.20eV_NbSi209","PHYSIC_RUN_Pair_Partition_Compton_study_-1000.00eV_NbSi209","PHYSIC_RUN_Pair_Partition_Compton_study_0.40eV_NbSi209"]#,"PHYSIC_RUN_Pair_Partition_Compton_study_0.00eV_NbSi209","PHYSIC_RUN_Pair_Partition_Compton_study_0.20eV_NbSi209","PHYSIC_RUN_Pair_Partition_Compton_study_-1000.00eV_NbSi209"]
List_Month = ["ALL_physics"]#"Mars2019","Avril2019","Octobre2019","Novembre2019","Avril2020","Mai2020","Juin2020"]#"Avril2019_acti"]#
for Folder in List_Folder:
        for month in List_Month:
            cmd = "./Smooth_spectrum --tree --Voltage 66 --Run-name "+Folder+"_"+month+" --output-name "+Folder+"_"+month+" --inputplot-list list_hist --inputfile ../"+Folder+"/"+month+".root -d NbSi209"
            print(cmd)
            os.system(cmd)


