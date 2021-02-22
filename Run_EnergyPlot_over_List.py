#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime

DirList = [0.4]#-1000.,0.,0.2]


for iDir in DirList:
	
	cmd = "./Skimmer --Parity 2 --Temp 20 -o PHYSIC_RUN_Pair_Partition_Compton_study -d NbSi209 --inputList List_Run_physics.txt --IonCut "+str(iDir)#+" --RejectedIoncut"
	print(cmd)
	os.system(cmd)

