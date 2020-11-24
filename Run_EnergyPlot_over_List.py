#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime

DirList = [0.,0.1,0.2,0.3,0.4]


for iDir in DirList:
	
	cmd = "./Skimmer --Parity 2 --Temp 20 -o PHYSIC_RUN_Pair_Partition_Ei_FID -d NbSi209 --inputList List_Run_physics.txt --IonCut "+str(iDir)
	print(cmd)
	os.system(cmd)

