#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime

DirList = [0.4]#-1000.,0.,0.2]


for iDir in DirList:
	
	cmd = "./Skimmer --Parity 2 --Temp 20 -o TEST_TO_DELETE -d NbSi209 --inputList List_Run_physics.txt --IonCut "+str(iDir)+" --prod prodj"#+" --RejectedIoncut"
	print(cmd)
	os.system(cmd)

