#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime

DirList   = ["0.00","0.20"]#"-1000.00",
EList   = ["0.00","200.00"]#"-1000.00",
MonthList = ["tc17a000-tc19a000","td02a000-td03a000-td04a000-td17a000-td18a000-td22a000","tj27a001-tj28a000-tj29a000","tk19a000-tk20a002-tk22a000-tk24a001-tk26a000","tl15a001-tl15a002-tl18a000-tl22a000-tl24a000","ud21a000-ud22a000-ud23a000-ud24a000-ud27a001-ud29a000","ue25a001-ue27a000-ue28a000-ue29a000-ue30a000","uf01a000-uf02a000-uf05a001-uf07a000-uf08a000-uf09a000"]
RootList = ["Mars2019.root","Avril2019.root","Octobre2019.root","Novembre2019.root","Decembre2019.root","Avril2020.root","Mai2020.root","Juin2020.root"]

it1 = 0 
for iDir in DirList:
    it = 0
    for month in MonthList:	
        cmd = "./Plot -i ../PHYSIC_RUN_Pair_Partition_Compton_study_"+iDir+"eV_NbSi209/"+RootList[it]+" --spectrum --tension Tension_to_plot.txt -t 20 -o "+iDir+"evcut"+month+"_compton_study -d NbSi209 --Ion "+EList[it1]+" --eff "+month
        it = it + 1
        print(cmd)
        os.system(cmd)
    it1 = it1 + 1

