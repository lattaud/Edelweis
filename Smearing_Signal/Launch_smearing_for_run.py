#!/usr/bin/python
#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen


List_Mass = ["15","20","30","40","50","100","500","1000"]
List_Eicut = ["0.00","100.00","200.00","300.00","400.00"]

for mass in List_Mass:
        for Cut in List_Eicut:
            cmd = "./Smearing --Resolution_phonon 0.1 --Fano 0.15 --v 0 --Percent 1 --Voltage 66 --efficiency-name tk20a002-tk22a000-tk26a000-tk19a000-tk24a001 -o Smeared_spectrum_impairpartion_physics_"+Cut+"eicut.root -d NbSi209 --inputplot-list List_mass.txt --inputfile Spectrum_migdal_Eer_updated.root  --Mass "+mass+" --Eicut "+Cut
            print(cmd)
            os.system(cmd)
