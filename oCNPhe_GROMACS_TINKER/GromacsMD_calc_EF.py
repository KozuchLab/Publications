
# The purpose of this script is to calculate electric fields for a solute in solvent with carbonyl vibrational probe. 
# usage: python calc_EF_GROMACS.py -a 1365  1366 -at C N -ch 0.2278 -0.3538 -t PYP_md

import sys, argparse, glob
import re
import itertools
import os
import numpy as np
from subprocess import call,Popen

### Argparser ###
def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 2)",nargs='+',type=int,required=True)
   parser.add_argument("-at","--atomtypes",help="Atomtypes for analysis (Max 2)",nargs='+',type=str,default=['C','O'])
   parser.add_argument("-ch","--charge",help="Charges on atoms 1 & 2. ",nargs=2,type=float,required=True)
   parser.add_argument("-t","--traj",help="Filename (without suffix .trr) of solvated trajectory",type=str,required=True)
   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args()
   return args

### Check if it's number ###
def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

### Write into log file ###
def log_write(string):
    logf = 'calc_EF_gromacs.log'
    f = open(logf, 'a')
    f.write(string+'\n')
    f.close
    print(string+'\n')

### Obtain coordinates and forces ###
def get_coords_forces():
    log_write("Reading trajectory %s.trr and the corresponding discharged solvent trajectory (*_0q.trr)." % (args.traj))
    trrfile = args.traj + '.trr'
    trrfile_0q = args.traj + '_0q.trr'
    tprfile = args.traj + '.tpr'
    tprfile_0q = args.traj + '_0q.tpr'

    f = open('probe.ndx', 'w')
    f.write('[ carbonyl ] \n \n   %d    %d  \n' %(args.atoms[0],args.atoms[1]))
    f.close
   
    f = open('calc_EF_gromacs.log','a')
    p=Popen('gmx traj -s ' + tprfile + ' -f  ' + trrfile + ' -n probe.ndx -ox co_coords.xvg', shell=True, stdout=f, stderr=f)
    p.wait()
    f.close
    f = open('calc_EF_gromacs.log','a')
    p=Popen('gmx traj -s ' + tprfile + ' -f  ' + trrfile + ' -n probe.ndx -ox co_forces.xvg', shell=True, stdout=f, stderr=f)
    p.wait()
    f.close
    f = open('calc_EF_gromacs.log','a')
    p=Popen('gmx traj -s ' + tprfile_0q + ' -f  ' + trrfile_0q + ' -n probe.ndx -ox co_forces_0q.xvg', shell=True, stdout=f, stderr=f)
    p.wait()
    f.close


#    call('gmx traj -s %s -f %s -n probe.ndx -ox co_coords.xvg' % (tprfile,trrfile), shell=True, stdout=None, stderr=None)
#    call('gmx traj -s %s -f %s -n probe.ndx -of co_forces.xvg' % (tprfile,trrfile), shell=True, stdout=None, stderr=None)
#    call('gmx traj -s %s -f %s -n probe.ndx -of co_forces_0q.xvg' % (tprfile_0q,trrfile_0q), shell=True, stdout=None, stderr=None)

### Main ###
global args
args = parse()
log_write('######################################################################')
log_write('# Script to calculate GROMACS-based electric field at a desired site #')
log_write('######################################################################')

#f = open('probe.ndx', 'w')
#f.write('[ carbonyl ] \n \n   %d    %d  \n' %(args.atoms[0],args.atoms[1]))
#f.close

#get_coords_forces()


# Get the relevant forces and coords from trajectories
f_x = open('co_coords.xvg' , 'r')
f_f = open('co_forces.xvg' , 'r')
f_f0q = open ('co_forces_0q.xvg'  , 'r')


# Calculate the the fields
fields = open('Electric_Fields_on_%s%s.txt' %(args.atomtypes[0],args.atomtypes[1]), 'w')
for line_x, line_f, line_f0q in zip(f_x, f_f, f_f0q):
        if isnumber( line_x.split()[0] ):
                #Get the time
                time_x = float( line_x.split()[0] )
                time_f0q = float( line_f0q.split()[0] )

                #Check the info
                if time_x != time_f0q:
                        print("time indeces do not match. Error!")
                        sys.exit()

                #Get the C and O coordinates
                #Calculate the CO bond length and CO unit vector
                x_C = np.array([float(x) for x in line_x.split()[1:4]])
                x_O = np.array([float(x) for x in line_x.split()[4:7]])
                COvec = x_O - x_C
                COlen = np.sqrt((COvec**2).sum())
                COunitvec = COvec/ COlen

     #Get the C and O forces
                f_C = np.array([float(x) for x in line_f.split()[1:4]])
                f_O = np.array([float(x) for x in line_f.split()[4:7]])

                #Get the C and O discharged forces
                f0q_C = np.array([float(x) for x in line_f0q.split()[1:4]])
                f0q_O = np.array([float(x) for x in line_f0q.split()[4:7]])

                #Calculate the electrostatic forces
                fe_C = f_C - f0q_C
                fe_O = f_O - f0q_O

                #Calculate the electrostatic force projection onto CO
                feproj_C = np.dot( fe_C, COunitvec )
                feproj_O = np.dot( fe_O, COunitvec )
                feproj_C_0q = np.dot( f0q_C, COunitvec )
                feproj_O_0q = np.dot( f0q_O, COunitvec )

                #Calculate the electric field and convert
                Eproj_C = (feproj_C /args.charge[0]) * 0.1036427
                Eproj_O = (feproj_O /args.charge[1]) * 0.1036427
                Eproj = (Eproj_C + Eproj_O)/2
                EprojDrop = (Eproj_O-Eproj_C)
                Eproj_C_0q = (feproj_C_0q /args.charge[0]) * 0.1036427
                Eproj_O_0q = (feproj_O_0q /args.charge[1]) * 0.1036427
                Eproj_0q = (Eproj_C_0q + Eproj_O_0q)/2
                EprojDrop_0q = (Eproj_O_0q - Eproj_C_0q)

                #Print to file
                fieldInfo = str( time_x ) + '\t' +str( Eproj_C ) + '\t' + str( Eproj_O ) + '\t' + str( Eproj )+ '\t' + str( EprojDrop ) + '\t' + str( COlen ) + '\t' + str( Eproj_0q )+ '\t' + str( EprojDrop_0q ) + '\n'
                fields.write ( fieldInfo )
f_x.close()
f_f.close()
f_f0q.close()
fields.close()

#Determine and write averages
fields = np.loadtxt(open('Electric_Fields_on_%s%s.txt' %(args.atomtypes[0],args.atomtypes[1])  , 'r'))
field_avg = np.mean(fields[:,3])
field_std = np.std(fields[:,3])
fielddrop_avg = np.mean(fields[:,4])
fielddrop_std = np.std(fields[:,4])
log_write(('Average electric field on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],field_avg,field_std))
log_write(('Average electric field drop on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],fielddrop_avg,fielddrop_std))
    
selffield_avg = np.mean(fields[:,6])
selffield_std = np.std(fields[:,6])
selffielddrop_avg = np.mean(fields[:,7])
selffielddrop_std = np.std(fields[:,7])
log_write(('Average self-electric field on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],selffield_avg,selffield_std))
log_write(('Average self-electric field drop on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],selffielddrop_avg,selffielddrop_std))

log_write('Electric Fields determined successfully.')
#########################

