###########################################################################################################
# This is a python-based analysis tool to evaluate structures, RMSFs, distances of Protein-Ligand complex #
# Written by Jacek Kozuch (Version 10, May 2020)                                                          #
###########################################################################################################
# Before starting do:                                                                                     
# gmx trjconv -s *.tpr -f *.trr -o *_center.xtc -center -pbc mol -ur compact                              
# gmx trjconv -s *.tpr -f *_center.xtc -o *_fit.xtc -fit rot+trans    
# Then do gmx hbond on your system:
# gmx make_ndx -f *.tpr -o z_index_hbond.ndx 
# 
# with group of acceptor and group of all without acceptor
# Then send with enough RAM (on server):
# gmx hbond -s .tpr -f *_fit.xtc -n z_index_hbond.ndx -num zz_Hbond_num.xvg -dist zz_Hbond_dist.xvg -hbn zz_Hbond_ndx.ndx -r 0.35 -shell 0.35 -a 30 << EOF
# 20
# 19
# 20
# EOF
###########################################################################################################

######################################################
# Needed Modules:                                    #
# module load chemistry gromacs py-numpy/1.14.3_py27 #
######################################################

import os
import numpy as np
import re
from subprocess import call,Popen,PIPE
import sys, argparse, glob
import itertools

#########################
######## Parse ##########
#########################

# example
# python z_gmx_hbond_extended_v3.py -p TLN_wat_md -d 0.35 -a 30 -accatm 12 -accref 11
# python z_gmx_hbond_extended_v3.py -p PYP_md -d 0.35 -a 30 -accref 406 -accatm 407 
# python z_gmx_hbond_extended_v3.py -p PYP_md -d 0.35 -a 30 -accatm 903 -accref 902
# python z_gmx_hbond_extended_v3.py -p PYP_md -d 0.35 -a 30 -accatm 1366 -accref 1365
# python z_gmx_hbond_extended_v3.py -p PYP_md -d 0.35 -a 30 -accatm 1436 -accref 1435


def parse():
   parser = argparse.ArgumentParser()
   parser.add_argument("-p","--prefix",     help="Filename of *.tpr and *.trr files of your trajectory",type=str,required=True)
   parser.add_argument("-d","--dist",       help="Max. H-bond D - A distance.",type=float,required=True)
   parser.add_argument("-a","--angle",      help="Max. H-bond H - D - A angle.",type=float,required=True)
   parser.add_argument("-accatm","--accatm",    help="For instance, for CN group: atom number of N atom.",type=int,required=True)   
   parser.add_argument("-accref","--accref",    help="For instance, for CN group: atom number of C atom.",type=int,required=True)
   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)
   args = parser.parse_args()
   return args


############################
### Check if it's number ###
############################
def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


#########################
### Create Index File ###
#########################

### Check if numbers in line and not header ###
def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

##########################
### Generate ndx Files ###
##########################

    # do all separately for Water and AAs at some point...
    # create index.ndx file with separate groups so that they can be used in gmx distance/pairdist and gmx angle
    # determine A - D distance for all pairs along entire trajectory
    # determine A - H - D angle ...
    # compare to num
    # separate non-hbonded, 1 hbonded, 2, 3... get for each number of frames/probability, ave/std distance and angles


def gen_hbond_ndx():
    check = 0
    line_test = 0
    HB_ndx_mat = []
    with open('zz_Hbond_ndx.ndx','r') as file_HB_ndx:
        for line in file_HB_ndx:
            line_test = line_test + 1
            if (check == 1 and '[ ' in line):
                check = 0
            if '[ donors_hydrogens' in line:
                check = 1
                continue
            if check == 1:
                HB_ndx = np.array([int(x) for x in line.split()[:]])
                if len(HB_ndx) == 2:
                    if HB_ndx_mat == []:
                        HB_ndx_mat = HB_ndx
                    else:
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx))
                if len(HB_ndx) == 4:
                    if HB_ndx_mat == []:
                        HB_ndx_mat = HB_ndx[0:2]
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                    else:
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[0:2]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                if len(HB_ndx) == 6:
                    if HB_ndx_mat == []:
                        HB_ndx_mat = HB_ndx[0:2]
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[4:6]))
                    else:
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[0:2]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[4:6]))
                if len(HB_ndx) == 8:
                    if HB_ndx_mat == []:
                        HB_ndx_mat = HB_ndx[0:2]
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[4:6]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[6:8]))
                    else:
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[0:2]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[2:4]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[4:6]))
                        HB_ndx_mat = np.vstack((HB_ndx_mat, HB_ndx[6:8]))
#    print(HB_ndx_mat)
                
    don_atom = np.transpose(HB_ndx_mat)[0]
    hyd_atom = np.transpose(HB_ndx_mat)[1]
    acc_atom = np.transpose(np.full((len(don_atom)),'%d' %(args.accatm), dtype = int))
    ref_atom = np.transpose(np.full((len(don_atom)),'%d' %(args.accref), dtype = int))

    # Create new ndx files:
    count = 0
    howmany = 0
    with open('zz_index_AD_dist.ndx','w') as index_AD_dist, open('zz_index_AH_dist.ndx','w') as index_AH_dist, open('zz_index_RAD_angle.ndx','w') as index_RAD_angle, open('zz_index_ADH_angle.ndx','w') as index_ADH_angle, open('zz_index_RAH_angle.ndx','w') as index_RAH_angle:
        index_RAD_angle.write('[ hbonds_RAD_angle' + str(count) + ' ] \n')
        index_RAH_angle.write('[ hbonds_RAH_angle' + str(count) + ' ] \n')
        index_ADH_angle.write('[ hbonds_ADH_angle' + str(count) + ' ] \n')

        for el_don, el_hyd, el_acc, el_ref in zip(don_atom, hyd_atom, acc_atom, ref_atom):

            index_AD_dist.write('[ hbonds_AD_dist_' + str(count) + ' ] \n')
            index_AH_dist.write('[ hbonds_AH_dist_' + str(count) + ' ] \n')

            index_AD_dist.write   ('  ' + str(el_acc) + '  ' + str(el_don) + '\n')
            index_AH_dist.write   ('  ' + str(el_acc) + '  ' + str(el_hyd) + '\n')
            index_RAD_angle.write ('  ' + str(el_ref) + '  ' + str(el_acc) + '  ' + str(el_don) + '\n')
            index_RAH_angle.write ('  ' + str(el_ref) + '  ' + str(el_acc) + '  ' + str(el_hyd) + '\n')
            index_ADH_angle.write ('  ' + str(el_acc) + '  ' + str(el_don) + '  ' + str(el_hyd) + '\n')
            count += 1
            
    don_atom = []
    hyd_atom = []
    acc_atom = []
    ref_atom = []
    HB_ndx_mat = []

    with open('in_for_distance.in', 'w') as dist_file_temp:
        for i in range(count):
            dist_file_temp.write(str(i) + ' \n')

            
    call('gmx distance -s %s.tpr -f %s_fit.xtc -n zz_index_AD_dist.ndx  -rmpbc -pbc -oall  zzz_AD_dist.xvg < in_for_distance.in' %(args.prefix, args.prefix), shell=True, stdout=None, stderr=None)
    call('gmx distance -s %s.tpr -f %s_fit.xtc -n zz_index_AH_dist.ndx  -rmpbc -pbc -oall  zzz_AH_dist.xvg < in_for_distance.in' %(args.prefix, args.prefix), shell=True, stdout=None, stderr=None)

    call('gmx angle              -f %s_fit.xtc -n zz_index_RAD_angle.ndx        -all -ov zzz_RAD_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)
    call('gmx angle              -f %s_fit.xtc -n zz_index_RAH_angle.ndx        -all -ov zzz_RAH_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)
    call('gmx angle              -f %s_fit.xtc -n zz_index_ADH_angle.ndx        -all -ov zzz_ADH_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)

    return howmany
    
def gen_hbond_ndx2(howmany):
    
    # if os.path.isfile('zzz_AD_dist.xvg') is False:
        
        print("INFO: Building zzz_AD_dist.xvg.")

        with open('zzz_AD_dist.xvg','w') as file_AD:
            file_AD.write('# This file was created Sat May 29 09:29:04 2021\n')
            file_AD.write('# Created by:\n')
            file_AD.write('#                      :-) GROMACS - gmx distance, 2018 (-:\n')
            file_AD.write('# \n')
            file_AD.write('# Executable:   /share/software/user/open/gromacs/2018/bin/gmx\n')
            file_AD.write('# Data prefix:  /share/software/user/open/gromacs/2018\n')
            file_AD.write('# Working dir:  /scratch/users/jkozuch/PYP/oTLN_wat\n')
            file_AD.write('# Command line:\n')
            file_AD.write('#   gmx distance -s TLN_wat_md.tpr -f TLN_wat_md_fit.xtc -n zz_index_AD_dist.ndx -rmpbc -pbc -oall zzz_AD_dist_temp_NUM.xvg - USING PYTHON SCRIPT\n')
            file_AD.write('# gmx distance is part of G R O M A C S:\n')
            file_AD.write('# \n')
            file_AD.write('# Gromacs Runs On Most of All Computer Systems\n')
            file_AD.write('# \n')
            file_AD.write('@    title \"Distance\"\n')
            file_AD.write('@    xaxis  label \"Time (ps)\"\n')
            file_AD.write('@    yaxis  label \"Distance (nm)\"\n')
            file_AD.write('@TYPE xy\n')
            
            
            for i in range(0, howmany):
                with open('zzz_AD_dist_temp_%d.xvg' % i,'r') as file_AD_num:
                    for line in file_AD_num:
                        if isnumber(line.split()[0]):
                            file_AD.write(line)

    # if os.path.isfile('zzz_AH_dist.xvg') is False:

        print("INFO: Building zzz_AH_dist.xvg.")
    
        with open('zzz_AH_dist.xvg','w') as file_AH:
            file_AH.write('# This file was created Sat May 29 09:29:04 2021\n')
            file_AH.write('# Created by:\n')
            file_AH.write('#                      :-) GROMACS - gmx distance, 2018 (-:\n')
            file_AH.write('# \n')
            file_AH.write('# Executable:   /share/software/user/open/gromacs/2018/bin/gmx\n')
            file_AH.write('# Data prefix:  /share/software/user/open/gromacs/2018\n')
            file_AH.write('# Working dir:  /scratch/users/jkozuch/PYP/oTLN_wat\n')
            file_AH.write('# Command line:\n')
            file_AH.write('#   gmx distance -s TLN_wat_md.tpr -f TLN_wat_md_fit.xtc -n zz_index_AH_dist.ndx -rmpbc -pbc -oall zzz_AH_dist_temp_NUM.xvg - USING PYTHON SCRIPT\n')
            file_AH.write('# gmx distance is part of G R O M A C S:\n')
            file_AH.write('# \n')
            file_AH.write('# Gromacs Runs On Most of All Computer Systems\n')
            file_AH.write('# \n')
            file_AH.write('@    title \"Distance\"\n')
            file_AH.write('@    xaxis  label \"Time (ps)\"\n')
            file_AH.write('@    yaxis  label \"Distance (nm)\"\n')
            file_AH.write('@TYPE xy\n')
            
            for i in range(0, howmany):
                with open('zzz_AH_dist_temp_%d.xvg' % i,'r') as file_AH_num:
                    for line in file_AH_num:
                        if isnumber(line.split()[0]):
                            file_AH.write(line)
                            
        call('rm *dist_temp*.xvg', shell=True)


#    p3 =  Popen('gmx angle              -f %s_fit.xtc -n zzz_index_RAD_angle.ndx        -all -ov zzz_RAD_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)
#    p4 =  Popen('gmx angle              -f %s_fit.xtc -n zzz_index_RAD_angle.ndx        -all -ov zzz_RAH_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)
#    p5 =  Popen('gmx angle              -f %s_fit.xtc -n zzz_index_RAD_angle.ndx        -all -ov zzz_ADH_angle.xvg' %(args.prefix), shell=True, stdout=None, stderr=None)
#    call('gmx distance -s %.tpr -f %.xtc -n zzz_index_AD_dist.ndx  -rmpbc -pbc -oall zzz_AD_dist.xvg   << EOT \n 0 \n EOT' % (args.prefix, args.prefix), shell=True, stdout=None, stderr=None)
#    call('gmx distance -s %.tpr -f %.xtc -n zzz_index_AH_dist.ndx  -rmpbc -pbc -oall zzz_AH_dist.xvg   << EOT \n 0 \n EOT' % (args.prefix, args.prefix), shell=True, stdout=None, stderr=None)
#    call('gmx angle    -s %.tpr -f %.xtc -n zzz_index_RAD_angle.ndx -rmpbc -pbc -all  zzz_RAD_angle.xvg << EOT \n 0 \n EOT' % (args.prefix, args.prefix), shell=True, stdout=None, stderr=None)
#    call('gmx angle    -s %.tpr -f %.xtc -n zzz_index_RAH_angle.ndx -rmpbc -pbc -all  zzz_RAH_angle.xvg << EOT \n 0 \n EOT' % (args.prefix, args.prefix), shell=True, stdout=None, stderr=None)
#    call('gmx angle    -s %.tpr -f %.xtc -n zzz_index_ADH_angle.ndx -rmpbc -pbc -all  zzz_ADH_angle.xvg << EOT \n 0 \n EOT' % (args.prefix, args.prefix), shell=True, stdout=None, stderr=None)

def analyze_Hbonds():
    HB_num_mat = []
    
    call("sed -i '/^#/d' zz_Hbond_num.xvg", shell=True)
    call("sed -i '/^@/d' zz_Hbond_num.xvg", shell=True)
    
    call("sed -i '/^#/d' zzz_AD_dist.xvg", shell=True)
    call("sed -i '/^@/d' zzz_AD_dist.xvg", shell=True)
    
    call("sed -i '/^#/d' zzz_AH_dist.xvg", shell=True)
    call("sed -i '/^@/d' zzz_AH_dist.xvg", shell=True)
    
    call("sed -i '/^#/d' zzz_ADH_angle.xvg", shell=True)
    call("sed -i '/^@/d' zzz_ADH_angle.xvg", shell=True)
    
    call("sed -i '/^#/d' zzz_RAD_angle.xvg", shell=True)
    call("sed -i '/^@/d' zzz_RAD_angle.xvg", shell=True)
    
    call("sed -i '/^#/d' zzz_RAH_angle.xvg", shell=True)
    call("sed -i '/^@/d' zzz_RAH_angle.xvg", shell=True)
    
    with open('zz_Hbond_num.xvg','r') as file_HB_num, open('zzz_AD_dist.xvg','r') as file_AD_dist, open('zzz_ADH_angle.xvg','r') as file_ADH_angle, open('zzz_AH_dist.xvg','r') as file_AH_dist, open('zzz_RAD_angle.xvg','r') as file_RAD_angle, open('zzz_RAH_angle.xvg','r') as file_RAH_angle,\
         open('zzz_HB.txt','w') as file_HB, open('Electric_Fields_on_CN.txt', 'r') as file_EF, open('zzz_HB_EF.txt', 'w') as file_HB_EF:

        HB_num_man = []
        frame = 0
        for line_HB_num, line_AD_dist, line_AH_dist, line_ADH_angle, line_RAD_angle, line_RAH_angle, line_EF in zip(file_HB_num, file_AD_dist, file_AH_dist, file_ADH_angle, file_RAD_angle, file_RAH_angle, file_EF):
            frame += 1
            if isnumber(line_HB_num.split()[0]):
                time  = line_HB_num.split()[0]
                HB_num = np.array([int(x) for x in line_HB_num.split()[1]])

            if isnumber(line_AD_dist.split()[0]):
                time  = line_AD_dist.split()[0]
                HB_AD_dist = np.array([float(x) for x in line_AD_dist.split()[1::]])

            if isnumber(line_AH_dist.split()[0]):
                time  = line_AH_dist.split()[0]
                HB_AH_dist = np.array([float(x) for x in line_AH_dist.split()[1::]])

            if isnumber(line_ADH_angle.split()[0]):
                time  = line_ADH_angle.split()[0]
                HB_ADH_angle = np.array([float(x) for x in line_ADH_angle.split()[2::]])

            if isnumber(line_RAD_angle.split()[0]):
                time  = line_RAD_angle.split()[0]
                HB_RAD_angle = np.array([float(x) for x in line_RAD_angle.split()[2::]])

            if isnumber(line_RAH_angle.split()[0]):
                time  = line_RAH_angle.split()[0]
                HB_RAH_angle = np.array([float(x) for x in line_RAH_angle.split()[2::]])
 
            if isnumber(line_EF.split()[0]):
                time  = line_EF.split()[0]
                EF = float(line_EF.split()[3])    
 
            # elements = len(HB_AD_dist)

            # for i in range(0, elements):
            #     if (float(HB_AD_dist[i]) <= float(args.dist) and float(HB_ADH_angle[i]) <= float(args.angle)) is False:
    
            #         HB_AD_dist[i]   = float(999)
            #         HB_ADH_angle[i] = float(999)
            #         HB_AH_dist[i]   = float(999)
            #         HB_RAD_angle[i] = float(999)
            #         HB_RAH_angle[i] = float(999)

            HB_AD_dist_sel   = []
            HB_ADH_angle_sel = []
            HB_AH_dist_sel   = []
            HB_RAD_angle_sel = []
            HB_RAH_angle_sel = []
            
            for dist, angle, AH_dist, RAD_angle, RAH_angle in zip(HB_AD_dist, HB_ADH_angle, HB_AH_dist, HB_RAD_angle, HB_RAH_angle):
                if float(dist) <= float(args.dist) and float(angle) <= float(args.angle):
                    HB_AD_dist_sel.append(dist)
                    HB_ADH_angle_sel.append(angle)
                    HB_AH_dist_sel.append(AH_dist)
                    HB_RAD_angle_sel.append(RAD_angle)
                    HB_RAH_angle_sel.append(RAH_angle)    
                    
            str_1 = str(len(HB_AD_dist_sel))
            str_2 = ''
            if len(HB_AD_dist_sel) > 0:
                for dist, angle in zip(HB_AD_dist_sel,HB_RAD_angle_sel):
                    str_2 = str_2 + str(round(dist,3)).rjust(9) + str(round(angle,3)).rjust(9) + '    ATOM  xxxxx  X   XXX X xxx  '
            file_HB.write(str_1 + str_2 + '\n')
            
            str_1 = str(frame).rjust(6)
            str_2 = str(round(EF,2)).rjust(9)
            str_3 = ''
            if len(HB_AD_dist_sel) > 0:
                for dist, angle in zip(HB_AD_dist_sel,HB_RAD_angle_sel):
                    str_3 = str_3 + str(round(dist,3)).rjust(9) + str(round(angle,3)).rjust(9) + '    ATOM  xxxxx  X   XXX X xxx  '
            file_HB_EF.write(str_1 + str_2 + str_3 + '\n')


############
### Main ###
############

global args
args = parse()

if os.path.isfile('zzz_AD_dist.xvg') is False or\
   os.path.isfile('zzz_AH_dist.xvg') is False or\
   os.path.isfile('zzz_ADH_angle.xvg') is False or\
   os.path.isfile('zzz_RAH_angle.xvg') is False or\
   os.path.isfile('zzz_RAD_angle.xvg') is False:

    call('rm zzz_*.*', shell=True, stdout=None, stderr=None)
    howmany = gen_hbond_ndx()
    print('################################################################################################')
    print("INFO: All data collected.")
    #print('################################################################################################')
    
    #gen_hbond_ndx2(howmany)
    #print('################################################################################################')
    print('INFO: Done generating ndx and xvg files.')
    print('################################################################################################')

print('################################################################################################')
print('INFO: Analyzing H-bonds.')
print('################################################################################################')

analyze_Hbonds()
print('################################################################################################')
print('INFO: Finished analyzing H-bond. Double check if you have the same number of H-bonds as in gmx hbonds.')
print('################################################################################################')















