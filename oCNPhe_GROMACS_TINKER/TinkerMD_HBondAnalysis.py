# Author Jacek Kozuch, kozuch.jacek@gmail.com (based on Richard Bradshaw's and Stephen Frieds scripts)
# Script to calculate the Amoeba-based electric field at a desired site
# Requirements: Tinker with modified analyze function (added option U: copy-paste option D and turn debug = .false.)

# usage example: 
# python TinkerMD_CN_Hbond.py -st PYP_F28oCNF_A_md -a 402 403
# python TinkerMD_CN_Hbond.py -st PYP_F92oCNF_A_md -a 1361 1362

import os
import sys, argparse, glob
import re
import itertools
import numpy as np
from subprocess import call,Popen
from datetime import datetime
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from os.path import exists as file_exists
#import mdtraj as md
import psutil

#################
### LogFile   ###
#################
def log_write(string):
    now    = datetime.now()
    dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
    with open('zzz_log.txt','a') as log:
        log.write(dt_str + ' ' + string + '\n')

#################
### Argparser ###
#################
def parse():
    log_write('')
    log_write('##############################################################')
    log_write('#                                                            #')
    log_write('# Title: TinkerMD_CN_Hbond                                   #')
    log_write('# Description: Analysis of H-bonding after Tinker MD         #')
    log_write('#                                                            #')
    log_write('# Copyright:            Copyright (c) Jacek Kozuch 2021-2022 #')
    log_write('#                                                            #')
    log_write('##############################################################')
    log_write('')
    log_write('')

    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 2). First C then N.",nargs='+',type=str,required=True)
#    parser.add_argument("-pa","--proteinatoms",help="Provide the range for the protein.",nargs='+',type=int)
#    parser.add_argument("-sa","--solventatoms",help="Provide the range for the solvent and ions.",nargs='+',type=int)
    parser.add_argument("-arc","--traj",help="Filename of Solvated Trajectory (Tinker arc).",type=str,required=True)
    parser.add_argument("-md","--maxdist",help="Maximum H-bond A-D distance.",type=float,default=3.5)
    parser.add_argument("-ma","--maxangle", help="Maximum H-bond A-D-H angle.",type=float,default=30)
    parser.add_argument("-be","--beginend", help="first and last frame.",nargs='+',type=int,default=[1,99999999999999])


    if len(sys.argv)==1:
       parser.print_help()
       sys.exit(1)  

    args = parser.parse_args()
    return args


args = parse()
bond_atoms   = args.atoms
arcname      = args.traj
arcfile      = args.traj + '.arc'
uindfile     = args.traj + '.uind'
#idx_prot     = range(args.proteinatoms[0], args.proteinatoms[1]+1)
#if args.solventatoms[1] == -1:
#    idx_solv     = range(args.solventatoms[0], args.solventatoms[1])
#else:
#    idx_solv     = range(args.solventatoms[0], args.solventatoms[1]+1)
max_dist    = args.maxdist
max_angle   = args.maxangle
atom_dict = {}

accep_file = 'zzz_accep.txt'
refer_file = 'zzz_refer.txt'
donor_file = 'zzz_donor.txt'
hydro_file = 'zzz_hydro.txt'

result_Hbond = 'zzz_HBONDS.txt'
result_RDAhist = 'zzz_RDA_hist.txt'
result_DAhist = 'zzz_DA_hist.txt'


#frames = 2500

R_atom = []
A_atom = [] 
D_atom = []
H_atom = []
pbcbox = []
D_info = []


##################
### Other Defs ###
##################


### Check if numbers in line and not header ###
def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def RAMusage():
    process = psutil.Process(os.getpid())
    RAM = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    return str(RAM)


def flatten(input):
    new_list = []
    for i in input:
        for j in i:
            new_list.append(j)
    return new_list

def dist(in_xyz_atom1, in_xyz_atom2):

    xyz_atom1 = np.array([float(i) for i in in_xyz_atom1])
    xyz_atom2 = np.array([float(i) for i in in_xyz_atom2])

    vector = xyz_atom2 - xyz_atom1
    if len(vector) == 1:
        vector = vector[0]
    distance = np.sqrt(np.dot(vector, vector))
    return distance

def check_closest(in_xyz_ref, in_xyz_check, in_pbcbox):
    
    pbcbox    = [float(i) for i in in_pbcbox]
    xyz_ref   = [float(i) for i in in_xyz_ref]
    xyz_check = [float(i) for i in in_xyz_check]

    mirrors =          [ [-pbcbox[0],  pbcbox[1],  pbcbox[2]],
                         [         0,  pbcbox[1],  pbcbox[2]], 
                         [ pbcbox[0],  pbcbox[1],  pbcbox[2]],
                         [-pbcbox[0],          0,  pbcbox[2]],
                         [         0,          0,  pbcbox[2]],
                         [ pbcbox[0],          0,  pbcbox[2]],
                         [-pbcbox[0], -pbcbox[1],  pbcbox[2]],
                         [         0, -pbcbox[1],  pbcbox[2]],
                         [ pbcbox[0], -pbcbox[1],  pbcbox[2]],
                         [-pbcbox[0],  pbcbox[1],          0],
                         [         0,  pbcbox[1],          0],
                         [ pbcbox[0],  pbcbox[1],          0],
                         [-pbcbox[0],          0,          0],
                         [         0,          0,          0],
                         [ pbcbox[0],          0,          0],
                         [-pbcbox[0], -pbcbox[1],          0],
                         [         0, -pbcbox[1],          0],
                         [ pbcbox[0], -pbcbox[1],          0],
                         [-pbcbox[0],  pbcbox[1], -pbcbox[2]],
                         [         0,  pbcbox[1], -pbcbox[2]], 
                         [ pbcbox[0],  pbcbox[1], -pbcbox[2]],
                         [-pbcbox[0],          0, -pbcbox[2]],
                         [         0,          0, -pbcbox[2]],
                         [ pbcbox[0],          0, -pbcbox[2]],
                         [-pbcbox[0], -pbcbox[1], -pbcbox[2]],
                         [         0, -pbcbox[1], -pbcbox[2]],
                         [ pbcbox[0], -pbcbox[1], -pbcbox[2]] ]

    mirratoms  = np.array(xyz_check) + np.array(mirrors)

    dists = []
    for xyz_atom in mirratoms:
        dists.append(dist(np.array(xyz_ref), xyz_atom))
    i_min  = np.where(dists == min(dists))
    result =  mirratoms[i_min].tolist()
    return next(iter(result))
  
def atomwithin(in_xyz_ref, in_xyz_check, in_pbcbox):

    pbcbox    = [float(i) for i in in_pbcbox]
    xyz_ref   = [float(i) for i in in_xyz_ref]
    xyz_check = [float(i) for i in in_xyz_check]

    global args

    xyz_atom = check_closest(xyz_ref, xyz_check, pbcbox)
    if dist(xyz_ref, xyz_atom) <= args.maxdist:
        return True

def if_H(idx):
    global atom_dict
    
    if atom_dict[idx] == 'H' or atom_dict[idx] == 'HN' or atom_dict[idx] == 'HO' or atom_dict[idx] == 'HS':
        return True

def if_any_H(idxs):
    global atom_dict
    for idx in idxs:
        if atom_dict[idx] == 'H' or atom_dict[idx] == 'HN' or atom_dict[idx] == 'HO' or atom_dict[idx] == 'HS':
            return True

def calc_angle(in_a1,in_a2,in_b1,in_b2):

    a1    = np.array([float(i) for i in in_a1])
    a2    = np.array([float(i) for i in in_a2])
    b1    = np.array([float(i) for i in in_b1])
    b2    = np.array([float(i) for i in in_b2])

    vec_a = a2 - a1
    vec_a = vec_a/np.sqrt(np.dot(vec_a,vec_a))
    vec_b = b2 - b1
    vec_b = vec_b/np.sqrt(np.dot(vec_b,vec_b))
    theta = np.arccos(np.dot(vec_a,vec_b))*180/np.pi
    return theta

def anglewithin(in_a1,in_a2,in_b1,in_b2):
    global args

    a1    = np.array([float(i) for i in in_a1])
    a2    = np.array([float(i) for i in in_a2])
    b1    = np.array([float(i) for i in in_b1])
    b2    = np.array([float(i) for i in in_b2])

    theta = calc_angle(a1,a2,b1,b2)
    if theta <= args.maxangle:
        return True

def calc_distance(in_a1, in_a2):
    global args

    a1    = np.array([float(i) for i in in_a1])
    a2    = np.array([float(i) for i in in_a2])
    vec   = a2 - a1
    return np.sqrt(np.dot(vec,vec)) 
       
def distancewithin(in_a1, in_a2):
    global args

    a1    = np.array([float(i) for i in in_a1])
    a2    = np.array([float(i) for i in in_a2])
    
    if calc_distance(a1, a2) <= 1.5:
       return True


def reshape(list):
    new_list = []
    for i in range(0,int(len(list)/3)):
        new_list.append(list[3*i:3*i+3])
    return new_list

def reshape2(list):
    new_list = []
    for i in range(0,int(len(list)/2)):
        new_list.append(list[2*i:2*i+2])
    return new_list



#################
### Main Defs ###
#################

def make_dict():
    global args, atom_dict
    arcfile      = args.traj + '.xyz'
    key   = []
    value = []
    with open(arcfile, 'r') as file:
        for line in file:
            if len(line.split()) > 1 and '.' not in line.split()[0]:
                key.append(line.split()[0])
                value.append(line.split()[1])
    atom_dict = dict(zip(key, value))

    log_write('INFO: Protein sequence read and seems all good. Will contine writing frames.')
    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())

    

def reduce_arc():
    global args, atom_dict, probe_file, donor_file, hydro_file, frames

    log_write('INFO: Looking for probe atoms.')

    pbcbox  = []
    xyz_ref = []
    xyz_acc = []
    acc     = []
    ref     = []
    frame  = 0

    arcfile      = args.traj + '.arc'
    xyzfile      = args.traj + '.xyz'

    with open(xyzfile, 'r') as file:
        for line in file:
            first_line = line
            break
#    with open(arcfile, 'r') as file:
#        for line in file:
#            if line == '\n':
#                break
#            else:
#                last_line = line
    check = False

    with open(arcfile, 'r') as file, open(accep_file, 'w') as out_file_A, open(refer_file, 'w') as out_file_R:
        for line in file:
                if line == first_line:

                    frame  = frame + 1
                    ref = []
                    acc = []
                    ref = [ frame ]
                    acc = [ frame ]
                    if (frame >= args.beginend[0]) and (frame <= args.beginend[1]):
                        check = True
                        print("Will be using frame " + str(frame))
                    else:
                        log_write('INFO: Avoiding frame ' + str(frame) + '.')
                        check = False                         
                elif "." in line.split()[0] and check == True:
                    log_write('INFO: Will be evaluating frame ' + str(frame) + '.')
                    ref = ref + [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
                    acc = acc + [float(line.split()[0]), float(line.split()[1]), float(line.split()[2])]
                    pbcbox.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
                elif line.split()[0] == args.atoms[0] and check == True:
                    ref = ref + [ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
                    xyz_ref.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
                    out_file_R.write("  ".join(str(x) for x in ref) + '\n')
                elif line.split()[0] == args.atoms[1] and check == True:
                    acc = acc + [ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
                    xyz_acc.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
                    out_file_A.write("  ".join(str(x) for x in acc) + '\n')

                if (frame > args.beginend[1]):
                    break

    log_write('INFO: Looking for donor and hydrogen atoms.')
    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())


    xyz_don = []
    xyz_hyd = []
    xyz_H   = []
    H_rem  = []
    frame  = 0
    frame2 = 0
    check = False
    check2 = False
    with open(arcfile, 'r') as file, open(donor_file, 'w') as out_file_D, open(hydro_file, 'w') as out_file_H:
        for line in file:
            if line == first_line:
                if check2 == True:
                    out_file_D.write("  ".join(str(x) for x in xyz_don) + '\n')
                    out_file_H.write("  ".join(str(x) for x in xyz_hyd) + '\n')                   
                frame = frame + 1
                frame2 = frame2 + 1
                xyz_don = []
                xyz_hyd = []
                xyz_don = [ frame ]
                xyz_hyd = [ frame ]
                check = True
                if (frame >= args.beginend[0]) and (frame <= args.beginend[1]):
                    check2 = True
                    print("Will be using frame " + str(frame))
                else:
                    log_write('INFO: Avoiding frame ' + str(frame) + '.')
                    frame2 = frame2 - 1
                    check2 = False     
            elif "." in line.split()[0]:
                continue
            elif line.split()[1] in ['N', 'O', 'S', 'NH', 'OH', 'SH'] and line.split()[0] not in args.atoms and if_any_H(line.split()[6:]) and check2 == True:
                xyz_atom = [ float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
                if atomwithin( xyz_acc[frame2-1], xyz_atom, pbcbox[frame2-1] ):
                    xyz_don = xyz_don + [ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
            elif line.split()[1] in ['H','HN','HO','HS']  and check2 == True:
                xyz_atom = [ float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
                if atomwithin( xyz_acc[frame2-1], xyz_atom, pbcbox[frame2-1] ):
                    xyz_hyd = xyz_hyd + [ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]) ]
            if (frame > args.beginend[1]):
                break

        out_file_D.write("  ".join(str(x) for x in xyz_don) + '\n')
        out_file_H.write("  ".join(str(x) for x in xyz_hyd) + '\n')  


    log_write('INFO: All atoms found!')
    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())

               
def collect_atoms():
    global args, atom_dict, probe_file, donor_file, hydro_file, frames, R_atom, A_atom, D_atom, D_info, H_atom, pbcbox


    log_write('INFO: Collect all relevant atom positions and check if they are within PBC, otherwise correct!')
 
    file_A  = open(accep_file, 'r')
    file_R  = open(refer_file, 'r')
    file_D  = open(donor_file, 'r')
    file_H  = open(hydro_file, 'r')

    pbcbox = []
    R_atom = []
    A_atom = [] 
    D_atom = []
    H_atom = []
    temp   = []
    temp_info = []
    D_info = []

    for line in file_A:
        pbcbox.append(line.split()[1:4])
        A_atom.append(line.split()[6:])
    for line, PBC, A in zip(file_R, pbcbox, A_atom):
        temp = check_closest(A, line.split()[6:], PBC)
        R_atom.append(temp)
        temp = []
    for line, PBC, A in zip(file_D, pbcbox, A_atom): 
        howmany = int((len(line.split()) - 1)/5)
        if howmany == 0:
            D_atom.append([])
            D_info.append([])
        else:
            for a in range(0, howmany):
                start =  3 + 5*a
                stop  =  6 + 5*a
                result = check_closest(A, line.split()[start:stop], PBC)
                temp.append(result)
                temp_info.append(line.split()[start-2:start])
            D_atom.append(flatten(temp))
            D_info.append(flatten(temp_info))
        temp = []
    for line, PBC, A in zip(file_H, pbcbox, A_atom): 
        howmany = int((len(line.split()) - 1)/5)
        if howmany == 0:
            H_atom.append([])
        else:
            for a in range(0, howmany):
                start =  3 + 5*a
                stop  =  6 + 5*a
                result = check_closest(A, line.split()[start:stop], PBC)
                temp.append(result)
            H_atom.append(flatten(temp))
        temp = []

    file_A.close()
    file_R.close()
    file_D.close()
    file_H.close()

    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())

        
def analyze_Hbonds():

    log_write('INFO: Analyze Hbonds!')
    global args, atom_dict, probe_file, donor_file, hydro_file, frames, R_atom, A_atom, D_atom, D_info, H_atom, pbcbox

    H_num     = []
    DA_dist   = []
    RDA_angle = []
    H_cone    = []
    D_which   = []

    H_temp    = []
    DA_temp   = []
    RDA_temp  = []
    info_temp = []
    
    print(len(R_atom))
    print(len(A_atom))
    print(len(D_atom))
    print(len(H_atom))
    print(len(D_info))
    
    count = 0

    for A, R, D, H, info in zip(A_atom, R_atom, D_atom, H_atom, D_info):
        count += 1
        for don, which in zip(reshape(D), reshape2(info)):
            for hyd in reshape(H):
                if anglewithin(don,A,don,hyd) and distancewithin(don,hyd) and dist(A, don) > dist(A, hyd):
                    H_temp.append(calc_angle(don,A,don,hyd))
                    DA_temp.append(calc_distance(don, A))
                    RDA_temp.append(calc_angle(A,R,A,don))
                    info_temp.append(' '.join(x for x in which))            

        if H_temp == 'None':
            H_cone.append([0])
            DA_dist.append([0])
            RDA_angle.append([0])
            D_which.append([0])
            # print(count)
        else:
            H_cone.append(H_temp)
            DA_dist.append(DA_temp)
            RDA_angle.append(RDA_temp)
            D_which.append(info_temp)
            # print(count)

        H_temp = []
        DA_temp = []
        RDA_temp = []
        info_temp = []

    for hbond in DA_dist:
        if any(x != 0 for x in hbond):
            H_num.append(len(hbond))
        else:
            H_num.append(0)

    DA_flat = flatten(DA_dist)
    DA_hist = np.histogram(DA_flat, 35, range=(0, 3.5), density = True)
    RDA_flat = flatten(RDA_angle)
    RDA_hist = np.histogram(RDA_flat, 36, range=(0, 180), density = True)
    
    print(len(H_num))
    print(len(DA_dist))
    print(len(RDA_angle))
    print(len(D_which))

    with open(result_Hbond, 'w') as file:
        for N, DA, RDA, which in zip(H_num, DA_dist, RDA_angle, D_which):
            file.write(str(N) + ' ' + ' '.join(str(np.round(x,2)).rjust(15) for x in DA) + ' ' + ' '.join(str(np.round(x,2)).rjust(15) for x in RDA) + ' ' + ' '.join(str(x).rjust(10) for x in which) + '\n')            

    with open(result_DAhist, 'w') as file:
        for bins, counts in zip(DA_hist[1], DA_hist[0]):
            file.write(str(np.round(bins,2)).rjust(6) + '    ' + str(counts/10)  + '\n') 

    with open(result_RDAhist, 'w') as file:
        for bins, counts in zip(RDA_hist[1], RDA_hist[0]):
            file.write(str(np.round(bins,2)).rjust(6) + '    ' + str(counts*5)  + '\n') 

    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())
    


### MAIN ###
make_dict()
if file_exists(accep_file) and file_exists(refer_file) and file_exists(donor_file) and file_exists(hydro_file):
    log_write('INFO: Reduced trajectory files found. Continue.')
else:
    log_write('INFO: Reduced trajectory files not found. Preparing')
    reduce_arc()
collect_atoms()
analyze_Hbonds()    

log_write('INFO: DONE!')
