# Author Jacek Kozuch, kozuch.jacek@gmail.com, 2022/08/05 
# Script to translate TINKER arc file into PDB

import time
import os
import sys, argparse, glob
import re
import itertools
import numpy as np
from subprocess import call,Popen
from datetime import datetime
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as Rot
from os.path import exists as file_exists
import psutil

#import mdtraj as md

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
    log_write('# Title: TinkerMD_xyz2pdb                                    #')
    log_write('# Description: Translating ARC File after Tinker MD          #')
    log_write('#                                                            #')
    log_write('# Copyright: Jacek Kozuch - 2022/06/05                       #')
    log_write('#                                                            #')
    log_write('##############################################################')
    log_write('')
    log_write('')

    parser = argparse.ArgumentParser()
    parser.add_argument("-arc","--arc",help="Filename of Solvated Trajectory (Tinker arc).",type=str,required=True)
    parser.add_argument("-xyz","--xyz",help="Filename of Solvated xyz file (Tinker xyz).",type=str,default=[])
    parser.add_argument("-f","--fit",help="Fit to CA, backbone heavy atoms, or specified atoms (ca, bb, spec).",type=str,default = 'ca')
    parser.add_argument("-al","--atomlist",help="Provide list of atomnumber that fit should be performed to.", type=str, nargs='+')
    parser.add_argument("-p","--precision",help="Provide precision of how many decimal places should be used.", type=int, default=3)
    parser.add_argument("-a","--atomsforwaters",help="Atom number to select waters.", type=int, nargs='+',default = [])
    parser.add_argument("-d","--distanceofwaters",help="Cut off distance to select waters in Angstrom.", type=float, nargs='+',default = [])
    parser.add_argument("-be","--beginend", help="first and last frame.",nargs='+',type=int,default=[1,99999999999999])


    if len(sys.argv)==1:
       parser.print_help()
       sys.exit(1)  

    args = parser.parse_args()
    return args


args = parse()
atom_dict = {}
resid_dict = {}
resid_dict_temp = {}
firstidxsolvent = 0
if args.xyz == []:
    args.xyz = args.arc 




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

def flatten(input):
    new_list = []
    for i in input:
        for j in i:
            new_list.append(j)
    return new_list

def RAMusage():
    process = psutil.Process(os.getpid())
    RAM = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    return str(RAM)

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
  
def precision(in_coord):

    global args
    coord = np.array([round(float(i),args.precision) for i in in_coord])
    return coord.tolist()

def findgeocenter(in_array, lastidx):
    
    global args
    sum_coord = np.array([0, 0, 0])
    for i in range(0, lastidx):
        npcoord = np.array([float(i) for i in in_array[i][2:5]])
        sum_coord = sum_coord + npcoord
    center = sum_coord/(lastidx-1)
    return center.tolist()
    
def centersystem(in_system,in_center,pbcbox,translate):

    global args, atom_dict
    new_system = []
    for atom in in_system:
        new_coords = np.array(check_closest(in_center,atom[2:5],pbcbox[0:3])) + np.array(translate)
        new_system.append([atom[0], atom[1], float(new_coords[0]), float(new_coords[1]), float(new_coords[2]), atom[5] ])

    return new_system

def get_rotmat(in_xyz,in_xyz_ref):
    global tstart  

    xyz = np.array(in_xyz)
    xyz_ref = np.array(in_xyz_ref)

    def RotMat(alpha,beta,gamma):

         Rz_alpha = np.array([[np.cos(alpha), -np.sin(alpha), 0            ], [ np.sin(alpha),  np.cos(alpha),  0            ], [ 0           , 0            , 1            ]])
         Ry_beta  = np.array([[np.cos(beta) ,  0            , np.sin(beta) ], [ 0            ,  1            ,  0            ], [-np.sin(beta), 0            , np.cos(beta) ]])
         Rx_gamma = np.array([[1            ,  0            , 0            ], [ 0            ,  np.cos(gamma), -np.sin(gamma)], [ 0           , np.sin(gamma), np.cos(gamma)]])
         
         return np.matmul(Rz_alpha,np.matmul(Ry_beta,Rx_gamma))    

    def rotate(xyz, alpha, beta, gamma, t_x, t_y, t_z):  
    
         R = RotMat(alpha,beta,gamma)
         return flatten(np.matmul(R,xyz.transpose()).transpose() - np.array([t_x, t_y, t_z]))

    popt, pcov = curve_fit(rotate, xyz, flatten(xyz_ref), [0, 0, 0, 0, 0, 0])

    RotMatrix = RotMat(popt[0], popt[1], popt[2])
    TransVec  = [popt[3], popt[4], popt[5]]

    return RotMatrix, TransVec


def fitsystem(in_current,in_first):

    global args, atom_dict, tstart, firstidxsolvent

    relevant_current = [] 
    relevant_first   = []

    if args.fit == 'ca':
        for atom_current, atom_first in zip(in_current[0:firstidxsolvent], in_first[0:firstidxsolvent]):
            if atom_current[1] == 'CA' and atom_current[5] in ['8', '2', '51','442','449']:
                relevant_current.append(atom_current[2:5])
                relevant_first.append(atom_first[2:5])

    elif args.fit == 'bb':
        for atom_current, atom_first in zip(in_current[0:last], in_first[0:last]):
            if atom_current[1] in ['N', 'CA', 'C', 'O'] and atom_current[5] in ['1', '2', '3', '5', '7', '8', '9', '11', '50', '51', '52', '53', '231', '233', '234','239','241','242','243','441','442','443','444','447','448','449','450']:
                relevant_current.append(atom_current[2:5])
                relevant_first.append(atom_first[2:5])

#    elif args.fit == 'spec':
         # needs to be added

    [RotMatrix, TransVec] = get_rotmat(relevant_current, relevant_first)

    current = (np.array(in_current).transpose()[2:5]).astype(float)
    new_current = np.matmul(RotMatrix, current).transpose() - np.array(TransVec)

    output = []
    for in_atom, atom in zip(in_current,new_current):
        output.append([in_atom[0], in_atom[1], atom[0], atom[1], atom[2], in_atom[5]])

    return output

def printPDB(out_file, frame, x_frame_centered):

    global resid_dict
    
    out_file.write('MODEL' + str(frame).rjust(9) + '\n')
    res = 0
    for atom in x_frame_centered[0:firstidxsolvent]:
        if atom[5] in ['1','7','51','231','239','441','447','494','522','543','']:
            res = res + 1
        coords = precision(atom[2:5])
        out_file.write('ATOM' + atom[0].rjust(7) + '  ' + atom[1].ljust(4) + resid_dict[str(res)] + ' A' + str(res).rjust(4) + '    ' + ''.join(str(x).rjust(8) for x in coords) + '  1.00  0.00' + '\n')

    out_file.write('TER '  + str(int(atom[0])+1).rjust(7) + '      ' + resid_dict[str(res)] + ' A' + str(res).rjust(4) + '\n')

    for atom in x_frame_centered[firstidxsolvent-1:]:
        if atom[5] in ['349','352','361']:
            res = res + 1
        coords = precision(atom[2:5])
        out_file.write('ATOM' + atom[0].rjust(7) + '  ' + atom[1].ljust(4) + 'SOL' + ' A' + str(res).rjust(4) + '    ' + ''.join(str(x).rjust(8) for x in coords) + '  1.00  0.00' + '\n')

    out_file.write('ENDMDL')
    if isinstance(int(frame)/100, int):
        log_write('INFO: Frame ' + str(frame) + ' written.')

def check_AAcode(line):

    resids = {   '2' : 'GLY',
                '13' : 'ALA',
                '15' : 'VAL',
                '19' : 'LEU',
                '25' : 'ILE',
                '33' : 'SER',
                '37' : 'THR',
                '43' : 'CYS',
                '48' : 'CYX',
                '51' : 'PRO',
               '241' : 'PRO',
                '61' : 'PHE',
                '70' : 'TYR',
                '80' : 'TZX',
                '89' : 'TRP',
               '106' : 'HIS',
               '117' : 'HID',
               '127' : 'HIE',
               '137' : 'ASP',
               '141' : 'ASH',
               '147' : 'ASN',
               '153' : 'GLU',
               '159' : 'GLH',
               '167' : 'GLN',
               '175' : 'MET',
               '182' : 'LYS',
               '192' : 'LYX',
               '202' : 'ARG',
               '213' : 'ORN',
               '401' : 'PCA',
               '417' : 'PCF',
               '427' : 'OCF',
               '442' : 'GLY',
               '452' : 'TYE',
               '462' : 'SPG',
               '494' : 'PEN',
               '522' : 'AVB',
               '543' : 'NOR',
               '579' : 'SAV'}
#               '349' : 'SOL',
#               '352' : 'NA+',
#               '361' : 'CL-'}
    
    if line.split()[5] in resids:
        return resids[line.split()[5]]
    else:
        return False



#################
### Main Defs ###
#################

def make_dict():
    global args, atom_dict, resid_dict, pbcbox, firstidxsolvent,resid_dict_temp

    arcfile    = args.xyz + '.xyz'
    key_atom   = []
    value_atom = []
    key_AA     = []
    value_AA   = []
    res = 0
    
    write = False
    with open(arcfile, 'r') as file:
        for line in file:
            line1 = line
            break

    with open(arcfile, 'r') as file:
        count = 1
        for line in file:
            line2 = line
            count += 1
            if count == 2:
                break
     
    with open(arcfile, 'r') as file:
        print('INFO: Process uses # MB RAM: ' + RAMusage())
        for line in file:
            if line != line1 and line != line2:
#            if len(line.split()) > 1 and '.' not in line.split()[0]:
                key_atom.append(line.split()[0])
                value_atom.append([ line.split()[1], line.split()[5] ])
                if firstidxsolvent == 0 and line.split()[5] == '349':
                    firstidxsolvent = int(line.split()[0])
                if check_AAcode(line) != False:
                    res = res + 1
                    key_AA.append(str(res))
                    value_AA.append(check_AAcode(line))

    atom_dict  = dict(zip(key_atom, value_atom))
    resid_dict = dict(zip(key_AA,   value_AA))
    print(resid_dict)
    log_write('INFO: Protein sequence read and seems all good. Will contine writing frames.')
    log_write('INFO: Process uses # MB RAM: ' + RAMusage())
    print('INFO: Process uses # MB RAM: ' + RAMusage())

def select_waters(current_frame):
    global args, atom_dict, pbcbox, firstidxsolvent, tstart,resid_dict,resid_dict_temp
    
    resid_dict_temp = resid_dict
    remove = []
    if args.atomsforwaters == []:
        return current_frame
    else:
        for sol in current_frame[firstidxsolvent:]:
            for atom,distance in zip(args.atomsforwaters,args.distanceofwaters):
                if sol[1] == 'O':
                    if dist(sol[2:5],current_frame[int(atom)-1][2:5]) > distance:
                        remove.append(int(sol[0])-1)
                        remove.append(int(sol[0]))
                        remove.append(int(sol[0])+1)
#                    else:
#                        print(dist(sol[2:5],current_frame[int(atom)][2:5]))
#                        print('This should be kept')

                if sol[1] in ['Na+','Cl-']:
                    if dist(sol[2:5],current_frame[int(atom)-1][2:5]) > distance:
                        remove.append(int(sol[0])-1)

    for idx in remove[::-1]:
        current_frame.pop(idx)
    return current_frame
        
def reduce_arc():
    global args, atom_dict, pbcbox, firstidxsolvent, tstart
    frame = 0
    frames = 0

    ref_frame        = [] 
    ref_center       = []
    current_frame    = []

    arcfile      = args.arc + '.arc'
    xyzfile      = args.xyz + '.xyz'

    with open(xyzfile, 'r') as file:
        for line in file:
            first_line = line
            break
    with open(arcfile, 'r') as file:
        for line in file:
            if line == '\n':
                break
            else:
                last_line = line
            

    with open(xyzfile, 'r') as file_xyz, open(arcfile, 'r') as file_arc:
        count_lines_xyz = 0
        count_lines_arc = 0
        for line_xyz in file_xyz:
            count_lines_xyz += 1
        for line_arc in file_arc:
            count_lines_arc += 1
        total_frames = int(count_lines_arc/count_lines_xyz)

 
    with open(xyzfile, 'r') as file:
        for line in file:
            if line == first_line:
                numberatoms = int(line.split()[0])
            elif '.' in line.split()[0]:
                pbcbox = [float(i) for i in line.split()]
            else:
                ref_frame.append([ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]), line.split()[5] ])

        ref_center = findgeocenter(ref_frame,firstidxsolvent)
        ref_frame_centered = centersystem(ref_frame[0:firstidxsolvent-1],ref_center,pbcbox[0:3], [0, 0, 0])
        print('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' read.')
        log_write('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' read.')
        print('INFO: Process uses # MB RAM: ' + RAMusage())
        log_write('INFO: Process uses # MB RAM: ' + RAMusage())



    with open(arcfile, 'r') as file, open('zzz_' + args.xyz + '.pdb', 'w') as out_file:
        for line in file:
            if line == first_line and frame == 0:
                numberatoms = int(line.split()[0])
                frame = frame + 1

            elif '.' in line.split()[0]:
                pbcbox = [float(i) for i in line.split()]

            elif line == first_line and (frame >= args.beginend[0]) and (frame <= args.beginend[1]):
                current_center = findgeocenter(current_frame,firstidxsolvent)
                translate = (np.array(current_center) - np.array(ref_center)).tolist()
                if args.atomsforwaters == []:
                    current_frame_centered = centersystem(current_frame[0:firstidxsolvent],current_center,pbcbox[0:3], translate)
                else:
                    current_frame_centered = centersystem(current_frame[0:],current_center,pbcbox[0:3], translate)
                current_frame_centered_with_waters = select_waters(current_frame_centered)
                current_frame_fit = fitsystem(current_frame_centered_with_waters,ref_frame_centered)
                printPDB(out_file, frame, current_frame_fit)
                print('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' printed.')
                log_write('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' printed.')
                print('INFO: Process uses # MB RAM: ' + RAMusage())
                log_write('INFO: Process uses # MB RAM: ' + RAMusage())


                numberatoms = int(line.split()[0])
                frame = frame + 1
                frames = frames + 1
                current_frame = []

            elif line == first_line:
                numberatoms = int(line.split()[0])
                print('INFO: ' + str(time.time()-tstart) + ' sec and avoiding frame ' + str(frame) + '.')
                log_write('INFO: ' + str(time.time()-tstart) + ' sec and avoiding frame ' + str(frame) + '.')

                frame = frame + 1
                frames = frames + 1
                current_frame = []

            else:
                try:
                    current_frame.append([ line.split()[0], line.split()[1], float(line.split()[2]), float(line.split()[3]), float(line.split()[4]), line.split()[5] ])
                except:
                    continue

            if (frame > args.beginend[1]):
                break

        if frame == total_frames and (frame <= args.beginend[1]):
                current_center = findgeocenter(current_frame,firstidxsolvent)
                translate = (np.array(current_center) - np.array(ref_center)).tolist()
                if args.atomsforwaters == []:
                    current_frame_centered = centersystem(current_frame[0:firstidxsolvent],current_center,pbcbox[0:3], translate)
                else:
                    current_frame_centered = centersystem(current_frame[0:],current_center,pbcbox[0:3], translate)
                current_frame_centered_with_waters = select_waters(current_frame_centered)
                current_frame_fit = fitsystem(current_frame_centered_with_waters,ref_frame_centered)
                printPDB(out_file, frame, current_frame_fit)
                print('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' printed.')
                log_write('INFO: ' + str(time.time()-tstart) + ' sec and frame ' + str(frame) + ' printed.')
                print('INFO: Process uses # MB RAM: ' + RAMusage())
                log_write('INFO: Process uses # MB RAM: ' + RAMusage())


                numberatoms = int(line.split()[0])
                frame = frame + 1
                frames = frames + 1
                current_frame = []


### MAIN ###
tstart = time.time()
log_write('INFO: Starting xyz2pdb of ' + args.arc + '.arc.')
make_dict()
log_write('INFO: Starting to reduce arc to protein only, aligning and printing as PDB:')
reduce_arc()

log_write('INFO: DONE!')
