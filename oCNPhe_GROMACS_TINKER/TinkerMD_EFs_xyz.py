# Author Jacek Kozuch, kozuch.jacek@gmail.com (based on Richard Bradshaw's and Stephen Frieds scripts)
# Script to calculate the Amoeba-based electric field at a desired site
# Requirements: Tinker with modified analyze function (added option U: copy-paste option D and turn debug = .false.)

# python TinkerMD_ProteinFields.py -st PYP_F92oCNF_md -t abs -a 1361 1362 -at C N -p 1.4150 1.1580 -f 2500
# python TinkerMD_ProteinFields.py -st PYP_F92oCNF_md -t abs -a 402 403 -at C N -p 1.4150 1.1580 -f 2500

# -arc PYP_F92oCNF_md -t abs -a 402 403 -at C N -p 1.4150 1.1580 -ax 396 -atx C


import os
import sys, argparse, glob
import re
import itertools
import numpy as np
from subprocess import call,Popen
from datetime import datetime
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit




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
    log_write('# Title: TinkerMD_ProteinFields                              #')
    log_write('# Description: Extraction of Electric Fields                 #')
    log_write('#                                                            #')
    log_write('# Copyright:            Copyright (c) Jacek Kozuch 2021-2022 #')
    log_write('#                                                            #')
    log_write('##############################################################')
    log_write('')
    log_write('')

    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--type",help="Type of analysis to run. 'abs' will calculate the absolute field exerted on a probe - manual referencing is required. 'env' will take the difference between the fields in the solvated and reference (probe only) phase at atoms 1 & 2. Defaults to 'env'. ",choices=['env','abs'],default='env')
    parser.add_argument("-a","--atoms",help="Atom numbers for analysis (Max 2)",nargs='+',type=int,required=True)
    parser.add_argument("-ax","--atomx",help="Atom number for definition of x-axis.",type=int,required=True)
    parser.add_argument("-at","--atomtypes",help="Atomtypes for analysis (Max 2)",nargs='+',type=str,default=['C','O'])
    parser.add_argument("-atx","--atomtypex",help="Atomtypes x axis definition",type=str,default='C')
    parser.add_argument("-p","--polarise",help="Polarisabilities for atoms 1 & 2. Defaults to carbonyl C and O, or nitrile C and N.",nargs=2,type=float)
    parser.add_argument("-sa","--soluteatoms",help="Number of atoms of your solute.",type=int)
    parser.add_argument("-arc","--traj",help="Filename of Solvated Trajectory (Tinker arc), or gasphase trajectory with '-t gas'. *.uind file of trajectory with same prefix needed it folder!",type=str,required=True)
    parser.add_argument("-rk","--refkey",help="Filename of Key file for calculating reference induced dipole (gas-phase calculation of solute)",type=str,default="ref.key")
    if len(sys.argv)==1:
       parser.print_help()
       sys.exit(1)  

    args = parser.parse_args()
    return args


def filenames(args):
    
    filename = {"prefix"    : args.traj,
                "arc"       : args.traj + '.arc',
                "uind"      : args.traj + '.uind',
                "arc_ref"   : args.refkey + '.arc',
                "uind_ref"  : args.refkey + '.uind',
                "key_ref"   : args.refkey + '.key',
                "xyz_1"     : 'xyz_'  + args.atomtypes[0] + str(args.atoms[0]) + '.txt',
                "xyz_2"     : 'xyz_'  + args.atomtypes[1] + str(args.atoms[1]) + '.txt',
                "xyz_3"     : 'xyz_'  + args.atomtypex    + str(args.atomx)    + '.txt',
                "uind_1"    : 'uind_' + args.atomtypes[0] + str(args.atoms[0]) + '.txt',
                "uind_2"    : 'uind_' + args.atomtypes[1] + str(args.atoms[1]) + '.txt',
                "uind_1_ref": 'uind_' + args.atomtypes[0] + str(args.atoms[0]) + '_ref.txt',
                "uind_2_ref": 'uind_' + args.atomtypes[1] + str(args.atoms[1]) + '_ref.txt'}
    
    return filename

def norm_vec(vec):
    
    vec_length = np.sqrt(np.dot(vec,vec))
    norm_vec = vec/vec_length
    
    return norm_vec

def gen_xyz_uind(args):

    files = filenames(args)    
    
    if os.path.isfile(files["xyz_1"]):
        log_write('INFO: ' + files["xyz_1"] + ' found.')    
    else:
        with open(files["arc"], 'r') as in_file, open(files["xyz_1"], 'w') as out_file:
            for line in in_file:
                if len(line.split()) > 1:
                    if str(line.split()[0]) == str(args.atoms[0]) and str(line.split()[1]) == str(args.atomtypes[0]):
                        out_file.write(line)
                        
    if os.path.isfile(files["xyz_2"]):
        log_write('INFO: ' + files["xyz_2"] + ' found.')    
    else:
        with open(files["arc"], 'r') as in_file, open(files["xyz_2"], 'w') as out_file:
            for line in in_file:
                if len(line.split()) > 1:
                    if str(line.split()[0]) == str(args.atoms[1]) and str(line.split()[1]) == str(args.atomtypes[1]):
                        out_file.write(line)
    
    if os.path.isfile(files["xyz_3"]):
        log_write('INFO: ' + files["xyz_3"] + ' found.')    
    else:
        with open(files["arc"], 'r') as in_file, open(files["xyz_3"], 'w') as out_file:
            for line in in_file:
                if len(line.split()) > 1:
                    if str(line.split()[0]) == str(args.atomx) and str(line.split()[1]) == str(args.atomtypex):
                        out_file.write(line)    

    if os.path.isfile(files["uind_1"]):
        log_write('INFO: ' + files["uind_1"] + ' found.')    
    else:
        with open(files["uind"], 'r') as in_file, open(files["uind_1"], 'w') as out_file:
            for line in in_file:
                if len(line.split()) > 1:
                    if str(line.split()[0]) == str(args.atoms[0]) and str(line.split()[1]) == str(args.atomtypes[0]):
                        out_file.write(line)
                        
    if os.path.isfile(files["uind_2"]):
        log_write('INFO: ' + files["uind_2"] + ' found.')    
    else:
        with open(files["uind"], 'r') as in_file, open(files["uind_2"], 'w') as out_file:
            for line in in_file:
                if len(line.split()) > 1:
                    if str(line.split()[0]) == str(args.atoms[1]) and str(line.split()[1]) == str(args.atomtypes[1]):
                        out_file.write(line)


def electric_fields(args):
    
    files = filenames(args)
    frame = 0
    pol_1 = args.polarise[0]
    pol_2 = args.polarise[1]
       
    if args.type == 'abs':
        with open(files["xyz_1"], 'r') as xyz_1_file,\
             open(files["xyz_2"], 'r') as xyz_2_file,\
             open(files["xyz_3"], 'r') as xyz_3_file,\
             open(files["uind_1"], 'r') as uind_1_file,\
             open(files["uind_2"], 'r') as uind_2_file,\
             open('zzz_RESULT.txt', 'w') as out_file:
                 
            for xyz_1, xyz_2, xyz_3, uind_1, uind_2 in zip(xyz_1_file,xyz_2_file,xyz_3_file,uind_1_file,uind_2_file):
                frame += 1
                x_1 = np.array([float(x) for x in xyz_1.split()[2:5]])
                x_2 = np.array([float(x) for x in xyz_2.split()[2:5]])
                x_3 = np.array([float(x) for x in xyz_3.split()[2:5]])
                u_1 = np.array([float(u) for u in uind_1.split()[2:5]])
                u_2 = np.array([float(u) for u in uind_2.split()[2:5]])
                               
                z_axis = norm_vec(x_2 - x_1)
                x_axis = x_3 - x_1
                y_axis = norm_vec(np.cross(z_axis, x_axis))
                x_axis = norm_vec(np.cross(y_axis, z_axis))                
                
                
                ef_1 = u_1/pol_1*299.79
                ef_2 = u_2/pol_2*299.79
                ef_avg  = (ef_1 + ef_2)/2
                ef_drop = ef_2 - ef_1
                
                ef_x = np.round(np.dot(ef_avg,x_axis),3)
                ef_y = np.round(np.dot(ef_avg,y_axis),3)
                ef_z = np.round(np.dot(ef_avg,z_axis),3)
                
                ef_drop_x = np.round(np.dot(ef_drop,x_axis),3)
                ef_drop_y = np.round(np.dot(ef_drop,y_axis),3)
                ef_drop_z = np.round(np.dot(ef_drop,z_axis),3)
                
                ef_off_axis = np.round(np.sqrt(ef_x**2 + ef_y**2),3)
                ef_drop_off_axis = np.round(np.sqrt(ef_drop_x**2 + ef_drop_y**2),3)
                
                out_file.write(str(frame).rjust(6)+\
                               str(ef_x).rjust(12)+\
                               str(ef_y).rjust(12)+\
                               str(ef_z).rjust(12)+\
                               str(ef_off_axis).rjust(12)+\
                               str(ef_drop_x).rjust(12)+\
                               str(ef_drop_y).rjust(12)+\
                               str(ef_drop_z).rjust(12)+\
                               str(ef_drop_off_axis).rjust(12)+'\n')




args = parse()
gen_xyz_uind(args)
electric_fields(args)





#################
### Rest ###
#################


### Check if numbers in line and not header ###
def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

### Check if polarizabilities are correct according to amoeba09.prm ##
def check_pol():
    if args.polarise is not None:
        pol_1 = args.polarise[0]
        pol_2 = args.polarise[1]
        log_write('Setting polarizabilities to specified values for %s - %s bond: %1.4f, %1.4f.' %(args.atomtypes[0],args.atomtypes[1],pol_1,pol_2))

    # Add later to just check in amoebabio18.prm for polarizabilities

    return pol_1, pol_2

### Strip solvent from trajectory and calculate induced dipoles for solute only ### 
def prep_ref():
    log_write("Reading trajectory %s as the solvated system and creating the gasphase system to obtain self-polarization (using %s)." % (args.traj, args.refkey))
    arcfile_sol = args.traj
    arcname_sol = args.traj.split('.')[0]
    uindfile_sol = args.traj.split('.')[0] + '.uind'
    keyfile_ref = args.refkey
    arcfile_ref = args.refkey.split('.')[0] + '.arc'
    uindfile_ref = args.refkey.split('.')[0] + '.uind'
    # Check here if ref.arc already exists - in case you are just repeating the run
    if os.path.isfile(uindfile_sol) is False:
        log_write('%s is missing. Cannot continue without calculating induced dipoles of the solvated system.' % uindfile_sol)
        quit()
    if os.path.isfile(keyfile_ref):
        log_write('Reading %s file for determination of induced dipoles for bare solute.' % keyfile_ref)
    else:
        log_write('Error: Key file for determination of gasphase induced dipoles is missing.')
        log_write('Copy *.key into folder and restart.')
        quit()
    with open(arcfile_sol,'r') as f:
        total_atoms = int(f.readline()) 
    if isnumber(total_atoms) is False:
        log_write('Error: Cannot determine total number of atoms in box.')
        quit()
    if os.path.isfile(arcfile_ref):
            if (len(open(arcfile_ref,'r').readlines()) >= ((args.soluteatoms+1)*args.frames-1)) and (len(open(arcfile_ref,'r').readlines()) <= ((args.soluteatoms+1)*args.frames+1)):
                log_write('Warning: %s (file without solvent) already in folder. Continue without stripping solvent from MD %s file.' % (arcfile_ref, arcfile_sol))
    if os.path.isfile(arcfile_ref) is False:
        log_write('Building refphase %s file.' % arcfile_ref)
        for num, line in enumerate(open(arcfile_sol,'r')):
            if isnumber(line):
                if int(line) == total_atoms: 
                    num_rem = num
                    with open(arcfile_ref,'a') as f: f.write('    %d\n' % args.soluteatoms)
            if (num > num_rem+1) and ((num - num_rem) < (args.soluteatoms+2)):
                with open(arcfile_ref,'a') as f: f.write(line)
    if os.path.isfile(uindfile_ref):
        if (len(open(uindfile_ref,'r').readlines()) >= ((args.soluteatoms+1)*args.frames-1)) and (len(open(uindfile_ref,'r').readlines()) <= ((args.soluteatoms+1)*args.frames+1)):
            log_write('Warning: %s file for gasphase already exists in folder. Continuing without building new *.uind file.' % uindfile_ref)
        else:
            log_write('Uind file in folder is corrupted. Computing new gas system %s file.' % uindfile_ref)
            call('analyze_uind %s U' % arcfile_ref, shell=True, stdout=None, stderr=None)
    else:
        log_write('Computing ref system %s file.' % uindfile_ref)
        call('analyze_uind %s U' % arcfile_ref, shell=True, stdout=None, stderr=None)




def count_lines(files):

    lines = []
    for file in files:
        with open(file, 'r') as f:
            lines.append(len(f.readlines()))
    if all(line == lines[0] for line in lines):
        return lines[0]
    else: 
        log_write('INFO: Something is wrong with length of the xyz and uind file - not same length')
        return min(lines)

### Calculated electric fields on specified bond ####
def calc_EF(inp_1,inp_2):

    global args, arcfile, uindfile, uindfile_ref, xyz_1, xyz_2, uind_1, uind_2, uind_1_ref, uind_2_ref, frames

    [pol_1,pol_2] = [inp_1,inp_2]

    num_lines = count_lines([xyz_1, xyz_2, uind_1, uind_2])
    frames    = num_lines
    print(num_lines)
 
    # Get the relevant coordinates and induced dipoles from trajectories
    x1 = open(xyz_1 , 'r').readlines()
    x2 = open(xyz_2 , 'r').readlines()
    u1 = open(uind_1 , 'r').readlines()
    u2 = open(uind_2 , 'r').readlines()
    if args.type == 'env':    
        u1_ref = open(uind_1_ref , 'r').readlines()
        u2_ref = open(uind_2_ref , 'r').readlines()
    elif args.type == 'abs':
        u1_ref = u1
        u2_ref = u2

    frame_num = np.arange(1,num_lines+1)

    # Calculate the the fields
    log_write('Calculating fields on specified bond.')
    fields = open('zzz_EFs_on_%s%s.txt' %(args.atomtypes[0],args.atomtypes[1])  , 'w')
    for line_frame, line_x1, line_x2, line_u1, line_u2, line_u1_ref, line_u2_ref in zip(frame_num, x1, x2, u1, u2, u1_ref, u2_ref):
#        if isnumber( line_x1.split()[0] ):
                #Get the coordinates - nomenclature of variables for CO bond:
                #Calculate the XY bond length and XY unit vector
                x_C = np.array([float(x) for x in line_x1.split()[2:5]])
                x_O = np.array([float(x) for x in line_x2.split()[2:5]])
                COvec = x_O - x_C
                COlen = np.sqrt((COvec**2).sum())
                COunitvec = COvec/ COlen

                #Get the C and O induced dipoles when solvated and not
                u_C = np.array([float(x) for x in line_u1.split()[2:5]])
                u_O = np.array([float(x) for x in line_u2.split()[2:5]])
                if args.type == 'env':   
                   u_C_ref = np.array([float(x) for x in line_u1_ref.split()[2:5]])
                   u_O_ref = np.array([float(x) for x in line_u2_ref.split()[2:5]])
                elif args.type == 'abs':
                   u_C_ref = np.zeros(3)
                   u_O_ref = np.zeros(3)
                #Calculate the electric field vectors when solvated and not 
                EF_C = u_C/pol_1*299.79
                EF_O = u_O/pol_2*299.79
                EF_C_ref = u_C_ref/pol_1*299.79
                EF_O_ref = u_O_ref/pol_2*299.79

                #Calculate the induced dipole projection onto CO
                EFproj_C = np.dot( EF_C, COunitvec ) - np.dot( EF_C_ref, COunitvec )
                EFproj_O = np.dot( EF_O, COunitvec ) - np.dot( EF_O_ref, COunitvec )
                EFproj_C_ref = np.dot( EF_C_ref, COunitvec )
                EFproj_O_ref = np.dot( EF_O_ref, COunitvec )

                #Calculate the electric field and convert
                EFproj = (EFproj_C + EFproj_O)/2
                EFproj_ref = (EFproj_C_ref + EFproj_O_ref)/2
                EFprojDrop = (EFproj_O-EFproj_C)
                EFprojDrop_ref = (EFproj_O_ref-EFproj_C_ref)

                #Print to file
                if args.type == 'env':
                    fieldInfo = str( line_frame ) + '\t' +str( EFproj_C ) + '\t' + str( EFproj_O ) + '\t' + str( EFproj )+ '\t' + str( EFprojDrop ) + '\t' + str( COlen ) + '\t' + str( EFproj_ref )+ '\t' + str( EFprojDrop_ref ) + '\n'
                    fields.write ( fieldInfo )
                elif args.type == 'abs':    
                    fieldInfo = str( line_frame ) + '\t' +str( EFproj_C ) + '\t' + str( EFproj_O ) + '\t' + str( EFproj )+ '\t' + str( EFprojDrop ) + '\t' + str( COlen ) + '\n'
                    fields.write ( fieldInfo )         
    fields.close()

### Do statistics on the results ###
def get_statistics():
    global frames
    fields = np.loadtxt(open('zzz_EFs_on_%s%s.txt' %(args.atomtypes[0],args.atomtypes[1])  , 'r'))

    bins_calc = int(int(frames)/10 + 1)
    
    f_bond = np.histogram(fields[:,3], bins = bins_calc, range = [-500, 501], density = True)
    f_drop = np.histogram(fields[:,3], bins = bins_calc, range = [-500, 501], density = True)

    x = f_bond[1][0:-1] + 10000/int(int(frames)/2)
    y = f_bond[0]

#    f_max = len(fields[:,3])
    f_max = 1
    f_avg = np.mean(fields[:,3])
    f_std = np.std(fields[:,3])
    f_med = np.median(fields[:,3])
    f_asy = f_avg - f_med
      
    def gauss_1(x, a_1, x_1, s_1):
        return a_1/(s_1*np.sqrt(2*np.pi))*np.exp(-(x-x_1)**2/(2*s_1**2))
    def gauss_2(x, a_1, x_1, s_1, a_2, x_2, s_2):
        return a_1/(s_1*np.sqrt(2*np.pi))*np.exp(-(x-x_1)**2/(2*s_1**2)) + a_2/(s_2*np.sqrt(2*np.pi))*np.exp(-(x-x_2)**2/(2*s_2**2))
    def gauss_3(x, a_1, x_1, s_1, a_2, x_2, s_2, a_3, x_3, s_3):
        return a_1/(s_1*np.sqrt(2*np.pi))*np.exp(-(x-x_1)**2/(2*s_1**2)) + a_2/(s_2*np.sqrt(2*np.pi))*np.exp(-(x-x_2)**2/(2*s_2**2)) + a_3/(s_3*np.sqrt(2*np.pi))*np.exp(-(x-x_3)**2/(2*s_3**2))
    def write_result(popt, perr, R2):
        log_write('INFO: Gaussian fit results on %s%s:' %(args.atomtypes[0],args.atomtypes[1]))
        log_write('INFO: Fitting parameters: ' + " ".join(str(round(element,2)) for element in popt))
        log_write('INFO: Standard errors   : ' + " ".join(str(round(element,2)) for element in perr))
        log_write('INFO: R2                : ' + str(R2))
    def fit_gauss(which, x, y, p0):
        if which == 'gauss_1':
            popt, pcov = curve_fit(gauss_1, x, y, p0)
            fit_hist = gauss_1(x, popt[0], popt[1], popt[2])
        elif which == 'gauss_2':
            popt, pcov = curve_fit(gauss_2, x, y, p0)
            fit_hist = gauss_2(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
        elif which == 'gauss_3':
            popt, pcov = curve_fit(gauss_3, x, y, p0)
            fit_hist = gauss_3(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
        perr = np.sqrt(np.diag(pcov))
        R2 = 1 - sum((y - fit_hist)**2)/sum((y - np.mean(y))**2)
        resid = y - fit_hist
        return popt, perr, R2, resid

    R2 = 0
    if R2 < 0.99:
        [popt, perr, R2, resid] = fit_gauss('gauss_1', x, y, [f_max, f_avg, f_std])
        write_result(popt, perr, R2)

    if R2 < 0.99:
        f_max = [popt[0]/2, popt[0]/2]
        f_avg = [popt[1] + f_asy, popt[1] - f_asy]
        f_std = [popt[2], popt[2]]
        [popt, perr, R2, resid] = fit_gauss('gauss_2', x, y, [f_max[0], f_avg[0], f_std[0], f_max[1], f_avg[1], f_std[1]])
        write_result(popt, perr, R2)

    if R2 < 0.99:
        f_max = [popt[0]/1.5, (popt[0]+popt[3])/3, popt[3]/1.5]
        f_avg = [popt[1], (popt[1]+popt[4])/2, popt[4]]
        f_std = [popt[2], (popt[2]+popt[5])/2, popt[5]]
        [popt, perr, R2, resid] = fit_gauss('gauss_3', x, y, [f_max[0], f_avg[0], f_std[0], f_max[1], f_avg[1], f_std[1], f_max[2], f_avg[2], f_std[2]])
        write_result(popt, perr, R2)

#    log_write(('Average electric field on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],field_avg,field_std))
#    log_write(('Average electric field drop on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],fielddrop_avg,fielddrop_std))
    
#    if args.type == 'env':
#        selffield_avg = np.mean(fields[:,6])
#        selffield_std = np.std(fields[:,6])
#        selffielddrop_avg = np.mean(fields[:,7])
#        selffielddrop_std = np.std(fields[:,7])
#        log_write(('Average self-electric field on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],selffield_avg,selffield_std))
#        log_write(('Average self-electric field drop on %s%s: %3.2f +/- %3.2f') %(args.atomtypes[0],args.atomtypes[1],selffielddrop_avg,selffielddrop_std))

# #### MAIN BELOW HERE ####
# [pol_1, pol_2] = check_pol()
# if args.type == 'env':
#     prep_ref()
# write_xyz_uind()
# calc_EF(pol_1,pol_2)
# get_statistics()    
# log_write('Electric Fields determined successfully.')
# #########################
    
    
