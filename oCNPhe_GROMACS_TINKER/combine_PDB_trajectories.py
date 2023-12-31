import os
import time
import numpy as np
import psutil
from datetime import datetime
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as Rot
files = []

############################################
# Filenames
files.append('zzz_PYP_F28oCNF_A_md_1.pdb')
files.append('zzz_PYP_F28oCNF_A_md_2.pdb')
files.append('zzz_PYP_F28oCNF_B_md_1.pdb')
files.append('zzz_PYP_F28oCNF_B_md_2.pdb')
out_filename = 'zzz_PYP_F28oCNF_md.pdb'

align_by = ['N', 'CA', 'C', 'O']
# align_by = ['CA']

# Debug
debug = True
to_debug = 2210

#############################################


global frame_for_fit
global p0_rottrans


start = time.time()

def write_log(string):
    if isinstance(string, list):
        string = ' '.join(str(x) for x in string)
    now    = datetime.now()
    dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
    print(dt_str + ' ' + string)
    with open('zzz_log.txt','a') as log:
        log.write(dt_str + ' ' + string + '\n')
        
def RAMusage():
    process = psutil.Process(os.getpid())
    RAM = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    return write_log('INFO: Using %s MB RAM.' %(str(RAM)))


def flatten(inp):
    new_list = []
    for i in inp:
        for j in i:
            new_list.append(j)
    return new_list

def get_rotmat(in_xyz,in_xyz_ref):
    global tstart  
    global frame_for_fit
    global p0_rottrans

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

    # print('NEXT ROUND')
    # print(p0_rottrans[frame_for_fit])
    p0 = p0_rottrans[frame_for_fit]
    popt, pcov = curve_fit(rotate, xyz, flatten(xyz_ref), p0)
    p0_rottrans[frame_for_fit] = popt

    RotMatrix = RotMat(popt[0], popt[1], popt[2])
    TransVec  = [popt[3], popt[4], popt[5]]

    return RotMatrix, TransVec


def fitsystem(in_current,in_first):

    [RotMatrix, TransVec] = get_rotmat(in_current, in_first)

    current = (np.array(in_current).transpose()).astype(float)
    new_current = np.matmul(RotMatrix, current).transpose() - np.array(TransVec)

    return new_current

def lin2xyz_coords(coords):
    
    new_coords = []
    
    for i in range(0, int(len(coords)/3)):
        new_coords.append([coords[3*i], coords[3*i+1], coords[3*i+2]])
        
    return new_coords
    
    
def xyz2lin_coords(coords):
    
    new_coords = []
    
    for atom in coords:
        new_coords.append(atom[0])
        new_coords.append(atom[1])
        new_coords.append(atom[2])
        
    return new_coords


def RotMat(alpha,beta,gamma):

     Rz_alpha = np.array([[np.cos(alpha), -np.sin(alpha), 0            ], [ np.sin(alpha),  np.cos(alpha),  0            ], [ 0           , 0            , 1            ]])
     Ry_beta  = np.array([[np.cos(beta) ,  0            , np.sin(beta) ], [ 0            ,  1            ,  0            ], [-np.sin(beta), 0            , np.cos(beta) ]])
     Rx_gamma = np.array([[1            ,  0            , 0            ], [ 0            ,  np.cos(gamma), -np.sin(gamma)], [ 0           , np.sin(gamma), np.cos(gamma)]])
     
     return np.matmul(Rz_alpha,np.matmul(Ry_beta,Rx_gamma))    

def rotate(xyz, alpha, beta, gamma, t_x, t_y, t_z):  
    
    # print(xyz)
    # print(alpha)
    # print(beta)
    # print(gamma)
    # print(t_x)
    # print(t_y)
    # print(t_z)

    R = RotMat(alpha,beta,gamma)
    return np.matmul(R,np.array(xyz)) - np.array([t_x, t_y, t_z])

def avg_struct(coord_list, *p):
    global frame_for_fit
    global p0_rottrans
       
    if len(p) == 1:
        p = p[0]
               
    frames = int(len(coord_list)/len(p))
          
    ref_coords = lin2xyz_coords(p)
      
    if p0_rottrans == []:
        for i in range(0, frames):
            p0_rottrans.append([0, 0, 0, 0, 0, 0])
    
    new_coords = []
    for i in range(0, frames):
        frame_for_fit = i
        temp_coords = []
        for j in range(0, len(p)):
            temp_coords.append(coord_list[len(p)*i+j])
        temp_coords = lin2xyz_coords(temp_coords)
        new_coords.append(xyz2lin_coords(fitsystem(temp_coords, ref_coords)))
        
    return np.mean(new_coords,axis=0)
    


#### MAIN ####

def which_atoms(align_by, info):

    atoms_all = []
    atoms_PRO = []
    
    for atom in align_by:
        if atom == 'N':
            atoms_all.append(atom)
        else:
            atoms_all.append(atom)
            atoms_PRO.append(atom) 
    if info == 'PRO':
        return atoms_PRO
    else:
        return atoms_all
    


remember = []
coords = []
atoms  = []
frames = 0
for file in files:
    with open(file, 'r') as in_file:
        previous_residue = 0
        collect = False
        for line in in_file:
            if 'MODEL' in line:
                # print(line)
                if debug == True and frames >= to_debug:
                    # frames -= 1
                    break
                frames += 1
                previous_residue = 0

            elif 'ENDMDL' in line:
                write_log('ENDMDL found')
            elif line.split()[3] == 'SOL':
                collect = False
                continue
            elif len(line.split()) > 5:
                if previous_residue < int(line.split()[5]):
                    if line.split()[3] == 'PRO':
                        atoms = which_atoms(align_by, 'PRO')
                    else:
                        atoms = which_atoms(align_by, '')
                    collect = True
                    previous_residue = int(line.split()[5])
                                       
            if collect == True:
                if len(line.split()) > 5:
                    if line.split()[2] in atoms:
                        coords.append(float(line.split()[6]))
                        coords.append(float(line.split()[7]))
                        coords.append(float(line.split()[8]))
                        atoms.remove(line.split()[2])
                        if atoms == []:
                            collect = False
                        remember.append([frames, line.split()[5], line.split()[3], line.split()[2]])
    if debug == True and frames >= to_debug:
        # frames -= 1
        break                    
                   
num_coords = int(len(coords)/frames)
avg_coords = coords[0:num_coords]
        
end = time.time()

write_log('Reading data took ' + str(round((end-start)/60,2)) + ' minutes.')
RAMusage()

p0_rottrans = []
    
write_log('Aligning frames until average converged.')

delta_rmsd = 1
rmsd = 1
while abs(delta_rmsd) > 0.000001:    
    
    new_avg_coords = avg_struct(coords, avg_coords)
    resid = avg_coords - new_avg_coords
    new_rmsd  = np.sqrt(np.mean(resid**2))
    delta_rmsd = rmsd - new_rmsd

    write_log('####')
    write_log('RSMD ' + str(rmsd))
    write_log('New RSMD ' + str(new_rmsd))
    write_log('Delta RSMD ' + str(delta_rmsd))
    RAMusage()
    
    rmsd = new_rmsd
    avg_coords = new_avg_coords
    
    
    
write_log('Writing new aligned PDB file.')
# print(p0_rottrans)

frames = 0
coords = []
with open(out_filename, 'w') as out_file:
    for file in files:
        with open(file, 'r') as in_file:
            for line in in_file:
                if 'MODEL' in line:
                    print(line)
                    if debug == True and frames >= to_debug:
                        # frames -= 1
                        break
                    frames += 1
                    out_file.write(line)
                elif 'ENDMDL' in line:
                    # out_file.write('ENDMDLMODEL       1 \n')
                    write_log('ENDMDL found')
                elif len(line.split()) > 5:
                    coords = np.array([float(x) for x in line.split()[6:9]])
                    p = p0_rottrans[frames-1]
                    new_coords = rotate(coords, p[0],p[1],p[2],p[3],p[4],p[5])
                    string = line[0:30] +\
                             str(round(new_coords[0],3)).rjust(8) +\
                             str(round(new_coords[1],3)).rjust(8) +\
                             str(round(new_coords[2],3)).rjust(8) +\
                             line[55:]
                    out_file.write(string)
        if debug == True and frames >= to_debug:
            # frames -= 1
            break           
    out_file.write('ENDMDL')


end = time.time()

write_log('FINISHED!')