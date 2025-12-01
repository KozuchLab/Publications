import json
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
from subprocess import run
import itertools
import matplotlib.pyplot as plt

def openjson(input_name):
    with open(input_name, 'r') as file:
        return json.load(file) 

def isint(string):
    try:
        first = int(string)
        return True
    except:
        return False
    
def flatlist(array):
    new_array = []
    for xs in array:
        for x in xs:
            new_array.append(x)
    return new_array

def addzeros(number):
    number = str(number)
    
    if len(number) == 1:
        number = '00' + number
    elif len(number) == 2:
        number = '0' + number
        
    return number

def delete_annoying_files(xyzfile,N):
    for i in range(1, N):
        if i < 10:
            string = '00'+str(i)
        elif i < 100:
            string = '0'+str(i)
        elif i >= 100:
            string = str(i)
        try:
            os.remove(xyzfile.replace('.xyz', '.'+string))
        except:
            pass
        try:
            os.remove(xyzfile+'_'+string)
        except:
            pass
        # run('rm *.%s* >./prompt_errors.txt' %(string),shell=True)

def number_of(input_xyz):
    with open(input_xyz, 'r') as in_file:
        for line in in_file:
            number_of_atoms = int(line)
            number_of_NMs = 3*number_of_atoms - 6
            break    

    return number_of_atoms, number_of_NMs

def xyz_linear(QM_xyz):

    np_QM_xyz = np.array(QM_xyz).flatten()

    return np_QM_xyz.tolist()


def xyz_Nby3(QM_xyz):

    length = int(len(QM_xyz)/3)

    new_QM_xyz = []
    for i in range(0, length):
        new_QM_xyz.append([QM_xyz[3*i], QM_xyz[3*i+1], QM_xyz[3*i+2]])
        
    return new_QM_xyz

def rotmat(old,new):
    np_old = np.array(old)
    np_new = np.array(new)
    
    [rot, rmsdrot] = R.align_vectors(np_old,np_new)
  
    return rot.as_matrix()

def rotate(old,rotmat):
    np_old = np.array(old)
    np_new = np.matmul(np_old,rotmat)
    return np_new.tolist()

def NMA_watch(all_gau_data, all_tink_data):
    with open('zzz_NMA_watch.txt', 'a') as out_file:
        string = 'qm_v'.rjust(8) + ' | ' + 'mm_v'.rjust(8) + ' | ' + 'qm_I'.rjust(8) + ' | ' + 'mm_I'.rjust(8) + '\n'
        out_file.write(string)
        for qm_v, mm_v, qm_I, mm_I, qm_dd, mm_dd in zip(all_gau_data['freqs'],\
                                                        all_tink_data['freqs'],\
                                                        all_gau_data['irints'],\
                                                        all_tink_data['irints'],\
                                                        all_gau_data['norm_tdms'],\
                                                        all_tink_data['norm_tdms']):
            try:
                str1 = str(round(qm_v)).rjust(8) + ' | '
                str2 = str(round(mm_v)).rjust(8) + ' | ' 
                str3 = str(round(qm_I)).rjust(8) + ' | ' 
                str4 = str(round(mm_I)).rjust(8) + ' | ' 
                str5 = ' '.join(str(x).rjust(8) for x in np.round(qm_dd,2)) + ' | '
                str6 = ' '.join(str(x).rjust(8) for x in np.round(mm_dd,2))
                string = str1 + str2 + str3 + str4 + str5 + str6 +  '\n'
            except:
                string = 'NANs found'
            out_file.write(string)
        out_file.flush() 
        
def NMA_plot(ini):
    
    def gaussian(x, sigma, x_0):
        y = (1/sigma/np.sqrt(2*np.pi))*np.exp(-(x-x_0)**2/2/sigma**2)     
        return y

    def lorentzian(x, sigma, x_0):
        y = (sigma/np.pi)*(1/((x - x_0)**2 + sigma**2))
        return y     

    def pseudovoigt(x, frac, sigma, x_0):
        sigma_G = sigma
        sigma_L = sigma * np.sqrt(2*np.log(2))
        y = frac * gaussian(x, sigma_G, x_0) + (1 - frac) * lorentzian(x, sigma_L, x_0)
        return y  
    
    def FWHM_to_sigma(FWHM):
        return FWHM/(2*np.sqrt(2*np.log(2)))
    
    FWHM = 12
    frac = 0.8
    wn_axis = np.arange(0,4001,1)
    sigma = FWHM_to_sigma(FWHM)
    qm_specs = []
    mm_specs = []
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        
        qm_temp = 0
        for A,w in zip(all_gau_data['irints'],all_gau_data['freqs']):
            qm_temp = qm_temp + A*pseudovoigt(wn_axis, frac, sigma, w)
        qm_specs.append(qm_temp)
        
        mm_temp = 0
        for A,w in zip(all_tink_data['irints'],all_tink_data['freqs']):
            mm_temp = mm_temp + A*pseudovoigt(wn_axis, frac, sigma, w)
        mm_specs.append(mm_temp)
        
    max_int = 1.1*max([max(max(x) for x in qm_specs), max(max(x) for x in mm_specs)])
    
    
    plt.figure(figsize=(8,10))
    for i in range(len(qm_specs)):
        qm_spec = np.subtract(qm_specs[i],i*max_int)
        mm_spec = np.subtract(mm_specs[i],i*max_int)
        which = ini['filenamejsons'][i]
        plt.plot(wn_axis,qm_spec, color = 'black', label = which.replace('.json',''))
        plt.plot(wn_axis,mm_spec, color = 'red')
    plt.xlim(ini['plot_range'][0], ini['plot_range'][1])
    plt.xlabel('Wavenumber / cm-1')
    plt.ylabel('Intensity / km/mol')
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0)
    plt.savefig('zzz_NMA_watch.png',bbox_inches="tight")
    plt.clf()
        
def write_pdb(filenames):
    callstring = 'xyzpdb %s -k %s' %(filenames['tink_min_aligned_xyz'], filenames['tink_temp_key'])
    run(callstring, shell=True)
    
def write_pdb_gau(filenames):
    callstring = 'xyzpdb %s -k %s' %(filenames['tink_min_xyz'], filenames['tink_temp_key'])
    run(callstring, shell=True)

def same_atom(atom_types):
    with open('filenames.json','r') as in_file:
        filenames = json.load(in_file)
        
    type_to_atom = {}
    with open(filenames['input_xyz'], 'r') as in_file:
        for line in in_file:
            if len(line.split()) > 1:
                type_to_atom[line.split()[5]] = line.split()[1]
            
    amoeba_type_1 = atom_types.split()[1]
    amoeba_type_2 = atom_types.split()[2]
    amoeba_type_3 = atom_types.split()[3]
    
    amoeba_atom_1 =  type_to_atom[amoeba_type_1]
    amoeba_atom_2 =  type_to_atom[amoeba_type_2]
    amoeba_atom_3 =  type_to_atom[amoeba_type_3]
    
    # print('XXXXXXXX')
    # print('Trying:' + atom_types)
    # print(amoeba_atom_1 + ' ' + amoeba_atom_2 + ' ' + amoeba_atom_3)
    if amoeba_atom_1 == amoeba_atom_3:
        # print('True')
        return True
    elif amoeba_atom_1 in  ['C','H','S'] and amoeba_atom_2 in  ['C'] and  amoeba_atom_3 in  ['C','H','S']:
        # print('True')
        return True
    else:
        # print('False')
        return False              
    
def int_coords(filenames,cart_coords,number_of_atoms):
        
    t2a = []
    all_bond_idxs = []
    all_angle_idxs = []
    all_torsion_idxs = []

    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            atomtype = int(line.split()[5])
            atomidx = int(int(line.split()[0])-1)
            t2a.append([atomtype, atomidx])

    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            for connect in line.split()[6:]:
                atom1 = int(int(line.split()[0])-1)
                atom2 = int(int(connect)-1)
                all_bond_idxs.append([atom1, atom2])

    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            for connect in line.split()[6:]:
                atom2 = int(int(line.split()[0])-1)
                atom3 = int(int(connect)-1)
                for bond in all_bond_idxs:
                    if atom2 == bond[1]:
                        if atom3 != bond[0]:
                            all_angle_idxs.append([bond[0],bond[1],atom3])   
                            
    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            for connect in line.split()[6:]:
                atom1 = int(int(connect)-1)
                atom2 = int(int(line.split()[0])-1)
                for angle in all_angle_idxs:
                    if atom2 == angle[0]:
                        if atom1 != angle[1]:
                            all_torsion_idxs.append([atom1,angle[0],angle[1],angle[2]])  
    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            for connect in line.split()[6:]:
                atom3 = int(int(line.split()[0])-1)
                atom4 = int(int(connect)-1)
                for angle in all_angle_idxs:
                    if atom3 == angle[1]:
                        if atom4 != angle[2]:
                            all_torsion_idxs.append([angle[0],angle[1],angle[2],atom4]) 

    set_bond_idxs = []
    for bond in all_bond_idxs:
        if bond[0] < bond[1]:
            new_bond = [bond[0], bond[1]]
            if new_bond not in set_bond_idxs:
                set_bond_idxs.append(new_bond)
        elif bond[1] < bond[0]:
            new_bond = [bond[1], bond[0]]
            if new_bond not in set_bond_idxs:
                set_bond_idxs.append(new_bond)
                
    set_angle_idxs = []
    for angle in all_angle_idxs:
        if angle[0] < angle[2]:
            new_angle = [angle[0], angle[1], angle[2]]
            if new_angle not in set_angle_idxs:
                set_angle_idxs.append(new_angle)
        elif angle[2] < angle[0]:
            new_angle = [angle[2], angle[1], angle[0]]
            if new_angle not in set_angle_idxs:
                set_angle_idxs.append(new_angle)
                
    set_torsion_idxs = []
    for torsion in all_torsion_idxs:
        if torsion[1] < torsion[2]:
            new_torsion = [torsion[0], torsion[1], torsion[2], torsion[3]]
            if new_torsion not in set_torsion_idxs:
                set_torsion_idxs.append(new_torsion)
        elif torsion[1] < torsion[2]:
            new_torsion = [torsion[3], torsion[2], torsion[1], torsion[0]]
            if new_torsion not in set_torsion_idxs:
                set_torsion_idxs.append(new_torsion)
                
    bonds = []
    angles = []

    for bond_idxs in set_bond_idxs:
        idx1 = bond_idxs[0]
        idx2 = bond_idxs[1]
        xyz1 = cart_coords[idx1]
        xyz2 = cart_coords[idx2]
        bond_vec = np.subtract(xyz1,xyz2)
        bond_len = np.sqrt(np.dot(bond_vec,bond_vec))
        bonds.append(bond_len)

    for angle_idxs in set_angle_idxs:
        idx1 = angle_idxs[0]
        idx2 = angle_idxs[1]
        idx3 = angle_idxs[2]
        xyz1 = cart_coords[idx1]
        xyz2 = cart_coords[idx2]
        xyz3 = cart_coords[idx3]
        bond_vec12 = np.subtract(xyz1,xyz2)
        bond_vec32 = np.subtract(xyz3,xyz2)
        bond_len12 = np.sqrt(np.sum(bond_vec12**2))
        bond_len32 = np.sqrt(np.sum(bond_vec32**2))
        bond_vec12 = bond_vec12/bond_len12
        bond_vec32 = bond_vec32/bond_len32
        angle = np.arccos(np.dot(bond_vec12,bond_vec32))*180/np.pi
        angles.append(angle)
        
    torsions = []
    for tor_idxs in set_torsion_idxs:
        idx1 = tor_idxs[0]
        idx2 = tor_idxs[1]
        idx3 = tor_idxs[2]
        idx4 = tor_idxs[3]
        xyz1 = cart_coords[idx1]
        xyz2 = cart_coords[idx2]
        xyz3 = cart_coords[idx3]
        xyz4 = cart_coords[idx4]
        bond_vec12 = np.subtract(xyz1,xyz2)
        bond_vec23 = np.subtract(xyz2,xyz3)
        bond_vec34 = np.subtract(xyz3,xyz4)
        bond_len12 = np.sqrt(np.sum(bond_vec12**2))
        bond_len23 = np.sqrt(np.sum(bond_vec23**2))
        bond_len34 = np.sqrt(np.sum(bond_vec34**2))
        bond21 = bond_vec12/bond_len12
        bond32 = bond_vec23/bond_len23
        bond43 = bond_vec34/bond_len34
        plane123 = np.cross(bond21,bond32)
        plane234 = np.cross(bond32,bond43)
        plane123_len = np.sqrt(np.dot(plane123,plane123))
        plane234_len = np.sqrt(np.dot(plane234,plane234))
        plane123 = plane123/plane123_len
        plane234 = plane234/plane234_len
        dihedral_cos = np.dot(plane123,plane234)
        plane_cross = np.cross(plane123,plane234)
        plane_cross = plane_cross 
        dihedral_sine = np.dot(plane_cross,bond32)
        torsion = np.arctan2(dihedral_cos, dihedral_sine)*180/np.pi
        torsion = -(torsion - 90)        
        torsions.append(torsion)
        
    int_coords = bonds + angles + torsions

    return int_coords

def invert_vectors(qm_dipders,mm_dipders):
    new_mm_dipders = []
    np_qm_dipders = np.array(qm_dipders)
    np_mm_dipders = np.array(mm_dipders)
    for qm, mm in zip(np_qm_dipders,np_mm_dipders):
        RMSD_p = np.sqrt(np.sum(np.subtract(qm, 1*mm)**2))
        RMSD_m = np.sqrt(np.sum(np.subtract(qm,-1*mm)**2))
        if RMSD_p <= RMSD_m:
            new_mm_dipders.append(1*mm)
        else:
            new_mm_dipders.append(-1*mm)
    
    return new_mm_dipders

def assign_by_displ(filenames):
    with open('ini.json', 'r') as file:
        ini = json.load(file)   
    with open(filenames['gau_json'], 'r') as file:
        all_gau_data = json.load(file)      
    with open(filenames['tink_json'], 'r') as file:
        all_tink_data = json.load(file)   
    number_of_atoms = all_gau_data['number_of_atoms']

    np_atom_masses = all_tink_data['masses']
    np_sqrt_masses = np.sqrt(np_atom_masses)
    weights = []
    # for a in np_atom_masses:
    #     weights.append([a])
    for a in np_sqrt_masses:
        weights.append([a])
        
    count = 0
        
    if ini['displacements'] == 'cartesians':
        QM_displ = np.array(all_gau_data['displs'])
        MM_displ = np.array(all_tink_data['displs'])
    elif ini['displacements'] == 'internals':
        QM_displ = []
        MM_displ = []
        qc = all_gau_data['coords']
        mc = all_tink_data['coords']
        qc_int = int_coords(filenames,qc)
        int_factors = []
        for qc_i in qc_int:
            if 0.5 < qc_i < 4.0:
                int_factors.append(1)
            else:
                int_factors.append(0.1)
        
        for qd,md in zip(all_gau_data['displs'],all_tink_data['displs']):
            
            QM_cart_dp = np.add(np.array(qc),0.05*np.array(qd))
            QM_cart_dm = np.subtract(np.array(qc),0.05*np.array(qd))
            MM_cart_dp = np.add(np.array(mc),0.05*np.array(md))
            MM_cart_dm = np.subtract(np.array(mc),0.05*np.array(md))
            QM_displ_int = np.subtract(int_coords(filenames,QM_cart_dp),int_coords(filenames,QM_cart_dm))/2/0.05
            MM_displ_int = np.subtract(int_coords(filenames,MM_cart_dp),int_coords(filenames,MM_cart_dm))/2/0.05
            
            QM_displ.append(np.multiply(QM_displ_int,int_factors))
            MM_displ.append(np.multiply(MM_displ_int,int_factors))
            
            # count = count +1
            # with open('LOG.LOG','a') as out_file:
            #     out_file.write('#### '+str(count)+'\n')
            #     for q_c,m_c,q,m in zip(int_coords(filenames,qc),int_coords(filenames,mc),QM_displ_int,MM_displ_int):
            #         out_file.write(str(q_c) + ' ' + str(m_c) + ' ' + str(q) + ' ' + str(m) + '\n')
                
    QM_freq  = np.array(all_gau_data['freqs'])
    MM_freq  = np.array(all_tink_data['freqs'])
    #QM_int  = np.array(all_gau_data['irints'])
    #MM_int  = np.array(all_tink_data['irints'])
    #QM_int  = np.divide(QM_int,np.mean(QM_int))
    #MM_int  = np.divide(QM_int,np.mean(MM_int))
    #QM_vec  = np.array(all_gau_data['norm_tdms'])
    #MM_vec  = np.array(all_tink_data['norm_tdms'])
    
    rmsd_matrix_v = []
    rmsd_matrix = []
    factor_matrix = []    
    for idx in range(0, len(MM_freq)):
        qm_x = QM_displ[idx]
        # if ini['displacements'] == 'cartesians':
        #     qm_x = np.multiply(qm_x, weights)
        qm_v = QM_freq[idx]
        #qm_I = QM_int[idx]
        #qm_vec = QM_vec[idx]

    
        rmsds = []
        rmsds_v = []
        factor = []

        for jdx in range(0, len(MM_freq)):
            mm_x_m =      MM_displ[jdx]
            mm_x_p = -1 * MM_displ[jdx]
            #mm_vec_p =      MM_vec[jdx]
            #mm_vec_m = -1 * MM_vec[jdx]
            
            # if ini['displacements'] == 'cartesians':
            #     mm_x_m = np.multiply(mm_x_m, weights)
            #     mm_x_p = np.multiply(mm_x_p, weights)
            
            mm_v   = MM_freq[jdx]
            #mm_I   = MM_int[jdx]
            
            # with 3*N coordinates that are evaluated and should correct on the 0.1 scalte
            # but just one frequency, which should be correct on the 1 cm-1 scale
            if jdx - len(MM_freq)/5 < idx < jdx - len(MM_freq)/5:
                if ini['displacements'] == 'cartesians':
                    rmsd_p = np.sqrt(np.sum((np.subtract(qm_x,mm_x_p))**2))            # this is = 2 is total mismatch
                    rmsd_m = np.sqrt(np.sum((np.subtract(qm_x,mm_x_m))**2))            # this is = 2 is total mismatch
                    # rmsd_v = 0.01*np.sqrt((np.subtract(qm_v,mm_v))**2)                 # this is = 2 if delta_v = 200
                    #rmsd_I = np.sqrt(np.sqrt((np.subtract(qm_I,mm_I))**2))
                    #rmsd_vec_p = np.sqrt(np.sum((np.subtract(qm_vec,mm_vec_p))**2)) 
                    #rmsd_vec_m = np.sqrt(np.sum((np.subtract(qm_vec,mm_vec_m))**2))
                
                elif ini['displacements'] == 'internals':    
                    rmsd_p = np.sqrt(np.average((np.subtract(qm_x,mm_x_p))**2))            # this is = 2 is total mismatch
                    rmsd_m = np.sqrt(np.average((np.subtract(qm_x,mm_x_m))**2))            # this is = 2 is total mismatch
                    # rmsd_v = 0.01*np.sqrt((np.subtract(qm_v,mm_v))**2)                 # this is = 2 if delta_v = 200
                    #rmsd_I = np.sqrt(np.sqrt((np.subtract(qm_I,mm_I))**2))
                    #rmsd_vec_p = np.sqrt(np.sum((np.subtract(qm_vec,mm_vec_p))**2)) 
                    #rmsd_vec_m = np.sqrt(np.sum((np.subtract(qm_vec,mm_vec_m))**2))
            else:
                rmsd_p = 999
                rmsd_m = 999
            
            rmsd_p = rmsd_p #+ rmsd_vec_p)/2 #+ rmsd_v #+ rmsd_I
            rmsd_m = rmsd_m #+ rmsd_vec_m)/2 #+ rmsd_v #+ rmsd_I
            
            if rmsd_p <= rmsd_m:
                rmsds.append(rmsd_p)
                factor.append(1)
            else:
                rmsds.append(rmsd_m)
                factor.append(-1)
            # rmsds_v.append(rmsd_v)
        rmsd_matrix_v.append(rmsds_v)
        rmsd_matrix.append(rmsds)
        factor_matrix.append(factor)   
        
    # print(rmsd_matrix[-1])
        
    if ini['displacements'] == 'cartesians':
        rmsd_threshold = 1
    elif ini['displacements'] == 'internals':
        rmsd_threshold = 1
        
    try:
        collect_lowrmsd_idxs = []
        for i in range(len(rmsd_matrix)):
            temp_lowrmsd_idxs = []
            for j in range(len(rmsd_matrix[0])):
                if rmsd_matrix[i][j] < rmsd_threshold and len(temp_lowrmsd_idxs) < 2:# and rmsd_matrix_v[i][j] < 2:
                    temp_lowrmsd_idxs.append(j)

            if temp_lowrmsd_idxs == []:
                temp_list = list(set(rmsd_matrix[i]))
                min_rmsds = [temp_list[0], temp_list[1]]#, temp_list[2]]
                for rmsd in min_rmsds:
                    for j in range(len(rmsd_matrix[0])):
                        if rmsd_matrix[i][j] == rmsd and rmsd < 1.5:
                            temp_lowrmsd_idxs.append(j)
            collect_lowrmsd_idxs.append(temp_lowrmsd_idxs)
        
        # remove again:
        for idx in range(len(collect_lowrmsd_idxs)):
            ele = []
            for jdx in collect_lowrmsd_idxs[idx]:
                ele.append(rmsd_matrix[idx][jdx])
            
            print(str(idx) + ' ' + ' '.join([str(x) for x in collect_lowrmsd_idxs[idx]]) + ' ' + ' '.join([str(x) for x in ele]))
        
        perms = list(itertools.product(*collect_lowrmsd_idxs))
        # print('####')
        # print(len(perms))
        # print(len(perms))
        new_perms = []
        for perm in perms:
            if len(perm) == len(set(perm)):
                new_perms.append(perm)
        total_rmsd = []
        # print(len(new_perms))
        for perm in new_perms:
            # print(perm)
            for idx in range(len(perm)):
                rmsd = 0
                jdx = perm[idx]
                rmsd = rmsd + rmsd_matrix[idx][jdx]
            total_rmsd.append(rmsd)
        perm_min = new_perms[np.argmin(total_rmsd)]
        np_assign = perm_min
        
    
    # for i in range(len(np_assign)):
    #     print(str(i) + ' ' + str(np_assign[i]))
    except:          
        min_i = 0
        min_j = 0
        global_min = 999999
        assignments = []
        for counter in range(0, len(MM_freq)):
            for i in range(len(rmsd_matrix)):
                for j in range(len(rmsd_matrix[0])):
                    # if rmsd_matrix_v[i][j] < 500:
                        if rmsd_matrix[i][j] < global_min:
                            global_min = rmsd_matrix[i][j]
                            min_i = i
                            min_j = j
            for i in range(len(rmsd_matrix)):
                rmsd_matrix[i][min_j] = 999999
            for j in range(len(rmsd_matrix)):
                rmsd_matrix[min_i][j] = 999999
            factor = factor_matrix[min_i][min_j]
            assignments.append([min_i,min_j,factor,global_min])
            global_min = 999999
        
        assignments = sorted(assignments)
        np_assign = (np.array(assignments).transpose())[1]
    
    a = [int(x) for x in np_assign]
    b = all_tink_data['freqs']
    c = all_tink_data['redmasses']
    d = all_tink_data['dipders']
    e = all_tink_data['displs']
    f = all_tink_data['norm_tdms']
    g = all_tink_data['irints']
    # d = all_tink_data['irints']
    # e = all_tink_data['irdipstr']
    # g = all_tink_data['tdms']

    
    a_new = []
    b_new = []
    c_new = []
    d_new = []
    e_new = []
    f_new = []
    g_new = []
    # h_new = []
    # i_new = []
    
    for idx in range(len(a)):
        for jdx in range(len(a)):
            if idx == a[jdx]:
                a_new.append(a[jdx])
                b_new.append(b[jdx])
                c_new.append(c[jdx])
                d_new.append(d[jdx])
                e_new.append(e[jdx])
                f_new.append(f[jdx])
                g_new.append(g[jdx])
                # g_new.append(g[jdx])
                # h_new.append(h[jdx])
                # i_new.append(i[jdx])
    
    all_tink_data['assignments'] = [int(x+1) for x in a]
    all_tink_data['freqs'] = b_new
    all_tink_data['redmasses'] = c_new
    all_tink_data['dipders'] = d_new
    all_tink_data['displs'] = e_new
    all_tink_data['norm_tdms']  = f_new
    all_tink_data['irints'] = g_new
    # all_tink_data['forceconsts'] = c_new
    # all_tink_data['irints'] = d_new
    # all_tink_data['irdipstr'] = e_new
    # all_tink_data['tdms'] = g_new
    
    with open(filenames['tink_json'], 'w') as file:
        json.dump(all_tink_data, file, indent = 4, ensure_ascii = False)