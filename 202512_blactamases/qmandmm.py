from auxfuncs import isint,flatlist,addzeros,rotmat,rotate,xyz_linear,xyz_Nby3,delete_annoying_files,same_atom,openjson
import shutil
import os
from subprocess import run
import numpy as np
import json

def is_gau_normal_termination(file):
    with open(file, 'r') as in_file:
        for line in in_file:
            if 'Normal termination' in line:
                return True
    return False

def wait_for_gau_normal_termination(file):
    wait = True
    while wait == True:
        if is_gau_normal_termination(file):
            wait = False
            
def is_tink_normal_termination(file):
    with open(file, 'r') as in_file:
        for line in in_file:
            if 'Normal Termination due to SmallGrad' in line:
                return True
            if 'Moments of Inertia and Principal Axes' in line:
                return True
            if 'Total Potential Energy' in line:
                return True
            if 'Vibrational Normal Mode' in line:
                return True

    return False

def wait_for_tink_normal_termination(file):
    wait = True
    while wait == True:
        if is_tink_normal_termination(file):
            wait = False            

# def read_gau_zmat(zmatfile,carts):
#     top_lines = []
#     with open(zmatfile, 'r') as in_file:
#         for line in in_file:        
#             try:
#                 if line.split()[0] in ['H','C','N','O','F','S','Cl']:
#                     top_lines.append(line)
#                 elif (line.split()[0]).split(',')[0] in ['H','C','N','O','F','S','Cl']:
#                     top_lines.append(line)
#             except:
#                 pass
        
#     num_lines = len(top_lines)
#     bottom_lines = []
#     with open(zmatfile, 'r') as in_file:
#         for line in in_file: 
#             # if 'Variables' in line:
#             #     bottom_lines.append(line)
#             if '=' in line:
#                 for idx in range(num_lines+1):
#                     if 'R'+str(idx)+'=' in line:
#                         bottom_lines.append(line)
#                     if 'A'+str(idx)+'=' in line:
#                         bottom_lines.append(line)
#                     if 'D'+str(idx)+'=' in line:
#                         bottom_lines.append(line)
                    
#     # for line in bottom_lines:
#     #     print(line)
    
#     lines = top_lines + [' Variables: \n'] + bottom_lines
                       
#     #check if there are similar dihedrals:
#     dihedral_connect = []
#     which_dihedral = []
#     new_top_lines = []
#     new_top_lines.append(top_lines[0])
#     new_top_lines.append(top_lines[1])
#     new_top_lines.append(top_lines[2])

#     for line in top_lines:
#         line1 = line.split()[0]
#         line2 = line1.split(',')
#         if len(line2) >= 7:
#             dihedral_connect.append([line2[1], line2[3], line2[5]])
#             which_dihedral.append(line2[6])
#         else:
#             dihedral_connect.append('none')
#             which_dihedral.append('none')            
            
#     new_dihedral_connect = []   
#     new_dihedrals = []
#     for idx in range(3,len(dihedral_connect)):
#         if dihedral_connect[idx] not in new_dihedral_connect:
#             new_dihedral_connect.append(dihedral_connect[idx])
#             new_top_lines.append(top_lines[idx])
            
#         else:
#             for jdx in range(1,idx):
#                 temp = [dihedral_connect[idx][0], dihedral_connect[idx][1], str(jdx)]
#                 if temp not in new_dihedral_connect:
#                     if str(jdx) != dihedral_connect[idx][0] and str(jdx) != dihedral_connect[idx][1]:
#                         xyz1 = carts[idx]
#                         xyz2 = carts[int(dihedral_connect[idx][0])-1]
#                         xyz3 = carts[int(dihedral_connect[idx][1])-1]
#                         xyz4 = carts[jdx-1]                      
#                         bond21 = np.subtract(xyz2,xyz1)
#                         bond32 = np.subtract(xyz3,xyz2)
#                         bond43 = np.subtract(xyz4,xyz3)
#                         bond21_len = np.sqrt(np.dot(bond21,bond21))
#                         bond32_len = np.sqrt(np.dot(bond32,bond32))
#                         bond43_len = np.sqrt(np.dot(bond43,bond43))
#                         bond21 = bond21/bond21_len
#                         bond32 = bond32/bond32_len
#                         bond43 = bond43/bond43_len
#                         plane123 = np.cross(bond21,bond32)
#                         plane234 = np.cross(bond32,bond43)
#                         plane123_len = np.sqrt(np.dot(plane123,plane123))
#                         plane234_len = np.sqrt(np.dot(plane234,plane234))
#                         plane123 = plane123/plane123_len
#                         plane234 = plane234/plane234_len
#                         dihedral_cos = np.dot(plane123,plane234)
                        
#                         plane_cross = np.cross(plane123,plane234)
#                         # plane_cross_len = np.sqrt(np.dot(plane_cross,plane_cross))
#                         plane_cross = plane_cross #/plane_cross_len
#                         dihedral_sine = np.dot(plane_cross,bond32)
#                         dihedral = np.arctan2(dihedral_cos, dihedral_sine)*180/np.pi
#                         dihedral = -(dihedral - 90)
#                         new_dihedral_connect.append(temp)
#                         break
#                     else:
#                         dihedral = 'None found'
#             if dihedral == 'None found':
#                 new_dihedral_connect.append(dihedral_connect[idx])
#                 new_top_lines.append(top_lines[idx])
#             else:
#                 line1 = top_lines[idx].split()[0]
#                 line2 = line1.split(',')
#                 new_line = ' '+\
#                            line2[0] + ',' +\
#                            line2[1] + ',' +\
#                            line2[2] + ',' +\
#                            line2[3] + ',' +\
#                            line2[4] + ',' +\
#                            str(jdx) + ',' +\
#                            line2[6] + ',' +\
#                            line2[7] + '\n'
#                 new_dihedrals.append([' ' + line2[6] + '=' + str(dihedral)])
#                 new_top_lines.append(new_line)
            
#     new_bottom_lines = []
    
#     # print(bottom_lines)
#     # print(new_dihedrals)
    
#     for idx in range(len(bottom_lines)):
#         found = False
#         for new_dih in new_dihedrals:
#             name1 = new_dih[0]
#             name2 = bottom_lines[idx]
#             name1 = name1.split('=')[0]
#             name2 = name2.split('=')[0]

#             if name2 == name1:
#                 # print(name1)
#                 # print(name2)

#                 new_bottom_lines.append(new_dih[0]+'\n')
#                 found = True
#                 break
#             else:
#                 found = False
#         if found == False:
#             new_bottom_lines.append(bottom_lines[idx])
                          
#     new_lines = new_top_lines + [' Variables: \n'] + new_bottom_lines
#     return new_lines

def read_gau_zmat(cartfile):
    # This needs fixing of equation for dihedrals, see above using tan2
    
    cartlines = read_gau_cart(cartfile)
    carts = read_gau_cart_coords(cartfile)
    
    bond_matrix = []
    angle_matrix = []
    dihedral_matrix = []
    
    bond_matrix = np.zeros((len(carts),len(carts)))
    for idx in range(len(carts)):
        for jdx in range(len(carts)):
            if idx == jdx:
                bond_matrix[idx,idx] = 999
            else:
                xyz1 = carts[idx]
                xyz2 = carts[jdx]
                bond = np.subtract(xyz1,xyz2)
                bond = np.sqrt(np.dot(bond,bond))
                bond_matrix[idx,jdx] = bond
                bond_matrix[jdx,idx] = bond
    angle_matrix = np.zeros((len(carts),len(carts),len(carts)))
    for idx in range(len(carts)):
        for jdx in range(len(carts)):
            for kdx in range(len(carts)):
                if idx == jdx or idx == kdx or jdx == kdx:
                    angle_matrix[idx,jdx,kdx] = 999
                else:
                    xyz1 = carts[idx]
                    xyz2 = carts[jdx]
                    xyz3 = carts[kdx]
                    bond12 = np.subtract(xyz1,xyz2)
                    bond32 = np.subtract(xyz3,xyz2)
                    bond12_len = np.sqrt(np.dot(bond12,bond12))
                    bond32_len = np.sqrt(np.dot(bond32,bond32))
                    bond12 = bond12/bond12_len
                    bond32 = bond32/bond32_len
                    angle = np.dot(bond12,bond32)
                    angle = np.arccos(angle)*180/np.pi
                    angle_matrix[idx,jdx,kdx] = angle
                    angle_matrix[kdx,jdx,idx] = angle
    dihedral_matrix = np.zeros((len(carts),len(carts),len(carts),len(carts)))
    for idx in range(len(carts)):
        for jdx in range(len(carts)):
            for kdx in range(len(carts)):
                for ldx in range(len(carts)):
                    if idx == jdx or idx == kdx or idx == ldx or jdx == kdx or jdx == ldx or kdx == ldx:
                        dihedral_matrix[idx,jdx,kdx,ldx] = 999
                    else:
                        xyz1 = carts[idx]
                        xyz2 = carts[jdx]
                        xyz3 = carts[kdx]
                        xyz4 = carts[ldx]
                        bond21 = np.subtract(xyz2,xyz1)
                        bond32 = np.subtract(xyz3,xyz2)
                        bond43 = np.subtract(xyz4,xyz3)
                        bond21_len = np.sqrt(np.dot(bond21,bond21))
                        bond32_len = np.sqrt(np.dot(bond32,bond32))
                        bond43_len = np.sqrt(np.dot(bond43,bond43))
                        bond21 = bond21/bond21_len
                        bond32 = bond32/bond32_len
                        bond43 = bond43/bond43_len
                        plane123 = np.cross(bond21,bond32)
                        plane234 = np.cross(bond32,bond43)
                        plane123_len = np.sqrt(np.dot(plane123,plane123))
                        plane234_len = np.sqrt(np.dot(plane234,plane234))
                        plane123 = plane123/plane123_len
                        plane234 = plane234/plane234_len
                        dihedral_cos = np.dot(plane123,plane234)
                        
                        plane_cross = np.cross(plane123,plane234)
                        # plane_cross_len = np.sqrt(np.dot(plane_cross,plane_cross))
                        plane_cross = plane_cross #/plane_cross_len
                        dihedral_sine = np.dot(plane_cross,bond32)
                        dihedral = np.arctan2(dihedral_cos, dihedral_sine)*180/np.pi
                        dihedral = -(dihedral - 90)
                        dihedral_matrix[idx,jdx,kdx,ldx] = dihedral
                        dihedral_matrix[ldx,kdx,jdx,idx] = dihedral




            
    adx = -1
    zmat_info = []
    for cart in cartlines:
        adx = adx + 1
        atom = cart.split()[0]
        zmat = [atom]
        
        #First atom:
        if adx == 0:
            zmat_info.append(zmat)
            continue
            
        #Second atom:
        if adx == 1:
            bond = bond_matrix[1][0]
            bond_matrix[0][1] = 999
            bond_matrix[1][0] = 999
            
            zmat = zmat + [0+1] + [bond]
            zmat_info.append(zmat)
            
        #Third atom:
        if adx == 2:
            bond = bond_matrix[2][0]
            bond_matrix[0][2] = 999
            bond_matrix[2][0] = 999
            
            angle = angle_matrix[2][0][1]
            angle_matrix[2][0][1] = 999
            angle_matrix[1][0][2] = 999
            
            zmat = zmat + [0+1] + [bond] + [1+1] + [angle]                
            zmat_info.append(zmat)
            
        #all other atoms:
        if adx > 2:
            bdx = np.argmin(bond_matrix[adx][0:adx])
            bond = bond_matrix[adx][bdx]
            bond_matrix[adx][bdx] = 999
            bond_matrix[bdx][adx] = 999
            
            for cdx in range(len(carts)):
                if cdx < adx:
                    if angle_matrix[adx][bdx][cdx] < 145:
                        angle = angle_matrix[adx][bdx][cdx]
                        break
                    elif angle_matrix[adx][bdx][cdx] < 175:
                        angle = angle_matrix[adx][bdx][cdx]
                        break
            angle_matrix[adx][bdx][cdx] = 999
            angle_matrix[cdx][bdx][adx] = 999
            
            for ddx in range(len(carts)):
                if ddx < adx:
                    if dihedral_matrix[adx][bdx][cdx][ddx] < 999:
                        dihedral = dihedral_matrix[adx][bdx][cdx][ddx]
                        break        
            dihedral_matrix[adx][bdx][cdx][ddx] = 999
            dihedral_matrix[ddx][cdx][bdx][adx] = 999
            for adx_temp in range(len(carts)):
                dihedral_matrix[adx_temp][bdx][cdx][ddx] = 999
                dihedral_matrix[ddx][cdx][bdx][adx_temp] = 999
            for ddx_temp in range(len(carts)):
                dihedral_matrix[ddx_temp][cdx][bdx][adx] = 999
                dihedral_matrix[adx][bdx][cdx][ddx_temp] = 999
            
            zmat = zmat + [bdx+1] + [bond] + [cdx+1] + [angle] + [ddx+1] + [dihedral]             
            zmat_info.append(zmat)
            
    zmat_top_lines = []
    zmat_bottom_lines = []
    idx = 0
    for zmat in zmat_info:
        idx = idx + 1 
        R_name = 'R'+str(idx)
        A_name = 'A'+str(idx)
        D_name = 'D'+str(idx)
        if idx == 1:
            zmat_top_lines.append(zmat[0])
        if idx == 2:
            string = [str(zmat[0]), str(zmat[1]), R_name]
            zmat_top_lines.append(','.join([str(x) for x in string]))
            zmat_bottom_lines.append(R_name+'='+str(zmat[2]))
        if idx == 3:
            string = [str(zmat[0]), str(zmat[1]), R_name, str(zmat[3]), A_name]
            zmat_top_lines.append(','.join([str(x) for x in string]))
            zmat_bottom_lines.append(R_name+'='+str(zmat[2]))
            zmat_bottom_lines.append(A_name+'='+str(zmat[4]))
        if idx > 3:
            string = [str(zmat[0]), str(zmat[1]), R_name, str(zmat[3]), A_name, str(zmat[5]), D_name]
            zmat_top_lines.append(','.join([str(x) for x in string]))
            zmat_bottom_lines.append(R_name+'='+str(zmat[2]))
            zmat_bottom_lines.append(A_name+'='+str(zmat[4]))
            zmat_bottom_lines.append(D_name+'='+str(zmat[6]))
    
    zmat_lines = zmat_top_lines + ['Variables:'] + zmat_bottom_lines
    
    new_zmat_lines = []
        
    for lines in zmat_lines:
        new_zmat_lines.append(lines + '\n')
            
    return new_zmat_lines

def read_gau_cart(file):
    lines = []
    with open(file, 'r') as in_file:
        for line in in_file:        
            try:
                if line.split()[0] in ['H','C','N','O','F','S','Cl']:
                    lines.append(line)
            except:
                pass
    # if lines == []:
    #     with open(file, 'r') as in_file:
    #         for line in in_file: 
    #             try:
    #                 if line.split()[0].split('(')[0] in ['H','C','N','O','F','S','Cl']:
    #                     line_new = line.split()[0].split('(')[0] + ' ' + ' '.join([x for x in line.split()[-3:]])
    #                     lines.append(line_new)
    #             except:
    #                 pass
    return lines

def read_gau_cart_coords(file):
    lines = []
    with open(file, 'r') as in_file:
        for line in in_file:        
            try:
                if line.split()[0] in ['H','C','N','O','F','S','Cl']:
                    lines.append([float(x) for x in line.split()[1:]])
            except:
                pass
    # if lines == []:
    #     with open(file, 'r') as in_file:
    #         for line in in_file: 
    #             try:
    #                 if line.split()[0].split('(')[0] in ['H','C','N','O','F','S','Cl']:
    #                     lines.append([float(x) for x in line.split()[-3:]])
    #             except:
    #                 pass
    return lines

def read_tink_xyz(file):
    tink_lines = []
    with open(file, 'r') as in_file:
        for line in in_file:        
            if len(line.split()) > 2:
                if isint(line.split()[0]) and line.split()[1] in ['H','C','N','O','F','S','Cl']:
                    tink_lines.append(line)
    return tink_lines

def read_tink_xyz_coords(file):
    tink_lines = []
    with open(file, 'r') as in_file:    
        for line in in_file:        
            if len(line.split()) > 2:
                if isint(line.split()[0]) and line.split()[1] in ['H','C','N','O','F','S','Cl']:
                    tink_lines.append([float(x) for x in line.split()[2:5]])
    return tink_lines

def generate_gau_init_files(filenames,ini):

    lines = read_tink_xyz(filenames['input_xyz'])
        
    with open(filenames['gau_init_com'], 'w') as out_file:
        out_file.write('%chk='+filenames['gau_init_chk']+'\n')
        out_file.write('%Mem=' + ini['mem'] + '\n')
        out_file.write('%nprocshared=' + ini['nproc'] + '\n')
        out_file.write('#p ' + ini['qm_level'] + ' nosymm opt=loose scf=(xqc,maxconventionalcycles=128) ')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('Init'+'\n')
        out_file.write('\n')
        out_file.write(ini['charge_spin']+'\n')
        for line in lines:
            out_file.write('   '.join(line.split()[1:5])+'\n')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('\n')
        
def generate_gau_cart_ef_files(filenames, suffs):

    for suff in suffs:
        org_file = filenames['gau_init_cart']
        ef_file  = filenames['gau_init_cart'].replace(filenames['gen_prefix'],filenames['gen_prefix']+suff)
        shutil.copy(org_file,ef_file)
        org_file = filenames['gau_init_zmat']
        ef_file  = filenames['gau_init_zmat'].replace(filenames['gen_prefix'],filenames['gen_prefix']+suff)
        shutil.copy(org_file,ef_file)

def generate_gau_opt_files(filenames,ini):

    carts = read_gau_cart_coords(filenames['gau_init_cart'])
    # lines = read_gau_zmat(filenames['gau_init_zmat'],carts)
    lines = read_gau_zmat(filenames['gau_init_cart'])
        
    with open(filenames['gau_opt_com'], 'w') as out_file:
        out_file.write('%chk='+filenames['gau_opt_chk']+'\n')
        out_file.write('%Mem=' + ini['mem'] + '\n')
        out_file.write('%nprocshared=' + ini['nproc'] + '\n')
        out_file.write('#p ' + ini['qm_level'] + ' nosymm opt=(z-matrix,maxcycles=256) int=ultrafine scf=(xqc,maxconventionalcycles=128) ')
        field_strength = int(50/5142.20674763/0.0001)
        if   '-50x' in filenames['gau_opt_com']:
            out_file.write('field=x+' + str(field_strength))
        elif '-50y' in filenames['gau_opt_com']:
            out_file.write('field=y+' + str(field_strength))
        elif '-50z' in filenames['gau_opt_com']:
            out_file.write('field=z+' + str(field_strength))
        elif '-100z' in filenames['gau_opt_com']:
            field_strength = 2*field_strength
            out_file.write('field=z+' + str(field_strength))
        elif '+50x' in filenames['gau_opt_com']:
            out_file.write('field=x-' + str(field_strength))
        elif '+50y' in filenames['gau_opt_com']:
            out_file.write('field=y-' + str(field_strength))
        elif '+50z' in filenames['gau_opt_com']:
            out_file.write('field=z-' + str(field_strength))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('Optimization'+'\n')
        out_file.write('\n')
        out_file.write(ini['charge_spin']+'\n')
        for line in lines:
            out_file.write(line)
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('\n')
        
        
def generate_gau_harm_files(filenames,ini):
    
    lines = read_gau_cart(filenames['gau_opt_cart'])
    
    with open(filenames['gau_harm_com'], 'w') as out_file:
        out_file.write('%chk='+filenames['gau_harm_chk']+'\n')
        out_file.write('%Mem=' + ini['mem'] + '\n')
        out_file.write('%nprocshared=' + ini['nproc'] + '\n')
        out_file.write('#p ' + ini['qm_level'] + ' nosymm freq=printderivatives int=ultrafine scf=(xqc,maxconventionalcycles=128) ')
        field_strength = int(50/5142.20674763/0.0001)
        if   '-50x' in filenames['gau_opt_com']:
            out_file.write('field=x+' + str(field_strength))
        elif '-50y' in filenames['gau_opt_com']:
            out_file.write('field=y+' + str(field_strength))
        elif '-50z' in filenames['gau_opt_com']:
            out_file.write('field=z+' + str(field_strength))
        elif '-100z' in filenames['gau_opt_com']:
            field_strength = 2*field_strength
            out_file.write('field=z+' + str(field_strength))
        elif '+50x' in filenames['gau_opt_com']:
            out_file.write('field=x-' + str(field_strength))
        elif '+50y' in filenames['gau_opt_com']:
            out_file.write('field=y-' + str(field_strength))
        elif '+50z' in filenames['gau_opt_com']:
            out_file.write('field=z-' + str(field_strength))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('Optimization'+'\n')
        out_file.write('\n')
        out_file.write(ini['charge_spin']+'\n')
        for line in lines:
            if len(line.split()) > 1:
                out_file.write('   '.join(line.split())+'\n')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('\n')

# def generate_gau_init_files(filenames,ini):

#     lines = read_tink_xyz(filenames['input_xyz'])
        
#     with open(filenames['gau_init_com'], 'w') as out_file:
#         out_file.write('%chk='+filenames['gau_init_chk']+'\n')
#         out_file.write('%Mem=1000MB'+'\n')
#         out_file.write('%nprocshared=4'+'\n')
#         out_file.write('#p mp2/6-31g** nosymm opt scf=(xqc,maxconventionalcycles=128) ')
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('Init'+'\n')
#         out_file.write('\n')
#         out_file.write(ini['charge_spin']+'\n')
#         for line in lines:
#             out_file.write('   '.join(line.split()[1:5])+'\n')
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('\n')

# def generate_gau_opt_files(filenames,ini):
    
#     lines = read_gau_cart(filenames['gau_opt_cart'])
        
#     with open(filenames['gau_opt_com'], 'w') as out_file:
#         out_file.write('%chk='+filenames['gau_opt_chk']+'\n')
#         out_file.write('%Mem=1000MB'+'\n')
#         out_file.write('%nprocshared=4'+'\n')
#         out_file.write('#p mp2/6-31g** nosymm opt int=ultrafine scf=(xqc,maxconventionalcycles=128) ')
#         field_strength = int(50/5142.20674763/0.0001)
#         if '-x' in filenames['gau_opt_com']:
#             out_file.write('field = x+' + str(field_strength))
#         elif '-y' in filenames['gau_opt_com']:
#             out_file.write('field = y+' + str(field_strength))
#         elif '-z' in filenames['gau_opt_com']:
#             out_file.write('field = z+' + str(field_strength))
#         elif '+x' in filenames['gau_opt_com']:
#             out_file.write('field = x-' + str(field_strength))
#         elif '+y' in filenames['gau_opt_com']:
#             out_file.write('field = y-' + str(field_strength))
#         elif '+z' in filenames['gau_opt_com']:
#             out_file.write('field = z-' + str(field_strength))
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('Optimization'+'\n')
#         out_file.write('\n')
#         out_file.write(ini['charge_spin']+'\n')
#         for line in lines:
#             out_file.write('   '.join(line.split()[1:5])+'\n')
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('\n')
        
# def generate_gau_cart_ef_files(filenames, suffs):
    
#     for suff in suffs:
#         org_file = filenames['gau_opt_cart']
#         ef_file  = filenames['gau_opt_cart'].replace(filenames['gen_prefix'],filenames['gen_prefix']+suff)
#         shutil.copy(org_file,ef_file)

        
# def generate_gau_harm_files(filenames,ini):
    
#     lines = read_gau_cart(filenames['gau_opt_cart'])
    
#     with open(filenames['gau_harm_com'], 'w') as out_file:
#         out_file.write('%chk='+filenames['gau_harm_chk']+'\n')
#         out_file.write('%Mem=1000MB'+'\n')
#         out_file.write('%nprocshared=4'+'\n')
#         out_file.write('#p mp2/6-31g** nosymm freq=printderivatives int=ultrafine scf=(xqc,maxconventionalcycles=128) ')
#         field_strength = int(50/5142.20674763/0.0001)
#         if '-x' in filenames['gau_opt_com']:
#             out_file.write('field=x+' + str(field_strength))
#         if '-y' in filenames['gau_opt_com']:
#             out_file.write('field=y+' + str(field_strength))
#         if '-z' in filenames['gau_opt_com']:
#             out_file.write('field=z+' + str(field_strength))
#         if '+x' in filenames['gau_opt_com']:
#             out_file.write('field=x-' + str(field_strength))
#         if '+y' in filenames['gau_opt_com']:
#             out_file.write('field=y-' + str(field_strength))
#         if '+z' in filenames['gau_opt_com']:
#             out_file.write('field=z-' + str(field_strength))
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('Optimization'+'\n')
#         out_file.write('\n')
#         out_file.write(ini['charge_spin']+'\n')
#         for line in lines:
#             out_file.write('   '.join(line.split())+'\n')
#         out_file.write('\n')
#         out_file.write('\n')
#         out_file.write('\n')



def submit_gau_init(filenames):
    if os.path.isfile(filenames['gau_init_out']) == False:
        callstring = 'g16 < %s > %s' %(filenames['gau_init_com'],filenames['gau_init_out'])
        run(callstring, shell=True)
        wait_for_gau_normal_termination(filenames['gau_init_out'])
    if os.path.isfile(filenames['gau_init_cart']) == False:
        callstring = 'newzmat -ichk -ocart -step 999 %s  %s' %(filenames['gau_init_chk'],filenames['gau_init_cart'])
        run(callstring, shell=True)
    if os.path.isfile(filenames['gau_init_zmat']) == False:
        callstring = 'newzmat -ichk -ozmat -step 999 -rebuildzmat %s  %s' %(filenames['gau_init_chk'],filenames['gau_init_zmat'])
        run(callstring, shell=True)

def submit_gau_opt(filenames):
    if os.path.isfile(filenames['gau_opt_out']) == False:
        callstring = 'g16 < %s > %s' %(filenames['gau_opt_com'],filenames['gau_opt_out'])
        run(callstring, shell=True)
        wait_for_gau_normal_termination(filenames['gau_opt_out'])
    if os.path.isfile(filenames['gau_opt_cart']) == False:
        callstring = 'newzmat -ichk -ocart -step 999 %s  %s' %(filenames['gau_opt_chk'],filenames['gau_opt_cart'])
        run(callstring, shell=True)
    if os.path.isfile(filenames['gau_opt_zmat']) == False:
        callstring = 'newzmat -ichk -ozmat -step 999 -rebuildzmat -noround %s  %s' %(filenames['gau_opt_chk'],filenames['gau_opt_zmat'])
        run(callstring, shell=True)
    if os.path.isfile(filenames['gau_opt_pdb']) == False:
        callstring = 'newzmat -ichk -opdb -step 999 -rebuildzmat -noround %s  %s' %(filenames['gau_opt_chk'],filenames['gau_opt_pdb'])
        run(callstring, shell=True)
        
def submit_gau_harm(filenames):
    if os.path.isfile(filenames['gau_harm_out']) == False:
        callstring = 'g16 < %s > %s' %(filenames['gau_harm_com'],filenames['gau_harm_out'])
        run(callstring, shell=True)
        wait_for_gau_normal_termination(filenames['gau_harm_out'])
    if os.path.isfile(filenames['gau_harm_fchk']) == False:
        callstring = 'formchk %s %s' %(filenames['gau_harm_chk'],filenames['gau_harm_fchk'])  
        run(callstring, shell=True)
    if os.path.isfile(filenames['gau_harm_cart']) == False:
        callstring = 'newzmat -ichk -ocart -step 999 %s  %s' %(filenames['gau_harm_chk'],filenames['gau_harm_cart'])
        run(callstring, shell=True)
    if os.path.isfile(filenames['gau_harm_zmat']) == False:
        callstring = 'newzmat -ichk -ozmat -step 999 -rebuildzmat -noround %s  %s' %(filenames['gau_harm_chk'],filenames['gau_harm_zmat'])
        run(callstring, shell=True)
        
def qm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor):
    atoms = []
    coords = []
    masses = []
    freqs = []
    irints = []
    irdipstr = []
    redmass = []
    forceconst = []
    dipders = []
    displs = []
    norm_tdms = []
    
    for line in read_gau_cart(filenames['gau_opt_cart']):        
        atoms.append(line.split()[0])
        coords.append([float(x) for x in line.split()[1:]])

            
    with open(filenames['gau_harm_out'], 'r') as in_file:
        for line in in_file:
            if 'Frequencies --' in line:
                freqs.append([float(x) for x in line.split()[-3:]])
            if 'Red. masses --' in line:
                redmass.append([float(x) for x in line.split()[-3:]])
            if 'Frc consts  --' in line:
                forceconst.append([float(x)*143.836 for x in line.split()[-3:]]) # in kcal/molA^2
            if 'IR Inten    --' in line:
                irints.append([float(x) for x in line.split()[-3:]])
            if 'Dipole derivative wrt mode' in line:
                dipders.append([float(x.replace('D','E')) for x in line.split()[-3:]])
                
    freqs = flatlist(freqs)
    redmass = flatlist(redmass)
    forceconst = flatlist(forceconst)
    irints = flatlist(irints)
    irdipstr = []
    for i,v in zip(irints,freqs):
        irdipstr.append(i*0.399212958/v) # in D^2
                
    dipders = flatlist(dipders)
    n = int(len(dipders)/3)
    dipders = np.reshape(dipders, (n,3)).tolist()
    tdms = []
    for dipder, v in zip(dipders, freqs):
        tdms.append(np.multiply(dipder,np.sqrt(0.399212958 / v)).tolist()) # in Debye
        
    for tdm in tdms:
        abs_tdm = np.sqrt(np.dot(tdm,tdm))
        norm_tdms.append((tdm/abs_tdm).tolist()) 
                
    with open(filenames['gau_harm_fchk'], 'r') as in_file:
        write = False
        for line in in_file:
            if 'Vib-E2' in line:
                break
            if write == True:
                masses.append([float(x) for x in line.split()])
            if 'Vib-AtMass' in line:
                write = True
                
    masses = flatlist(masses)

    with open(filenames['gau_harm_fchk'], 'r') as in_file:
        write = False
        for line in in_file:
            if 'ClPar MaxAn' in line:
                break
            if 'Anharmonic Number of Normal Modes' in line:
                break
            if write == True:
                displs.append([float(x) for x in line.split()])
            if 'Vib-Modes' in line:
                write = True
                
    displs = flatlist(displs)
    n = int(len(displs)/3/number_of_atoms)
    m = number_of_atoms
    displs = np.reshape(displs, (n,m,3)).tolist()
    
    all_gau_data = {'number_of_atoms':number_of_atoms,
                    'atoms':atoms,
                    'coords':coords,
                    'masses':masses,
                    'number_of_NMs':number_of_NMs,
                    'freqs':[v*scaling_factor for v in freqs],
                    'irints':irints,
                    'redmasses':redmass,
                    'dipders':dipders,
                    'displs':displs,
                    'norm_tdms':norm_tdms}#,
                    # 'irdipstr':irdipstr,
                    # 'tdms':tdms,
                    # 'forceconsts':forceconst}
    
    with open(filenames['gau_json'], 'w') as file:
        json.dump(all_gau_data, file, indent = 4, ensure_ascii = False)
    
    
def generate_tink_min_files(filenames, number_of_atoms):
    tink_lines = read_tink_xyz(filenames['input_xyz'])          
    gau_lines = read_gau_cart(filenames['gau_opt_cart'])
    
    with open(filenames['tink_min_xyz'], 'w') as out_file:
        out_file.write(str(number_of_atoms).rjust(6)+'\n')
        for gl,tl in zip(gau_lines,tink_lines):
            str1 = tl.split()[0].rjust(6)
            str2 = tl.split()[1].rjust(4)
            str3 = ' '.join([x.rjust(18) for x in gl.split()[1:]])
            str4 = ' '.join([x.rjust(6) for x in tl.split()[5:]])
            out_file.write(str1 + str2 + str3 + str4 + '\n')
            
def submit_tink_min(filenames):
    try:
        run('rm '+filenames['tink_min_xyz']+'_2',shell=True)
    except:
        pass
    
    callstring = 'minimize %s -k %s 0.1 > %s' %(filenames['tink_min_xyz'],filenames['tink_temp_key'],filenames['tink_min_out'])
    run(callstring, shell=True)
    wait_for_tink_normal_termination(filenames['tink_min_out'])
    
def align_tink_xyz(filenames,number_of_atoms):
    tink_lines = read_tink_xyz(filenames['input_xyz'])          
    tink_coords = read_tink_xyz_coords(filenames['tink_min_xyz']+'_2')          
    gau_coords = read_gau_cart_coords(filenames['gau_opt_cart'])  
    all_gau_data = openjson(filenames['gau_json']) 
    
    # print(tink_coords)
    # print(gau_coords)
    
    np_QM_xyz = np.array(gau_coords)
    np_MM_xyz = np.array(tink_coords)
    np_atom_masses = np.array([[x] for x in all_gau_data['masses']])
    
    # QM_center = np.average(np_QM_xyz,0)
    # MM_center = np.average(np_MM_xyz,0)
    
    # np_QM_xyz = np_QM_xyz - QM_center
    # np_MM_xyz = np_MM_xyz - MM_center
      
    QM_mass_center = np.sum(np.multiply(np_QM_xyz,np_atom_masses),axis=0)/np.sum(np_atom_masses)
    MM_mass_center = np.sum(np.multiply(np_MM_xyz,np_atom_masses),axis=0)/np.sum(np_atom_masses)
        
    np_QM_xyz = np_QM_xyz - QM_mass_center
    np_MM_xyz = np_MM_xyz - MM_mass_center
    
    scaled_QM_xyz = np.multiply(np_QM_xyz, np_atom_masses)
    scaled_MM_xyz = np.multiply(np_MM_xyz, np_atom_masses)
    
    Rmat = rotmat(scaled_MM_xyz,scaled_QM_xyz)
    MM_xyz_new = rotate(np_MM_xyz,Rmat) + QM_mass_center    
    
    new_tink_coords = np.round(MM_xyz_new,8).tolist()
    
    with open(filenames['tink_min_aligned_xyz'], 'w') as out_file:
        out_file.write(str(number_of_atoms).rjust(6)+'\n')
        for tc,tl in zip(new_tink_coords,tink_lines):
            str1 = tl.split()[0].rjust(6)
            str2 = tl.split()[1].rjust(4)
            str3 = ' '.join([str(x).rjust(14) for x in tc])
            str4 = ' '.join([x.rjust(6) for x in tl.split()[5:]])
            out_file.write(str1 + str2 + str3 + str4 + '\n')
            
def mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,which):
    atoms = []
    coords = []
    masses = []
    freqs = []
    irints = []
    irdipstr = []
    redmass = []
    forceconst = []
    dipders = []
    displs = []
    tdms = []
    norm_tdms = []
    
    with open('ini.json', 'r') as file:
        ini = json.load(file)
    with open(filenames['gau_json'], 'r') as file:
        all_gau_data = json.load(file)      
        
    atoms = all_gau_data['atoms']
    masses = all_gau_data['masses']
    coords = read_tink_xyz_coords(filenames['tink_min_aligned_xyz'])
    
    # if which == 'qm_based':
    
    #     for nm in range(0,number_of_NMs):
    #         suffix = '_'+addzeros(nm+1) 
                           
    #         nrgs = []
    #         if os.path.isfile(filenames['tink_nrg_out']+suffix):
    #             with open(filenames['tink_nrg_out']+suffix) as in_file:
    #                 for line in in_file:
    #                     if 'Total Potential Energy' in line:
    #                         nrgs.append(float(line.split()[-2]))
    #         dipoles = []
    #         if os.path.isfile(filenames['tink_dip_out']+suffix):
    #             with open(filenames['tink_dip_out']+suffix) as in_file:
    #                 for line in in_file:
    #                     if 'Dipole X,Y,Z-Components' in line:
    #                         dipoles.append([float(x) for x in line.split()[-3:]])   
                            
    #         if nrgs != []:
    #             if ini['displacements'] == 'cartesians':
    #                 stepsizes = [-0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20]
    #             elif ini['displacements'] == 'internals':
    #                 stepsizes = [-0.10, -0.05, 0.00, 0.05, 0.10]
    #             # print(nrgs)
    #             nrgs_fit = np.polyfit(stepsizes,nrgs,3)
    #             forceconst.append(2*nrgs_fit[-3])
    #             # print(2*nrgs_fit[-3])
    #             freq = np.sqrt(2*nrgs_fit[-3]*4.184*1e16/all_gau_data['redmasses'][nm])
    #             freq = 1000*freq/(2*np.pi*299792458)
    #             freqs.append(freq)
    
    #         if dipoles != []:        
    #             dip_fit = np.polyfit(stepsizes,dipoles,3)
    #             dipder = dip_fit[-2]                                     # in D/Ang
    #             dipder = dipder * 0.52917721092                          # in D/Bohr
    #             dipder = dipder * 0.39343                                # in e
    #             dipder = dipder / np.sqrt(all_gau_data['redmasses'][nm]) # in e/sqrt(amu)
    #             dipder = dipder * 31.223068                              # in sqrt(km/mol)
    
    #             dipders.append(dipder.tolist())
    #             irints.append(np.dot(dipder,dipder))
    #             tdm = np.multiply(dipder,np.sqrt(0.399212958 / all_gau_data['freqs'][nm])).tolist() # in Debye
    #             tdms.append(tdm)
    #             irdipstr.append(np.dot(tdm,tdm))
    #             abs_tdm = np.sqrt(np.dot(tdm,tdm))
    #             norm_tdms.append((tdm/abs_tdm).tolist())
                
    # elif which == 'tinker_vibrate':

    #     if os.path.isfile(filenames['tink_vib_out']):
        
    #         with open(filenames['tink_vib_out'], 'r') as in_file:
    #             write = False
    #             for line in in_file:
    #                 if 'Vibrational Normal Mode' in line:
    #                     write = True
    #                 if write == True:
    #                     if 'Vibrational Normal Mode' in line:
    #                         freq_temp = float(line.split()[-2])
    #                         freqs.append(freq_temp)
    #                     if len(line.split()) == 4:
    #                         first = line.split()[0]
    #                         if isint(first):
    #                             displs.append([float(x) for x in line.split()[-3:]])
    #         freqs = freqs[-number_of_NMs:]
    #         displs = flatlist(displs)
    #         n = int(len(displs)/3/number_of_atoms)
    #         m = number_of_atoms
    #         displs = np.reshape(displs, (n,m,3)).tolist()
    #         displs = displs[-number_of_NMs:]
                        
       
    #     for nm in range(6,number_of_NMs+6):
    #         # +6 or -6 are here to account for vibrate also giving translational modes
    #         stepsizes = [] 
    #         suffix = '_'+addzeros(nm+1)  
    #         replxyz = addzeros(nm+1)               
    #         nma_displ = displs[nm-6]
    #         if os.path.isfile(filenames['tink_vib_xyz'].replace('xyz',replxyz)):
    #             with open(filenames['tink_vib_xyz'].replace('xyz',replxyz),'r') as in_file:
    #                 coords_for_stepsize = []
    #                 temp_coords = []
    #                 for line in in_file:
    #                     if len(line.split()) > 1:
    #                         temp_coords.append([float(x) for x in line.split()[2:5]])
    #                     elif len(line.split()) == 1:
    #                         coords_for_stepsize.append(temp_coords)
    #                         temp_coords = []
    #                 coords_for_stepsize.append(temp_coords)
    #             coords_for_stepsize.pop(0)
    #             for c in coords_for_stepsize:
    #                 tink_displ = np.subtract(c,coords)
    #                 tink_displ = xyz_linear(tink_displ)
    #                 nma_displ = xyz_linear(nma_displ)
    #                 stepsize = []
    #                 for d1,d2 in zip(tink_displ,nma_displ):
    #                     try:
    #                         d1_div_d2 = d1/d2   
    #                         if d1_div_d2 == np.inf or d1_div_d2 == -np.inf or np.isnan(d1_div_d2):
    #                             continue
    #                         else:
    #                             stepsize.append(d1/d2)
    #                     except:
    #                         continue
    #                 stepsizes.append(np.mean(stepsize))
                               
    #         nrgs = []
    #         if os.path.isfile(filenames['tink_nrg_out']+suffix):
    #             with open(filenames['tink_nrg_out']+suffix) as in_file:
    #                 for line in in_file:
    #                     if 'Total Potential Energy' in line:
    #                         nrgs.append(float(line.split()[-2]))
    #         dipoles = []
    #         if os.path.isfile(filenames['tink_dip_out']+suffix):
    #             with open(filenames['tink_dip_out']+suffix) as in_file:
    #                 for line in in_file:
    #                     if 'Dipole X,Y,Z-Components' in line:
    #                         dipoles.append([float(x) for x in line.split()[-3:]])   
                            
    #         if nrgs != []:
    #             nrgs_fit = np.polyfit(stepsizes,nrgs,4)
    #             forceconst.append(2*nrgs_fit[-3])
    #             rm = (1000/(2*np.pi*299792458))**2 * (2*nrgs_fit[-3]*4.184*1e16) / freqs[nm-6]**2
    #             redmass.append(rm)
    
    #         if dipoles != []:        
    #             dip_fit = np.polyfit(stepsizes,dipoles,3)
    #             dipder = dip_fit[-2]                                     # in D/Ang
    #             dipder = dipder * 0.52917721092                          # in D/Bohr
    #             dipder = dipder * 0.39343                                # in e
    #             dipder = dipder / np.sqrt(rm) # in e/sqrt(amu)
    #             dipder = dipder * 31.223068                              # in sqrt(km/mol)
    
    #             dipders.append(dipder.tolist())
    #             irints.append(np.dot(dipder,dipder))
    #             tdm = np.multiply(dipder,np.sqrt(0.399212958 / rm)).tolist() # in Debye
    #             tdms.append(tdm)
    #             irdipstr.append(np.dot(tdm,tdm))
    #             abs_tdm = np.sqrt(np.dot(tdm,tdm))
    #             norm_tdms.append((tdm/abs_tdm).tolist())                
                
    if which == 'manual_nma':

        if os.path.isfile(filenames['tink_vib_out']):
        
            with open(filenames['tink_vib_out'], 'r') as in_file:
                write = False
                for line in in_file:
                    if 'Vibrational Normal Mode' in line:
                        write = True
                    if write == True:
                        if 'Vibrational Normal Mode' in line:
                            freq_temp = float(line.split()[-2])
                            freqs.append(freq_temp)
                        if len(line.split()) == 4:
                            first = line.split()[0]
                            if isint(first):
                                displs.append([float(x) for x in line.split()[-3:]])
            freqs = freqs[6:]
            displs = flatlist(displs)
            n = int(len(displs)/3/number_of_atoms)
            m = number_of_atoms
            displs = np.reshape(displs, (n,m,3)).tolist()
            displs = displs[6:]
             
            with open(filenames['tink_vib_out'], 'r') as in_file:
                write = True
                for line in in_file:
                    if 'Vibrational Normal Mode' in line:
                        write = False
                        break
                    if write == True:
                        if len(line.split()) > 1:
                            if isint(line.split()[0]):
                                rm = float(line.split()[-1])
                                redmass.append(rm)
                                dipder = [float(x) for x in line.split()[2:5]]
                                dipders.append(dipder)
                                
            dipders = dipders[6:]
            redmass = redmass[6:]
                        
        for nm in range(len(freqs)):                                
            if dipders != []:  
                for dipder,rm in zip(dipders,redmass):
                    irints.append(np.dot(dipder,dipder))
                    tdm = np.multiply(dipder,np.sqrt(0.399212958 / rm)).tolist() # in Debye
                    tdms.append(tdm)
                    irdipstr.append(np.dot(tdm,tdm))
                    abs_tdm = np.sqrt(np.dot(tdm,tdm))
                    norm_tdms.append((tdm/abs_tdm).tolist())                
                           
    delete_annoying_files(filenames['tink_dip_out'],number_of_NMs+7)
    delete_annoying_files(filenames['tink_nrg_out'],number_of_NMs+7)
    delete_annoying_files(filenames['tink_vib_xyz'],number_of_NMs+7)
               
    # all_tink_data = {'number_of_atoms':number_of_atoms,
                     # 'atoms':atoms,
                     # 'coords':coords,
                     # 'masses':masses,
                     # 'number_of_NMs':number_of_NMs,
                     # 'freqs':[v*scaling_factor for v in freqs],
                     # 'forceconsts':forceconst,
                     # 'irints':irints,
                     # 'irdipstr':irdipstr,
                     # 'redmasses':redmass,
                     # 'dipders':dipders,
                     # 'tdms':tdms,
                     # 'displs':displs,
                     # 'norm_tdms':norm_tdms}
    
    all_tink_data = {'number_of_atoms':number_of_atoms,
                     'atoms':atoms,
                     'coords':coords,
                     'masses':masses,
                     'number_of_NMs':number_of_NMs,
                     'freqs':[v for v in freqs],
                     'irints':irints,
                     'redmasses':redmass,
                     'dipders':dipders,
                     'displs':displs,
                     'norm_tdms':norm_tdms}#,
                     # 'tdms':tdms,
                     # 'irdipstr':irdipstr,
                     # 'forceconsts':forceconst}
    
    # if os.path.isfile(filenames['tink_json']):
    #     last_file = False
    #     i = 0
    #     while last_file == False:
    #         i = i + 1
    #         if os.path.isfile(filenames['tink_json']+'_'+str(i)):
    #             continue
    #         else:
    #             shutil.copy(filenames['tink_json'],filenames['tink_json']+'_'+str(i))
            
    
    with open(filenames['tink_json'], 'w') as file:
        json.dump(all_tink_data, file, indent = 4, ensure_ascii = False)
        
def readkey(keyfile):

    MM_p = {'vdw': [],
            'bond': [],
            'angle': [],
            'strbnd': [],
            'improper': [],
            'opbend': [],
            'torsion': []
            }

    with open(keyfile, 'r') as in_file:
        write = False
        for line in in_file:
            if 'PARAMETERS TO OPTIMIZE' in line:
                write = True
            if 'DO NOT CHANGE' in line:
                write = False
            if write == True:
                if len(line.split()) > 1:
                    if line.split()[0] in ['vdw', 'bond', 'angle', 'strbnd', 'improper']:
                        info = ' '.join(x for x in line.split()[0:-2])
                        param = [float(x) for x in line.split()[-2:]]

                        if 'vdw' in line.split()[0]:
                            MM_p['vdw'].append([info, param])
                        elif 'bond' in line.split()[0]:
                            MM_p['bond'].append([info, param])
                        elif 'angle' in line.split()[0]:
                            MM_p['angle'].append([info, param])
                        elif 'strbnd' in line.split()[0]:
                            MM_p['strbnd'].append([info, param])
                        elif 'improper' in line.split()[0]:
                            MM_p['improper'].append([info, param])
                    if line.split()[0] in ['torsion']:
                        info = ' '.join(x for x in line.split()[0:5])
                        p1 = float(line.split()[5])
                        p2 = float(line.split()[8])
                        p3 = float(line.split()[11])
                        a1 = float(line.split()[6])
                        a2 = float(line.split()[9])
                        a3 = float(line.split()[12])
                        # if len(line.split()) > 14:
                        #     p4 = float(line.split()[14])
                        #     a4 = float(line.split()[15])
                        # else:
                        #     p4 = 0
                        #     a4 = 0
                        # if len(line.split()) > 17:
                        #     p5 = float(line.split()[17])
                        #     a5 = float(line.split()[18])
                        # else:
                        #     p5 = 0
                        #     a5 = 0
                        # if len(line.split()) > 20:
                        #     p6 = float(line.split()[20])
                        #     a6 = float(line.split()[21])
                        # else:
                        #     p6 = 0
                        #     a6 = 0
                        # param = [p1, p2, p3, p4, p5, p6]
                        # angles = [a1, a2, a3, a4, a5, a6]
                        param = [p1, p2, p3]
                        angles = [a1, a2, a3]
                        MM_p['torsion'].append([info, param, angles])

    return MM_p  

def init_b0(filenames,cart_coords,number_of_atoms):
        
    t2as = []
    all_bond_idxs = []
    all_angle_idxs = []

    with open(filenames['input_xyz'], 'r') as xyz_file:
        for line in xyz_file:        
            try:
                if int(line) == number_of_atoms:
                    continue
            except:
                pass
            atomtype = int(line.split()[5])
            atomidx = int(int(line.split()[0])-1)
            t2as.append([atomtype, atomidx])

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

    xyz_bond_ats = []
    xyz_angle_ats = []
    for bond_idxs in set_bond_idxs:
        temp_bond_ats = []
        idx1 = bond_idxs[0]
        idx2 = bond_idxs[1]
        for t2a in t2as:
            if idx1 == t2a[1]:
                temp_bond_ats.append(t2a[0])
        for t2a in t2as:
            if idx2 == t2a[1]:
                temp_bond_ats.append(t2a[0])
        xyz_bond_ats.append(temp_bond_ats)
    for angle_idxs in set_angle_idxs:
        temp_angle_ats = []
        idx1 = angle_idxs[0]
        idx2 = angle_idxs[1]
        idx3 = angle_idxs[2]
        for t2a in t2as:
            if idx1 == t2a[1]:
                temp_angle_ats.append(t2a[0])
        for t2a in t2as:
            if idx2 == t2a[1]:
                temp_angle_ats.append(t2a[0])
        for t2a in t2as:
            if idx3 == t2a[1]:
                temp_angle_ats.append(t2a[0])
        xyz_angle_ats.append(temp_angle_ats)

    key_bond_ats = []
    key_angle_ats = []
    with open(filenames['input_key'], 'r') as key_file:  
        for line in key_file:
            if 'bond' in line:
                key_bond_ats.append([int(line.split()[1]), int(line.split()[2])])
    with open(filenames['input_key'], 'r') as key_file:  
        for line in key_file:
            if 'angle' in line:
                key_angle_ats.append([int(line.split()[1]), int(line.split()[2]), int(line.split()[3])])
                
    key_bonds = []  
    key_angles = []
    for key_at in key_bond_ats:
        temp_bonds = []
        for xyz_at, bond in zip(xyz_bond_ats,bonds):
            if key_at == xyz_at or key_at == xyz_at[::-1]:
                temp_bonds.append(bond)
        key_bonds.append(np.mean(temp_bonds))
        
    for key_at in key_angle_ats:
        temp_angles = []
        for xyz_at, angle in zip(xyz_angle_ats,angles):
            if key_at == xyz_at or key_at == xyz_at[::-1]:
                temp_angles.append(angle)
        key_angles.append(np.mean(temp_angles))                
                            
    return key_bonds + key_angles

def b0_params(MM_params):

    MM_p = []
    b_top = [] 
    b_low = []
    if MM_params['bond'] != []:
        for element in MM_params['bond']:
            MM_p.append(element[1][1])
            b_top.append(2.2)
            b_low.append(0.8)
    if MM_params['angle'] != []:
        for element in MM_params['angle']:
            MM_p.append(element[1][1])
            b_top.append(180)
            b_low.append(0)
    if MM_params['torsion'] != []:
        for element in MM_params['torsion']:
            for el in element[1]:
                MM_p.append(el)
                b_top.append( 30)
                b_low.append(-30)

                # if element[1][0] == 0:
                #     b_top.append( 0.00000001)
                #     b_low.append(-0.00000001)
                # b_top.append( 30)
                # b_low.append(-30)

    return MM_p, [b_low, b_top]

def k_params(MM_params):

    MM_p = []
    b_top = [] 
    b_low = []
    if MM_params['bond'] != []:
        for element in MM_params['bond']:
            MM_p.append(element[1][0])
            b_top.append(np.inf)
            b_low.append(0.5)
    if MM_params['angle'] != []:
        for element in MM_params['angle']:
            MM_p.append(element[1][0])
            b_top.append(np.inf)
            b_low.append(10)
    if MM_params['strbnd'] != []:
        for element in MM_params['strbnd']:
            if same_atom(element[0]):
                MM_p.append(element[1][0])
                b_top.append(np.inf)
                b_low.append(-np.inf)
            elif same_atom(element[0]) == False:
                MM_p.append(element[1][0])
                b_top.append(np.inf)
                b_low.append(-np.inf)
                MM_p.append(element[1][1])
                b_top.append(np.inf)
                b_low.append(-np.inf)
    if MM_params['improper'] != []:
        for element in MM_params['improper']:
            MM_p.append(element[1][0])
            b_top.append(np.inf)
            b_low.append(0)
    if MM_params['torsion'] != []:
        for element in MM_params['torsion']:
            for el in element[1]:
                MM_p.append(el)
                if element[1][0] == 0:
                    b_top.append( 0.00000001)
                    b_low.append(-0.00000001)
                else:
                    b_top.append( 30)
                    b_low.append(-30)

    return MM_p, [b_low, b_top]

def k_params_tor(MM_params):

    MM_p = []
    b_top = [] 
    b_low = []
    if MM_params['torsion'] != []:
        for element in MM_params['torsion']:
            MM_p.append(element[1][0])
            if element[1][0] == 0:
                b_top.append( 0.00000001)
                b_low.append(-0.00000001)
            b_top.append( 30)
            b_low.append(-30)

    return MM_p, [b_low, b_top]

def writekey(old_keyfile, new_keyfile, old_MM_param, param, which):
    shutil.copy(new_keyfile, new_keyfile+'throwaway')
    with open(new_keyfile+'throwaway', 'r') as in_file, open(new_keyfile, 'w') as out_file, open('zzz_watch_param.out', 'w') as watch_file:
        write = True
        for line in in_file:
            if 'EXTERNAL-FIELD' in line:
                continue
            if 'PARAMETERS TO OPTIMIZE' in line:
                write = False
            if 'DO NOT CHANGE' in line:
                write = True
            if 'RESTRAIN' in line:
                # print('found  ' + line)
                if 'vibrate' in which:
                    if '#' in line:
                        out_file.write(line)
                    else:
                        out_file.write('#' + line)
                    continue
                if 'minimize' in which:
                    if '#' in line:
                        out_file.write(line.replace('#',''))
                    else:
                        out_file.write(line)
                    continue
            elif write == True:
                out_file.write(line)

        if any('-50x' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD -50.0   0.0   0.0 \n')
        if any('-50y' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD   0.0 -50.0   0.0 \n')
        if any('-50z' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD   0.0   0.0 -50.0  \n')
        if any('-100z' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD   0.0   0.0 -100.0  \n')
        if any('+50x' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD  50.0   0.0   0.0 \n')
        if any('+50y' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD   0.0  50.0   0.0 \n')
        if any('+50z' in ele for ele in which):
            out_file.write('EXTERNAL-FIELD   0.0   0.0  50.0 \n')


        out_file.write('# PARAMETERS TO OPTIMIZE \n')
        if 'fit_xyz' in which:
            if 'vdw' in old_MM_param:
                length_vdw = len(old_MM_param['vdw'])
                for i in range(0, length_vdw):
                    string_1 = old_MM_param['vdw'][i][0]
                    string_2 = old_MM_param['vdw'][i][1][0]
                    string_3 = old_MM_param['vdw'][i][1][1]
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                length_vdw = len(old_MM_param['vdw'])*2
            else:
                length_vdw = 0
            if 'bond' in old_MM_param:
                length_bond = len(old_MM_param['bond'])
                for i in range(0, length_bond):
                    string_1 = old_MM_param['bond'][i][0]
                    string_2 = str(old_MM_param['bond'][i][1][0])
                    string_3 = str(param[i+length_vdw])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
            else:
                length_bond = 0
            if 'angle' in old_MM_param:
                length_angle = len(old_MM_param['angle'])
                for i in range(0, length_angle):
                    string_1 = old_MM_param['angle'][i][0]
                    string_2 = str(old_MM_param['angle'][i][1][0])
                    string_3 = str(param[i+length_vdw+length_bond])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
            else:
                length_angle = 0
            if 'strbnd' in old_MM_param:
                length_strbnd = len(old_MM_param['strbnd'])                   
                for i in range(0, length_strbnd):
                    string_1 = old_MM_param['strbnd'][i][0]
                    string_2 = str(old_MM_param['strbnd'][i][1][0])
                    string_3 = str(old_MM_param['strbnd'][i][1][0])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
            else:
                length_strbnd = 0      
            if 'torsion' in old_MM_param:
                length_torsion = len(old_MM_param['torsion'])
                for i in range(0, length_torsion):
                    string_1 = old_MM_param['torsion'][i][0]
                    string_2 = str(param[  3*i+length_vdw+length_bond+length_angle])
                    angle_2 = str(old_MM_param['torsion'][i][2][0])
                    string_3 = str(param[1+3*i+length_vdw+length_bond+length_angle])
                    angle_3 = str(old_MM_param['torsion'][i][2][1])
                    string_4 = str(param[2+3*i+length_vdw+length_bond+length_angle])
                    angle_4 = str(old_MM_param['torsion'][i][2][2])
                    # string_5 = str(param[3+i+length_vdw+length_bond+length_angle])
                    # angle_5 = str(old_MM_param['angle'][i][2][3])
                    # string_6 = str(param[4+i+length_vdw+length_bond+length_angle])
                    # angle_6 = str(old_MM_param['angle'][i][2][4])
                    # string_7 = str(param[5+i+length_vdw+length_bond+length_angle])
                    # angle_7 = str(old_MM_param['angle'][i][2][5])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' +'\n') #\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' + '\n')#\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
            else:
                length_torsion = 0  
                
        if 'fit_freq' in which:
            if 'vdw' in old_MM_param:
                length_vdw = len(old_MM_param['vdw'])
                for i in range(0, length_vdw):
                    string_1 = old_MM_param['vdw'][i][0]
                    string_2 = str(param[2*i])
                    string_3 = str(param[2*i+1])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                length_vdw = len(old_MM_param['vdw'])*2
            else:
                length_vdw = 0
            if 'bond' in old_MM_param:
                length_bond = len(old_MM_param['bond'])
                for i in range(0, length_bond):
                    string_1 = old_MM_param['bond'][i][0]
                    string_2 = str(param[i+length_vdw])
                    string_3 = str(old_MM_param['bond'][i][1][1])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
            else:
                length_bond = 0
            if 'angle' in old_MM_param:
                length_angle = len(old_MM_param['angle'])
                for i in range(0, length_angle):
                    string_1 = old_MM_param['angle'][i][0]
                    string_2 = str(param[i+length_vdw+length_bond])
                    string_3 = str(old_MM_param['angle'][i][1][1])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
            else:
                length_angle = 0
            if 'strbnd' in old_MM_param:
                length_strbnd = len(old_MM_param['strbnd'])
                
                for ele in old_MM_param['strbnd']:
                    if same_atom(ele[0]) == False:
                        length_strbnd = length_strbnd + 1
                length_strbnd = int(length_strbnd)
                
                i = 0
                i2 = 0
                for ele in old_MM_param['strbnd']:
                    if same_atom(ele[0]):
                        string_1 = old_MM_param['strbnd'][i][0]
                        string_2 = str(param[i+length_vdw+length_bond+length_angle])
                        string_3 = str(param[i+length_vdw+length_bond+length_angle])
                        i2 = i2 + 1
                        out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                        watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                        i = i + 1
                        
                    elif same_atom(ele[0]) == False:
                        string_1 = old_MM_param['strbnd'][i][0]
                        string_2 = str(param[i2+length_vdw+length_bond+length_angle])
                        i2 = i2 + 1
                        string_3 = str(param[i2+length_vdw+length_bond+length_angle])
                        i2 = i2 + 1
                        out_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                        watch_file.write(string_1 + ' ' + string_2 + ' ' + string_3 + '\n')
                        i = i + 1

            else:
                length_strbnd = 0              
            if 'torsion' in old_MM_param:
                length_torsion = len(old_MM_param['torsion'])
                for i in range(0, length_torsion):
                    string_1 = old_MM_param['torsion'][i][0]
                    string_2 = str(param[  3*i+length_vdw+length_bond+length_angle+length_strbnd])
                    angle_2 = str(old_MM_param['torsion'][i][2][0])
                    string_3 = str(param[1+3*i+length_vdw+length_bond+length_angle+length_strbnd])
                    angle_3 = str(old_MM_param['torsion'][i][2][1])
                    string_4 = str(param[2+3*i+length_vdw+length_bond+length_angle+length_strbnd])
                    angle_4 = str(old_MM_param['torsion'][i][2][2])
                    # string_5 = str(param[3+i+length_vdw+length_bond+length_angle+length_strbnd])
                    # angle_5 = str(old_MM_param['angle'][i][2][3])
                    # string_6 = str(param[4+i+length_vdw+length_bond+length_angle+length_strbnd])
                    # angle_6 = str(old_MM_param['angle'][i][2][4])
                    # string_7 = str(param[5+i+length_vdw+length_bond+length_angle+length_strbnd])
                    # angle_7 = str(old_MM_param['angle'][i][2][5])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' +'\n') #\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' + '\n')#\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
            else:
                length_torsion = 0     
                
        if 'pucker_tor' in which:
            if 'torsion' in old_MM_param:
                length_torsion = len(old_MM_param['torsion'])
                for i in range(0, length_torsion):
                    string_1 = old_MM_param['torsion'][i][0]
                    string_2 = str(param[  3*i])
                    angle_2 = str(old_MM_param['torsion'][i][2][0])
                    string_3 = str(param[1+3*i])
                    angle_3 = str(old_MM_param['torsion'][i][2][1])
                    string_4 = str(param[2+3*i])
                    angle_4 = str(old_MM_param['torsion'][i][2][2])
                    # string_5 = str(param[3+i])
                    # angle_5 = str(old_MM_param['angle'][i][2][3])
                    # string_6 = str(param[4+i])
                    # angle_6 = str(old_MM_param['angle'][i][2][4])
                    # string_7 = str(param[5+i])
                    # angle_7 = str(old_MM_param['angle'][i][2][5])
                    out_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' +'\n') #\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
                    watch_file.write(string_1 + ' ' + string_2 + ' ' + angle_2 + ' ' + '1' + ' ' +\
                                                    string_3 + ' ' + angle_3 + ' ' + '2' + ' ' +\
                                                    string_4 + ' ' + angle_4 + ' ' + '3' + ' ' + '\n')#\
                                                    # string_5 + ' ' + angle_5 + ' ' + '4' + ' ' +\
                                                    # string_6 + ' ' + angle_6 + ' ' + '5' + ' ' +\
                                                    # string_7 + ' ' + angle_7 + ' ' + '6' + ' ' + '\n')
                
        watch_file.flush()
        out_file.flush()
        os.remove(new_keyfile+'throwaway')
        
def rewritekey(keyfile, list_change, list_notchange):
    with open(keyfile, 'r') as in_file, open(keyfile+'_1', 'w') as out_file_1, open(keyfile+'_2', 'w') as out_file_2:
        param_dict = {'bond' : [],
                      'angle': [],
                      'strbnd': [],
                      'torsion': []}       
        donotchange_params = True

        for line in in_file:
            if 'PARAMETERS TO OPTIMIZE' in line:
                donotchange_params = False
            if 'DO NOT CHANGE' in line:
                donotchange_params = True     
            if donotchange_params == True:
                if 'bond' in line:
                    param_dict['bond'].append(line)
                elif 'angle' in line:
                    param_dict['angle'].append(line)
                elif 'strbnd' in line:
                    param_dict['strbnd'].append(line)
                else:
                    out_file_1.write(line)   
            if donotchange_params == False:
                if 'bond' in line:
                    param_dict['bond'].append(line)
                elif 'angle' in line:
                    param_dict['angle'].append(line)
                elif 'strbnd' in line:
                    param_dict['strbnd'].append(line)
                elif 'torsion' in line:
                    param_dict['torsion'].append(line)
                else:
                    out_file_2.write(line)  
                    
        for param in list_notchange:
            if param_dict[param] != []:
                for line in param_dict[param]:
                    out_file_1.write(line)
        for param in list_change:
            if param_dict[param] != []:
                for line in param_dict[param]:
                    out_file_2.write(line)

            
        out_file_1.flush()
        out_file_2.flush()

                    
    with open(keyfile, 'w') as out_file, open(keyfile+'_1', 'r') as in_file_1, open(keyfile+'_2', 'r') as in_file_2:
        for line in in_file_1:
            out_file.write(line) 
        for line in in_file_2:
            out_file.write(line) 
        out_file.flush()

    # os.remove(keyfile+'_1')
    # os.remove(keyfile+'_2')
        
def add_restrainplane(keyfile,additional):
    shutil.copy(keyfile, keyfile+'throwaway')
    
    
    with open(keyfile+'throwaway', 'r') as in_file, open(keyfile, 'w') as out_file:
        for line in in_file:
            if 'RESTRAIN-PLANE' in line:
                return True

    with open(keyfile+'throwaway', 'r') as in_file, open(keyfile, 'w') as out_file:
        write = False
        write2 = True
        for line in in_file:
            if 'atom' in line:
                write = True
            if write == True and write2 == True:
                out_file.write('RESTRAIN-PLANE X 1 0.0 \n')
                out_file.write('RESTRAIN-PLANE Y 1 0.0 \n')
                out_file.write('RESTRAIN-PLANE Z 1 0.0 \n')
                out_file.write('RESTRAIN-PLANE X 2 0.0 \n')
                out_file.write('RESTRAIN-PLANE Y 2 0.0 \n')
                out_file.write('RESTRAIN-PLANE Y 3 0.0 \n')
                if additional != []:
                    for element in additional:
                        out_file.write(element + ' \n')
                out_file.write('\n')
                write2 = False
            out_file.write(line)
                
    os.remove(keyfile+'throwaway')