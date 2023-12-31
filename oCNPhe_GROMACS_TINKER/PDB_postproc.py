# Version 2023/09/19

import numpy as np

################################################################
# INPUT:

    

in_file_prefix = 'PYP_md_out'
C_atom = 406
N_atom = 407

# in_file_prefix = 'zzz_PYP_F92oCNF_md'
# C_atom = 1361
# N_atom = 1362

water_cutoff = 0.1 # threshold wrt to maximum probability to cutoff
water_grid_step = 0.5 #in Angstroms
water_grid_size = 10 # for 10 A, it will make grid of  [-5...+5],[-5...+5],[-5...+5]
water_grid_unit = 2 # with 1 A, volumes of 1 x 1 x 1 A^3 will be evaluated; this should be >= 2*water_grid_step

maxHBdist  = 3.5
maxHBangle = 30

debug_mode = False
debug_frames = 1000

    
################################################################

in_pdb = in_file_prefix+'.pdb'
out_pdb_first  = '01_'+in_file_prefix+'_first.pdb'
out_pdb_prot   = '01_'+in_file_prefix+'_prot_separated.pdb'
out_pdb_wat    = '01_'+in_file_prefix+'_wat_separated.pdb'

out_hbond_prot = '02_'+in_file_prefix+'_hbond_prot.out'
out_hbond_wat  = '02_'+in_file_prefix+'_hbond_wat.out'
out_hbond_none = '02_'+in_file_prefix+'_hbond_none.out'
out_hbond_prot_wo = '02_'+in_file_prefix+'_hbond_prot_wo.out'
out_hbond_wat_wo  = '02_'+in_file_prefix+'_hbond_wat_wo.out'
out_hbond_none_wo = '02_'+in_file_prefix+'_hbond_none_wo.out'
out_hbond_prot_info = '02_'+in_file_prefix+'_hbond_prot_info.out'
out_hbond_wat_info  = '02_'+in_file_prefix+'_hbond_wat_info.out'
out_hbond_none_info = '02_'+in_file_prefix+'_hbond_none_info.out'


out_pdb_none_only   = '11_'+in_file_prefix+'_none_only.pdb'
out_pdb_prot_only   = '11_'+in_file_prefix+'_prot_only.pdb'
out_pdb_wat_only    = '11_'+in_file_prefix+'_wat_only.pdb'
out_pdb_prot_wat    = '11_'+in_file_prefix+'_prot_wat.pdb'

out_hbond_prot_only = '12_'+in_file_prefix+'_hbond_prot_only.out'
out_hbond_wat_only  = '12_'+in_file_prefix+'_hbond_wat_only.out'
out_hbond_prot_wat  = '12_'+in_file_prefix+'_hbond_prot_wat.out'
out_hbond_prot_only_wo = '12_'+in_file_prefix+'_hbond_prot_only_wo.out'
out_hbond_wat_only_wo  = '12_'+in_file_prefix+'_hbond_wat_only_wo.out'
out_hbond_prot_wat_wo  = '12_'+in_file_prefix+'_hbond_prot_wat_wo.out'
out_hbond_prot_only_info = '12_'+in_file_prefix+'_hbond_prot_only_info.out'
out_hbond_wat_only_info  = '12_'+in_file_prefix+'_hbond_wat_only_info.out'
out_hbond_prot_wat_info  = '12_'+in_file_prefix+'_hbond_prot_wat_info.out'


out_pdb_all_avg     = '20_'+in_file_prefix+'_all_avg.pdb'
out_pdb_all_wat     = '20_'+in_file_prefix+'_all_avg_H2O.pdb'
out_pdb_all_rmsf    = '20_'+in_file_prefix+'_all_rmsf.pdb'
out_pdb_all_Bfactor = '20_'+in_file_prefix+'_all_avg_Bfactor.pdb'
out_pdb_all_Zscore  = '20_'+in_file_prefix+'_all_avg_Zscore.pdb'

out_pdb_none_avg     = '21_'+in_file_prefix+'_none_avg.pdb'
out_pdb_none_wat     = '21_'+in_file_prefix+'_none_avg_H2O.pdb'
out_pdb_none_rmsf    = '21_'+in_file_prefix+'_none_rmsf.pdb'
out_pdb_none_Bfactor = '21_'+in_file_prefix+'_none_avg_Bfactor.pdb'
out_pdb_none_Zscore  = '21_'+in_file_prefix+'_none_avg_Zscore.pdb'

out_pdb_prot_avg     = '22_'+in_file_prefix+'_prot_avg.pdb'
out_pdb_prot_rmsf    = '22_'+in_file_prefix+'_prot_rmsf.pdb'
out_pdb_prot_wat     = '22_'+in_file_prefix+'_prot_avg_H2O.pdb'
out_pdb_prot_Bfactor = '22_'+in_file_prefix+'_prot_avg_Bfactor.pdb'
out_pdb_prot_Zscore  = '22_'+in_file_prefix+'_prot_avg_Zscore.pdb'

out_pdb_wat_avg      = '23_'+in_file_prefix+'_wat_avg.pdb'
out_pdb_wat_wat      = '23_'+in_file_prefix+'_wat_avg_H2O.pdb'
out_pdb_wat_rmsf     = '23_'+in_file_prefix+'_wat_rmsf.pdb'
out_pdb_wat_Bfactor  = '23_'+in_file_prefix+'_wat_avg_Bfactor.pdb'
out_pdb_wat_Zscore   = '23_'+in_file_prefix+'_wat_avg_Zscore.pdb'

out_pdb_prot_wat_avg     = '24_'+in_file_prefix+'_prot_wat_avg.pdb'
out_pdb_prot_wat_wat     = '24_'+in_file_prefix+'_prot_wat_avg_H2O.pdb'
out_pdb_prot_wat_rmsf    = '24_'+in_file_prefix+'_prot_wat_rmsf.pdb'
out_pdb_prot_wat_Bfactor = '24_'+in_file_prefix+'_prot_wat_avg_Bfactor.pdb'
out_pdb_prot_wat_Zscore  = '24_'+in_file_prefix+'_prot_wat_avg_Zscore.pdb'

chain_info = True
if chain_info == False:
    corr = -1
else:
    corr = 0

print('INFO: Separating Protein from Water for separate PDB files.')
C_coords = []
N_coords = []
Dp_coords = []
Hp_coords = []
Dw_coords = []
Hw_coords = []
Dp_identity = []
Dw_identity = []
count = -1
with open(in_pdb, 'r') as in_file, open(out_pdb_prot, 'w') as out_prot, open(out_pdb_wat, 'w') as out_wat, open(out_pdb_first, 'w') as out_first:
    for line in in_file:
        temp_CN = []

        if 'MODEL' in line:
            count += 1
            out_wat.write(line)
        if debug_mode == True:
            if count > debug_frames:
                break
        if count == 0:
            out_first.write(line)
        if 'OCF' in line:
            out_wat.write(line)
        if 'SOL' in line:
            out_wat.write(line)
        else:
            out_prot.write(line)
            
        if 'OCF' in line:
            if line.split()[1] == str(C_atom):
                C_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])
            elif line.split()[1] == str(N_atom):
                N_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])     
        if 'SOL' in line:
            if line.split()[2] in ['O', 'OW']:
                Dw_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])
                Dw_identity.append([count, line[0:27]])
            if line.split()[2] in ['H', 'HW1', 'HW2']:
                Hw_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])     
        elif 'MODEL' not in line and 'END' not in line and line.split()[1] != str(N_atom):
            if line.split()[2] in ['O', 'OH', 'N', 'S', 'OG', 'ND1', 'ND2', 'OG1', 'NZ', 'NE2', 'OE1', 'OE2', 'NH1', 'NH2']:
                Dp_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])
                Dp_identity.append([count, line[0:27]])
            if line.split()[2] in ['H', 'HN', 'HO', 'HD1', 'HE2', 'HG', '1HD2', '2HD2', 'HG1', 'HZ1', 'HZ2', 'HZ3', '1HE2', '2HE2', 'HH', 'HE2', '1HH1', '1HH2', '2HH1', '2HH2']:
                Hp_coords.append([count, [float(x) for x in line.split()[6+corr:9+corr]]])      
     
print('INFO: Analyzing H-bonds.')
C_list = []
N_list = []
which_Dp_list = []
which_Dw_list = []
Dp_list = []
Dw_list = []
Hp_list = []
Hw_list = []
for el in C_coords:
    Dp_list.append([])
    Dw_list.append([])
    Hp_list.append([])
    Hw_list.append([])
    which_Dp_list.append([])
    which_Dw_list.append([])

for C_el, N_el in zip(C_coords, N_coords):
    C_list.append(C_el[1])
    N_list.append(N_el[1])

def if_within(atom, ref_atom, max_dist):
    
    if max_dist == []:
        max_dist = 3.5
    dist_vec = np.array(atom) - np.array(ref_atom)
    dist = np.sqrt(np.dot(dist_vec, dist_vec))
    if dist <= max_dist:
        return True
    
for Dp_el,Dp_id in zip(Dp_coords,Dp_identity):
    frame = Dp_el[0]
    # print(Dp_el)
    if len(Dp_el) == 1:
        # Dp_list[frame].append([])
        print('No element in frame ' + str(frame))
    elif if_within(N_list[frame], Dp_el[1], maxHBdist):
        Dp_list[frame].append(Dp_el[1])
        # print(Dp_id[1])
        which_Dp_list[frame].append(Dp_id[1])
        
for Dw_el,Dw_id in zip(Dw_coords,Dw_identity):
    frame = Dw_el[0]
    if len(Dw_el) == 1:
        # Dw_list[frame].append([])
        print('No element in frame ' + str(frame))
    elif if_within(N_list[frame], Dw_el[1], maxHBdist):
        Dw_list[frame].append(Dw_el[1])
        # print(Dw_id[1])
        which_Dw_list[frame].append(Dw_id[1])
        
for Hp_el in Hp_coords:
    frame = Hp_el[0]
    if len(Hp_el) == 1:
        # Dp_list[frame].append([])
        print('No element in frame ' + str(frame))
    elif if_within(N_list[frame], Hp_el[1], maxHBdist):
        Hp_list[frame].append(Hp_el[1])
        
for Hw_el in Hw_coords:
    frame = Hw_el[0]
    if len(Hw_el) == 1:
        # Dw_list[frame].append([])
        print('No element in frame ' + str(frame))
    elif if_within(N_list[frame], Hw_el[1], maxHBdist):
        Hw_list[frame].append(Hw_el[1])
    
# print(which_Dp_list)
# print(which_Dw_list)   

def if_hbond(donor, acceptor, hydrogen, max_dist, max_angle):
    if max_dist == []:
        max_dist = 3.5
    if max_angle == []:
        max_angle = 30
    DA_vec = np.array(donor) - np.array(acceptor)
    DA_dist = np.sqrt(np.dot(DA_vec, DA_vec))
    norm_DA_vec = DA_vec/DA_dist
    
    HA_vec = np.array(hydrogen) - np.array(acceptor) 
    HA_dist = np.sqrt(np.dot(HA_vec, HA_vec))    
    # norm_HA_vec = HA_vec/HA_dist
    
    DH_vec = np.array(donor) - np.array(hydrogen) 
    DH_dist = np.sqrt(np.dot(DH_vec, DH_vec))    
    norm_DH_vec = DH_vec/DH_dist
    
    angle = np.arccos(np.dot(norm_DA_vec,norm_DH_vec))*180/np.pi
    
    if DH_dist <= 1.8:
        if HA_dist <= DA_dist:
            if abs(angle) <= max_angle:
                return True   
            
    return False
        
def dist_angle(C, N, D, H):

    CN_vec = np.array(N) - np.array(C)
    CN_dist = np.sqrt(np.dot(CN_vec, CN_vec))    
    norm_CN_vec = CN_vec/CN_dist
    
    DN_vec = np.array(N) - np.array(D)
    DN_dist = np.sqrt(np.dot(DN_vec, DN_vec))    
    norm_DN_vec = DN_vec/DN_dist
    
    angle = np.arccos(np.dot(norm_CN_vec,norm_DN_vec))*180/np.pi
    
    return DN_dist, angle
    
    
# collect info and write

hbond_water = []

for C, N, D, H, D_id in zip(C_list, N_list, Dw_list, Hw_list, which_Dw_list):
    temp_hbond = []
    if D == []:
        hbond_water.append(False)
        continue
    elif H == []:
        hbond_water.append(False)
        continue
    for D_sub, D_id_sub in zip(D, D_id):
        for H_sub in H:   
            # print(C)
            # print(N)
            # print(D_sub)
            # print(H_sub)
            if if_hbond(D_sub, N, H_sub, maxHBdist, maxHBangle):
                dist, angle = dist_angle(C, N, D_sub, H_sub)
                temp_hbond.append([round(dist,3), round(angle,3), D_id_sub])

    if temp_hbond == []:
        hbond_water.append(False)
    else:
        hbond_water.append(temp_hbond)
                
# print(hbond_water)

hbond_protein = []

for C, N, D, H, D_id in zip(C_list, N_list, Dp_list, Hp_list, which_Dp_list):
    temp_hbond = []
    if D == []:
        hbond_protein.append(False)
        continue
    elif H == []:
        hbond_protein.append(False)
        continue
    for D_sub, D_id_sub in zip(D, D_id):
        for H_sub in H:   
            # print(C)
            # print(N)
            # print(D_sub)
            # print(H_sub)
            if if_hbond(D_sub, N, H_sub, maxHBdist, maxHBangle):
                dist, angle = dist_angle(C, N, D_sub, H_sub)
                temp_hbond.append([round(dist,3), round(angle,3), D_id_sub])

    if temp_hbond == []:
        hbond_protein.append(False)
    else:
        hbond_protein.append(temp_hbond)
                
# print(hbond_protein)
# collect info [prot, wat]
hbond_info = []
with open(out_hbond_prot, 'w') as prot_out, open(out_hbond_wat, 'w') as wat_out, open(out_hbond_none, 'w') as none_out, open(out_hbond_prot_only, 'w') as prot_only_out, open(out_hbond_wat_only, 'w') as wat_only_out, open(out_hbond_prot_wat, 'w') as prot_wat_out,\
     open(out_hbond_prot_wo, 'w') as prot_out_wo, open(out_hbond_wat_wo, 'w') as wat_out_wo, open(out_hbond_none_wo, 'w') as none_out_wo, open(out_hbond_prot_only_wo, 'w') as prot_only_out_wo, open(out_hbond_wat_only_wo, 'w') as wat_only_out_wo, open(out_hbond_prot_wat_wo, 'w') as prot_wat_out_wo,\
     open(out_hbond_prot_info, 'w') as prot_out_info, open(out_hbond_wat_info, 'w') as wat_out_info, open(out_hbond_none_info, 'w') as none_out_info, open(out_hbond_prot_only_info, 'w') as prot_only_out_info, open(out_hbond_wat_only_info, 'w') as wat_only_out_info, open(out_hbond_prot_wat_info, 'w') as prot_wat_out_info:
    frame = 1
    for prot, wat in zip(hbond_protein,hbond_water):
        
        if prot != False:
            for prot_sub in prot:
                prot_out.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + '    ' + prot_sub[2]  + ' \n')
                prot_out_wo.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + ' \n')
                prot_out_info.write(str(frame).rjust(6) + ' ' + prot_sub[2]  + ' \n')
        if wat != False:
            for wat_sub in wat:
                wat_out.write(str(frame).rjust(6) + ' ' + str(wat_sub[0]).rjust(8) + ' ' + str(wat_sub[1]).rjust(8) + '    ' + wat_sub[2] + ' \n')
                wat_out_wo.write(str(frame).rjust(6) + ' ' + str(wat_sub[0]).rjust(8) + ' ' + str(wat_sub[1]).rjust(8) + ' \n')
                wat_out_info.write(str(frame).rjust(6) + ' ' + wat_sub[2] + ' \n')
        if prot == False and wat == False:
            none_out.write(str(frame) + ' None  \n')
            none_out_wo.write(str(frame) + ' \n')
            hbond_info.append([False, False])
        if prot != False and wat == False:
            for prot_sub in prot:
                prot_only_out.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + '    ' + prot_sub[2] + ' \n')
                prot_only_out_wo.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + ' \n')
                prot_only_out_info.write(str(frame).rjust(6) + ' ' + prot_sub[2] + ' \n')
            hbond_info.append([True, False])
        if prot == False and wat != False:
            for wat_sub in wat:
                wat_only_out.write(str(frame).rjust(6) + ' ' + str(wat_sub[0]).rjust(8) + ' ' + str(wat_sub[1]).rjust(8) + '    ' + wat_sub[2] + ' \n')
                wat_only_out_wo.write(str(frame).rjust(6) + ' ' + str(wat_sub[0]).rjust(8) + ' ' + str(wat_sub[1]).rjust(8) + ' \n')
                wat_only_out_info.write(str(frame).rjust(6) + ' ' + wat_sub[2] + ' \n')
            hbond_info.append([False, True]) 
        if prot != False and wat != False:
            for prot_sub in prot:
                prot_wat_out.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + '    ' + prot_sub[2] + ' \n')
                prot_wat_out_wo.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + ' \n')
                prot_wat_out_info.write(str(frame).rjust(6) + ' ' + prot_sub[2] + ' \n')
            for wat_sub in wat:
                prot_wat_out.write(str(frame).rjust(6) + ' ' + str(wat_sub[0]).rjust(8) + ' ' + str(wat_sub[1]).rjust(8) + '    ' + wat_sub[2] + ' \n')
                prot_wat_out_wo.write(str(frame).rjust(6) + ' ' + str(prot_sub[0]).rjust(8) + ' ' + str(prot_sub[1]).rjust(8) + ' \n')
                prot_wat_out_info.write(str(frame).rjust(6) + ' ' + prot_sub[2] + ' \n')
            hbond_info.append([True, True])  
        
        frame += 1
                    
# print(hbond_info)
            
count = -1
with open(in_pdb, 'r') as in_file, open(out_pdb_none_only, 'w') as out_none, open(out_pdb_prot_only, 'w') as out_prot, open(out_pdb_wat_only, 'w') as out_wat, open(out_pdb_prot_wat, 'w') as out_prot_wat:
    for line in in_file:

        if 'MODEL' in line:
            count += 1
            out_none.write(line)
            out_prot.write(line)
            out_wat.write(line)
            out_prot_wat.write(line)
        if debug_mode == True:
            if count > debug_frames:
                break
            
        if hbond_info[count] == [False, False]:
            out_none.write(line)
        elif hbond_info[count] == [True, False]:
            out_prot.write(line)            
        elif hbond_info[count] == [False, True]:
            out_wat.write(line)  
        elif hbond_info[count] == [True, True]:
            out_prot_wat.write(line)  
                
     
print('INFO: Generating average PDB files and RMSFs.')
# AVERAGING all
##############################################################################
temp_for_avg = []
n_to_divide_by = 0
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append([0,0,0])
            
# print(temp_for_avg)
    
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_first, 'r') as in_first, open(out_pdb_all_avg, 'w') as out_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
            n_to_divide_by += 1

        if 'MODEL' not in line:
              if 'TER' not in line:
                try:                      
                    atom = int(line.split()[1])-1
                    # print(line)
                    temp_for_avg[atom][0] = temp_for_avg[atom][0] + float(line.split()[6])
                    temp_for_avg[atom][1] = temp_for_avg[atom][1] + float(line.split()[7])
                    temp_for_avg[atom][2] = temp_for_avg[atom][2] + float(line.split()[8])
                except:
                    print('WARNING: atom number out of range')     
     
    avg_list = np.array(temp_for_avg)/n_to_divide_by
    for line in in_first:
        line_out = []
        if 'MODEL' in line:
            out_file.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                line_out.append((line.split()[0]).ljust(6))
                line_out.append((line.split()[1]).rjust(5))
                line_out.append('  ' + (line.split()[2]).ljust(4))
                line_out.append((line.split()[3]).ljust(3))
                line_out.append((line.split()[4]).rjust(2))
                line_out.append((line.split()[5]).rjust(4))
                line_out.append('    ' + str(round(avg_list[atom][0],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][1],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][2],3)).rjust(8))    
                line_out.append((line.split()[9]).rjust(6))
                line_out.append((line.split()[10]).rjust(6))
                out_file.write(''.join(x for x in line_out) + '\n')
            except:
                print('WARNING in line: ' + line)
                
temp_for_avg = []
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append(0)
print(avg_list)
            
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_all_avg, 'r') as in_avg, open(out_pdb_all_rmsf, 'w') as out_file, open(out_pdb_all_Bfactor, 'w') as out_file_B, open(out_pdb_all_Zscore, 'w') as out_file_Z:
    for line in in_file:
        if 'MODEL' in line:
            count += 1    

        if 'MODEL' not in line:
              if 'TER' not in line:
                try:                  
                    atom = int(line.split()[1])-1                  
                    dist_vec = np.array(avg_list[atom]) - np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                    dist2 = np.dot(dist_vec,dist_vec)
                    # print(dist2)
                    temp_for_avg[atom] = temp_for_avg[atom] + dist2
                except:
                    print('WARNING: atom number out of range')

    avg_list = np.array(temp_for_avg)/n_to_divide_by
    RMSF_list = np.sqrt(avg_list)
    B_factor_list = [round(x,2) for x in (8*np.pi**2*RMSF_list**2/3)]
    B_factor_avg = np.mean(B_factor_list)
    B_factor_std = np.std(B_factor_list)
    Zscore_list = [round(x,2) for x in (B_factor_list - B_factor_avg)/B_factor_std]
    
    for line in in_avg:
        line_out = []
        line_out_B = []
        line_out_Z = []
        if 'MODEL' in line:
            out_file.write(line)
            out_file_B.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                # line_out.append((line.split()[0]).ljust(6))
                # line_out.append((line.split()[1]).rjust(5))
                # line_out.append('  ' + (line.split()[2]).ljust(4))
                # line_out.append((line.split()[3]).ljust(3))
                # line_out.append((line.split()[4]).rjust(2))
                # line_out.append((line.split()[5]).rjust(4))
                line_out.append(line[0:60])
                line_out.append(str(round(RMSF_list[atom],2)).rjust(7))
                out_file.write(' '.join(x for x in line_out) + '\n')
                
                line_out_B.append(line[0:60])
                line_out_B.append(str(B_factor_list[atom]).rjust(7))
                out_file_B.write(' '.join(x for x in line_out_B) + '\n')
                
                line_out_Z.append(line[0:60])
                line_out_Z.append(str(Zscore_list[atom]).rjust(7))
                out_file_Z.write(' '.join(x for x in line_out_Z) + '\n')
            except:
                print('WARNING in line: ' + line)      
                
N_frames_all = n_to_divide_by


# AVERAGING non-h-bonded
##############################################################################
temp_for_avg = []
n_to_divide_by = 0
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append([0,0,0])
            
# print(temp_for_avg)
    
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_first, 'r') as in_first, open(out_pdb_none_avg, 'w') as out_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
            if hbond_info[count] == [False, False]:
                n_to_divide_by += 1
        if hbond_info[count] == [False, False]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                try:                      
                    atom = int(line.split()[1])-1
                    # print(line)
                    temp_for_avg[atom][0] = temp_for_avg[atom][0] + float(line.split()[6])
                    temp_for_avg[atom][1] = temp_for_avg[atom][1] + float(line.split()[7])
                    temp_for_avg[atom][2] = temp_for_avg[atom][2] + float(line.split()[8])
                except:
                    print('WARNING: atom number out of range')     
     
    avg_list = np.array(temp_for_avg)/n_to_divide_by
    for line in in_first:
        line_out = []
        if 'MODEL' in line:
            out_file.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                line_out.append((line.split()[0]).ljust(6))
                line_out.append((line.split()[1]).rjust(5))
                line_out.append('  ' + (line.split()[2]).ljust(4))
                line_out.append((line.split()[3]).ljust(3))
                line_out.append((line.split()[4]).rjust(2))
                line_out.append((line.split()[5]).rjust(4))
                line_out.append('    ' + str(round(avg_list[atom][0],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][1],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][2],3)).rjust(8))    
                line_out.append((line.split()[9]).rjust(6))
                line_out.append((line.split()[10]).rjust(6))
                out_file.write(''.join(x for x in line_out) + '\n')
            except:
                print('WARNING in line: ' + line)
                
temp_for_avg = []
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append(0)
print(avg_list)
            
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_none_avg, 'r') as in_avg, open(out_pdb_none_rmsf, 'w') as out_file, open(out_pdb_none_Bfactor, 'w') as out_file_B, open(out_pdb_none_Zscore, 'w') as out_file_Z:
    for line in in_file:
        if 'MODEL' in line:
            count += 1    
        if hbond_info[count] == [False, False]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                try:                  
                    atom = int(line.split()[1])-1                  
                    dist_vec = np.array(avg_list[atom]) - np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                    dist2 = np.dot(dist_vec,dist_vec)
                    # print(dist2)
                    temp_for_avg[atom] = temp_for_avg[atom] + dist2
                except:
                    print('WARNING: atom number out of range')

    avg_list = np.array(temp_for_avg)/n_to_divide_by
    RMSF_list = np.sqrt(avg_list)
    B_factor_list = [round(x,2) for x in (8*np.pi**2*RMSF_list**2/3)]
    B_factor_avg = np.mean(B_factor_list)
    B_factor_std = np.std(B_factor_list)
    Zscore_list = [round(x,2) for x in (B_factor_list - B_factor_avg)/B_factor_std]
    
    for line in in_avg:
        line_out = []
        line_out_B = []
        line_out_Z = []
        if 'MODEL' in line:
            out_file.write(line)
            out_file_B.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                # line_out.append((line.split()[0]).ljust(6))
                # line_out.append((line.split()[1]).rjust(5))
                # line_out.append('  ' + (line.split()[2]).ljust(4))
                # line_out.append((line.split()[3]).ljust(3))
                # line_out.append((line.split()[4]).rjust(2))
                # line_out.append((line.split()[5]).rjust(4))
                line_out.append(line[0:60])
                line_out.append(str(round(RMSF_list[atom],2)).rjust(7))
                out_file.write(' '.join(x for x in line_out) + '\n')
                
                line_out_B.append(line[0:60])
                line_out_B.append(str(B_factor_list[atom]).rjust(7))
                out_file_B.write(' '.join(x for x in line_out_B) + '\n')
                
                line_out_Z.append(line[0:60])
                line_out_Z.append(str(Zscore_list[atom]).rjust(7))
                out_file_Z.write(' '.join(x for x in line_out_Z) + '\n')
            except:
                print('WARNING in line: ' + line)      
                
N_frames_none = n_to_divide_by

##############################################################################
# AVERAGING protein-h-bonded
temp_for_avg = []
n_to_divide_by = 0
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append([0,0,0])
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_first, 'r') as in_first, open(out_pdb_prot_avg, 'w') as out_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
            if hbond_info[count] == [True, False]:
                n_to_divide_by += 1
        if hbond_info[count] == [True, False]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                try:                      
                    atom = int(line.split()[1])-1
                    # print(line)
                    temp_for_avg[atom][0] = temp_for_avg[atom][0] + float(line.split()[6])
                    temp_for_avg[atom][1] = temp_for_avg[atom][1] + float(line.split()[7])
                    temp_for_avg[atom][2] = temp_for_avg[atom][2] + float(line.split()[8])
                except:
                    print('WARNING: atom number out of range')     
    # print(temp_for_avg)
    # print(n_to_divide_by)
    avg_list = np.array(temp_for_avg)/n_to_divide_by
    for line in in_first:
        line_out = []
        if 'MODEL' in line:
            out_file.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                line_out.append((line.split()[0]).ljust(6))
                line_out.append((line.split()[1]).rjust(5))
                line_out.append('  ' + (line.split()[2]).ljust(4))
                line_out.append((line.split()[3]).ljust(3))
                line_out.append((line.split()[4]).rjust(2))
                line_out.append((line.split()[5]).rjust(4))
                line_out.append('    ' + str(round(avg_list[atom][0],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][1],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][2],3)).rjust(8))    
                line_out.append((line.split()[9]).rjust(6))
                line_out.append((line.split()[10]).rjust(6))
                out_file.write(''.join(x for x in line_out) + '\n')
            except:
                print('WARNING in line: ' + line)    


temp_for_avg = []
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append(0)
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_prot_avg, 'r') as in_avg, open(out_pdb_prot_rmsf, 'w') as out_file, open(out_pdb_prot_Bfactor, 'w') as out_file_B, open(out_pdb_prot_Zscore, 'w') as out_file_Z:
    for line in in_file:
        if 'MODEL' in line:
            count += 1    
        if hbond_info[count] == [True, False]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                  
                try:                  
                    atom = int(line.split()[1])-1
                    dist_vec = np.array(avg_list[atom]) - np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                    dist2 = np.dot(dist_vec,dist_vec)
                    temp_for_avg[atom] = temp_for_avg[atom] + dist2
                except:
                    print('WARNING: atom number out of range')

    avg_list = np.array(temp_for_avg)/n_to_divide_by
    RMSF_list = np.sqrt(avg_list)
    B_factor_list = [round(x,2) for x in (8*np.pi**2*RMSF_list**2/3)]
    B_factor_avg = np.mean(B_factor_list)
    B_factor_std = np.std(B_factor_list)
    Zscore_list = [round(x,2) for x in (B_factor_list - B_factor_avg)/B_factor_std]
    
    for line in in_avg:
        line_out = []
        line_out_B = []
        line_out_Z = []
        if 'MODEL' in line:
            out_file.write(line)
            out_file_B.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                # line_out.append((line.split()[0]).ljust(6))
                # line_out.append((line.split()[1]).rjust(5))
                # line_out.append('  ' + (line.split()[2]).ljust(4))
                # line_out.append((line.split()[3]).ljust(3))
                # line_out.append((line.split()[4]).rjust(2))
                # line_out.append((line.split()[5]).rjust(4))
                line_out.append(line[0:60])
                line_out.append(str(round(RMSF_list[atom],2)).rjust(7))
                out_file.write(' '.join(x for x in line_out) + '\n')
                
                line_out_B.append(line[0:60])
                line_out_B.append(str(B_factor_list[atom]).rjust(7))
                out_file_B.write(' '.join(x for x in line_out_B) + '\n')
                
                line_out_Z.append(line[0:60])
                line_out_Z.append(str(Zscore_list[atom]).rjust(7))
                out_file_Z.write(' '.join(x for x in line_out_Z) + '\n')
            except:
                print('WARNING in line: ' + line)        

N_frames_prot = n_to_divide_by

##############################################################################
# AVERAGING wat-h-bonded
temp_for_avg = []
n_to_divide_by = 0
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append([0,0,0])    
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_first, 'r') as in_first, open(out_pdb_wat_avg, 'w') as out_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
            if hbond_info[count] == [False, True]:
                n_to_divide_by += 1
        if hbond_info[count] == [False, True]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                try:                      
                    atom = int(line.split()[1])-1
                    # print(line)
                    temp_for_avg[atom][0] = temp_for_avg[atom][0] + float(line.split()[6])
                    temp_for_avg[atom][1] = temp_for_avg[atom][1] + float(line.split()[7])
                    temp_for_avg[atom][2] = temp_for_avg[atom][2] + float(line.split()[8])
                except:
                    print('WARNING: atom number out of range')     
     
    avg_list = np.array(temp_for_avg)/n_to_divide_by
    for line in in_first:
        line_out = []
        if 'MODEL' in line:
            out_file.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                line_out.append((line.split()[0]).ljust(6))
                line_out.append((line.split()[1]).rjust(5))
                line_out.append('  ' + (line.split()[2]).ljust(4))
                line_out.append((line.split()[3]).ljust(3))
                line_out.append((line.split()[4]).rjust(2))
                line_out.append((line.split()[5]).rjust(4))
                line_out.append('    ' + str(round(avg_list[atom][0],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][1],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][2],3)).rjust(8))    
                line_out.append((line.split()[9]).rjust(6))
                line_out.append((line.split()[10]).rjust(6))
                out_file.write(''.join(x for x in line_out) + '\n')
            except:
                print('WARNING in line: ' + line)            
                
temp_for_avg = []
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append(0)
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_wat_avg, 'r') as in_avg, open(out_pdb_wat_rmsf, 'w') as out_file, open(out_pdb_wat_Bfactor, 'w') as out_file_B, open(out_pdb_wat_Zscore, 'w') as out_file_Z:
    for line in in_file:
        if 'MODEL' in line:
            count += 1    
        if hbond_info[count] == [False, True]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                  
                try:                  
                    atom = int(line.split()[1])-1
                    dist_vec = np.array(avg_list[atom]) - np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                    dist2 = np.dot(dist_vec,dist_vec)
                    temp_for_avg[atom] = temp_for_avg[atom] + dist2
                except:
                    print('WARNING: atom number out of range')

    avg_list = np.array(temp_for_avg)/n_to_divide_by
    RMSF_list = np.sqrt(avg_list)
    B_factor_list = [round(x,2) for x in (8*np.pi**2*RMSF_list**2/3)]
    B_factor_avg = np.mean(B_factor_list)
    B_factor_std = np.std(B_factor_list)
    Zscore_list = [round(x,2) for x in (B_factor_list - B_factor_avg)/B_factor_std]
    
    for line in in_avg:
        line_out = []
        line_out_B = []
        line_out_Z = []
        if 'MODEL' in line:
            out_file.write(line)
            out_file_B.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                # line_out.append((line.split()[0]).ljust(6))
                # line_out.append((line.split()[1]).rjust(5))
                # line_out.append('  ' + (line.split()[2]).ljust(4))
                # line_out.append((line.split()[3]).ljust(3))
                # line_out.append((line.split()[4]).rjust(2))
                # line_out.append((line.split()[5]).rjust(4))
                line_out.append(line[0:60])
                line_out.append(str(round(RMSF_list[atom],2)).rjust(7))
                out_file.write(' '.join(x for x in line_out) + '\n')
                
                line_out_B.append(line[0:60])
                line_out_B.append(str(B_factor_list[atom]).rjust(7))
                out_file_B.write(' '.join(x for x in line_out_B) + '\n')
                
                line_out_Z.append(line[0:60])
                line_out_Z.append(str(Zscore_list[atom]).rjust(7))
                out_file_Z.write(' '.join(x for x in line_out_Z) + '\n')
            except:
                print('WARNING in line: ' + line)      

N_frames_wat = n_to_divide_by

##############################################################################        
# AVERAGING protein&water-h-bonded
temp_for_avg = []
n_to_divide_by = 0
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append([0,0,0])        
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_first, 'r') as in_first, open(out_pdb_prot_wat_avg, 'w') as out_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
            if hbond_info[count] == [True, True]:
                n_to_divide_by += 1
        if hbond_info[count] == [True, True]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                try:                      
                    atom = int(line.split()[1])-1
                    # print(line)
                    temp_for_avg[atom][0] = temp_for_avg[atom][0] + float(line.split()[6])
                    temp_for_avg[atom][1] = temp_for_avg[atom][1] + float(line.split()[7])
                    temp_for_avg[atom][2] = temp_for_avg[atom][2] + float(line.split()[8])
                except:
                    print('WARNING: atom number out of range')     
     
    avg_list = np.array(temp_for_avg)/n_to_divide_by
    for line in in_first:
        line_out = []
        if 'MODEL' in line:
            out_file.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                line_out.append((line.split()[0]).ljust(6))
                line_out.append((line.split()[1]).rjust(5))
                line_out.append('  ' + (line.split()[2]).ljust(4))
                line_out.append((line.split()[3]).ljust(3))
                line_out.append((line.split()[4]).rjust(2))
                line_out.append((line.split()[5]).rjust(4))
                line_out.append('    ' + str(round(avg_list[atom][0],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][1],3)).rjust(8))
                line_out.append(         str(round(avg_list[atom][2],3)).rjust(8))    
                line_out.append((line.split()[9]).rjust(6))
                line_out.append((line.split()[10]).rjust(6))
                out_file.write(''.join(x for x in line_out) + '\n')
            except:
                print('WARNING in line: ' + line)

temp_for_avg = []
with open(out_pdb_first, 'r') as in_first:
    for line in in_first:
        if 'MODEL' in line:
            continue
        elif 'SOL' in line:
            continue
        else:
            temp_for_avg.append(0)
count = -1
with open(out_pdb_prot, 'r') as in_file, open(out_pdb_prot_wat_avg, 'r') as in_avg, open(out_pdb_prot_wat_rmsf, 'w') as out_file, open(out_pdb_prot_wat_Bfactor, 'w') as out_file_B, open(out_pdb_prot_wat_Zscore, 'w') as out_file_Z:
    for line in in_file:
        if 'MODEL' in line:
            count += 1    
        if hbond_info[count] == [True, True]:
            if 'MODEL' not in line:
              if 'TER' not in line:
                  
                try:                  
                    atom = int(line.split()[1])-1
                    dist_vec = np.array(avg_list[atom]) - np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                    dist2 = np.dot(dist_vec,dist_vec)
                    temp_for_avg[atom] = temp_for_avg[atom] + dist2
                except:
                    print('WARNING: atom number out of range')

    avg_list = np.array(temp_for_avg)/n_to_divide_by
    RMSF_list = np.sqrt(avg_list)
    B_factor_list = [round(x,2) for x in (8*np.pi**2*RMSF_list**2/3)]
    B_factor_avg = np.mean(B_factor_list)
    B_factor_std = np.std(B_factor_list)
    Zscore_list = [round(x,2) for x in (B_factor_list - B_factor_avg)/B_factor_std]
    
    for line in in_avg:
        line_out = []
        line_out_B = []
        line_out_Z = []
        if 'MODEL' in line:
            out_file.write(line)
            out_file_B.write(line)
            continue
        elif 'SOL' in line:
            continue
        elif 'TER' in line:
            continue
        else:
            try:
                atom = int(line.split()[1])-1
                # line_out.append((line.split()[0]).ljust(6))
                # line_out.append((line.split()[1]).rjust(5))
                # line_out.append('  ' + (line.split()[2]).ljust(4))
                # line_out.append((line.split()[3]).ljust(3))
                # line_out.append((line.split()[4]).rjust(2))
                # line_out.append((line.split()[5]).rjust(4))
                line_out.append(line[0:60])
                line_out.append(str(round(RMSF_list[atom],2)).rjust(7))
                out_file.write(' '.join(x for x in line_out) + '\n')
                
                line_out_B.append(line[0:60])
                line_out_B.append(str(B_factor_list[atom]).rjust(7))
                out_file_B.write(' '.join(x for x in line_out_B) + '\n')
                
                line_out_Z.append(line[0:60])
                line_out_Z.append(str(Zscore_list[atom]).rjust(7))
                out_file_Z.write(' '.join(x for x in line_out_Z) + '\n')
            except:
                print('WARNING in line: ' + line)      

N_frames_prot_wat = n_to_divide_by

##############################################################################

grid_axis = np.arange(-water_grid_size/2, water_grid_size/2+water_grid_step, water_grid_step)
unit_length = water_grid_unit/2

print('INFO: Generating average water files for ALL.')
#WATER POSITONS:
N_avg = []
with open(out_pdb_all_avg, 'r') as in_file:
    for line in in_file:
        if len(line.split()) > 4:
            if line.split()[1] == str(N_atom):
                N_avg = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])  

grid = []
counter = []
for x in grid_axis:
    for y in grid_axis:
        for z in grid_axis:
            grid.append([x, y, z])
            counter.append(0)

count = -1
with open(out_pdb_wat, 'r') as in_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
        if debug_mode == True:
            if count > 10:
                break
        if 'SOL' in line:
            if line.split()[2] in ['O','OW']: 
                temp_sol = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                temp_vec = temp_sol - N_avg
                temp_dist_list = []
                for i in range(0, len(grid)):
                    if (temp_vec[0] >= grid[i][0] - unit_length) and temp_vec[0] < grid[i][0] + unit_length:
                        if (temp_vec[1] >= grid[i][1] - unit_length) and temp_vec[1] < grid[i][1] + unit_length:
                            if (temp_vec[2] >= grid[i][2] - unit_length) and temp_vec[2] < grid[i][2] + unit_length:
                                counter[i] += 1
                                

count = 0
print(max(counter))
with open(out_pdb_all_wat, 'w') as out_file:
    out_file.write('MODEL 1 \n')
    for sub_counter, sub_grid in zip(counter,grid):
        if sub_counter > max(counter)*water_cutoff:
            occupancy = round(sub_counter/N_frames_none,2)
            count += 1
            line_out = []
            line_out.append(('ATOM').ljust(6))
            line_out.append((str(count)).rjust(5))
            line_out.append('  ' + ('O').ljust(4))
            line_out.append(('SOL').ljust(3))
            line_out.append(('A').rjust(2))
            line_out.append((str(count)).rjust(4))
            line_out.append('    ' + str(round(sub_grid[0]+N_avg[0],3)).rjust(8))
            line_out.append(         str(round(sub_grid[1]+N_avg[1],3)).rjust(8))
            line_out.append(         str(round(sub_grid[2]+N_avg[2],3)).rjust(8))    
            line_out.append((str(occupancy)).rjust(6))
            line_out.append(('0.00').rjust(6))
            # print(''.join(x for x in line_out) + '\n')
            out_file.write(''.join(x for x in line_out) + '\n')       



print('INFO: Generating average water files for NONE.')
#WATER POSITONS:
N_avg = []
with open(out_pdb_none_avg, 'r') as in_file:
    for line in in_file:
        if len(line.split()) > 4:
            if line.split()[1] == str(N_atom):
                N_avg = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])  

grid = []
counter = []
for x in grid_axis:
    for y in grid_axis:
        for z in grid_axis:
            grid.append([x, y, z])
            counter.append(0)

count = -1
with open(out_pdb_wat, 'r') as in_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
        if debug_mode == True:
            if count > 10:
                break
        if 'SOL' in line and hbond_info[count] == [False, False]:
            if line.split()[2] in ['O','OW']: 
                temp_sol = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                temp_vec = temp_sol - N_avg
                temp_dist_list = []
                for i in range(0, len(grid)):
                    if (temp_vec[0] >= grid[i][0] - unit_length) and temp_vec[0] < grid[i][0] + unit_length:
                        if (temp_vec[1] >= grid[i][1] - unit_length) and temp_vec[1] < grid[i][1] + unit_length:
                            if (temp_vec[2] >= grid[i][2] - unit_length) and temp_vec[2] < grid[i][2] + unit_length:
                                counter[i] += 1
                                

count = 0
print(max(counter))
with open(out_pdb_none_wat, 'w') as out_file:
    out_file.write('MODEL 1 \n')
    for sub_counter, sub_grid in zip(counter,grid):
        if sub_counter > max(counter)*water_cutoff:
            occupancy = round(sub_counter/N_frames_none,2)
            count += 1
            line_out = []
            line_out.append(('ATOM').ljust(6))
            line_out.append((str(count)).rjust(5))
            line_out.append('  ' + ('O').ljust(4))
            line_out.append(('SOL').ljust(3))
            line_out.append(('A').rjust(2))
            line_out.append((str(count)).rjust(4))
            line_out.append('    ' + str(round(sub_grid[0]+N_avg[0],3)).rjust(8))
            line_out.append(         str(round(sub_grid[1]+N_avg[1],3)).rjust(8))
            line_out.append(         str(round(sub_grid[2]+N_avg[2],3)).rjust(8))    
            line_out.append((str(occupancy)).rjust(6))
            line_out.append(('0.00').rjust(6))
            # print(''.join(x for x in line_out) + '\n')
            out_file.write(''.join(x for x in line_out) + '\n')         
    
    
    
    
    
print('INFO: Generating average water files for PROT.')    
N_avg = []
with open(out_pdb_prot_avg, 'r') as in_file:
    for line in in_file:
        if len(line.split()) > 4:
            if line.split()[1] == str(N_atom):
                N_avg = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])  

grid = []
counter = []
for x in grid_axis:
    for y in grid_axis:
        for z in grid_axis:
            grid.append([x, y, z])
            counter.append(0)
              
count = -1
with open(out_pdb_wat, 'r') as in_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
        if debug_mode == True:
            if count > 10:
                break
        if 'SOL' in line and hbond_info[count] == [True, False]:
            if line.split()[2] in ['O','OW']: 
                temp_sol = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                temp_vec = temp_sol - N_avg
                temp_dist_list = []
                for i in range(0, len(grid)):
                    if (temp_vec[0] >= grid[i][0] - unit_length) and temp_vec[0] < grid[i][0] + unit_length:
                        if (temp_vec[1] >= grid[i][1] - unit_length) and temp_vec[1] < grid[i][1] + unit_length:
                            if (temp_vec[2] >= grid[i][2] - unit_length) and temp_vec[2] < grid[i][2] + unit_length:
                                counter[i] += 1
                                

count = 0
print(max(counter))
with open(out_pdb_prot_wat, 'w') as out_file:
    out_file.write('MODEL 1 \n')
    for sub_counter, sub_grid in zip(counter,grid):
        if sub_counter > max(counter)*water_cutoff:
            occupancy = round(sub_counter/N_frames_prot,2)
            count += 1
            line_out = []
            line_out.append(('ATOM').ljust(6))
            line_out.append((str(count)).rjust(5))
            line_out.append('  ' + ('O').ljust(4))
            line_out.append(('SOL').ljust(3))
            line_out.append(('A').rjust(2))
            line_out.append((str(count)).rjust(4))
            line_out.append('    ' + str(round(sub_grid[0]+N_avg[0],3)).rjust(8))
            line_out.append(         str(round(sub_grid[1]+N_avg[1],3)).rjust(8))
            line_out.append(         str(round(sub_grid[2]+N_avg[2],3)).rjust(8))    
            line_out.append((str(occupancy)).rjust(6))
            line_out.append(('0.00').rjust(6))
            # print(''.join(x for x in line_out) + '\n')
            out_file.write(''.join(x for x in line_out) + '\n')      

            

            
            
print('INFO: Generating average water files for WATER.')            
N_avg = []
with open(out_pdb_wat_avg, 'r') as in_file:
    for line in in_file:
        if len(line.split()) > 4:
            if line.split()[1] == str(N_atom):
                N_avg = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])  

grid = []
counter = []
for x in grid_axis:
    for y in grid_axis:
        for z in grid_axis:
            grid.append([x, y, z])
            counter.append(0)
              
count = -1
with open(out_pdb_wat, 'r') as in_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
        if debug_mode == True:
            if count > 10:
                break
        if 'SOL' in line and hbond_info[count] == [False, True]:
            if line.split()[2] in ['O','OW']: 
                temp_sol = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                temp_vec = temp_sol - N_avg
                temp_dist_list = []
                for i in range(0, len(grid)):
                    if (temp_vec[0] >= grid[i][0] - unit_length) and temp_vec[0] < grid[i][0] + unit_length:
                        if (temp_vec[1] >= grid[i][1] - unit_length) and temp_vec[1] < grid[i][1] + unit_length:
                            if (temp_vec[2] >= grid[i][2] - unit_length) and temp_vec[2] < grid[i][2] + unit_length:
                                counter[i] += 1
                                

count = 0
print(max(counter))
with open(out_pdb_wat_wat, 'w') as out_file:
    out_file.write('MODEL 1 \n')
    for sub_counter, sub_grid in zip(counter,grid):
        if sub_counter > max(counter)*water_cutoff:
            occupancy = round(sub_counter/N_frames_wat,2)
            count += 1
            line_out = []
            line_out.append(('ATOM').ljust(6))
            line_out.append((str(count)).rjust(5))
            line_out.append('  ' + ('O').ljust(4))
            line_out.append(('SOL').ljust(3))
            line_out.append(('A').rjust(2))
            line_out.append((str(count)).rjust(4))
            line_out.append('    ' + str(round(sub_grid[0]+N_avg[0],3)).rjust(8))
            line_out.append(         str(round(sub_grid[1]+N_avg[1],3)).rjust(8))
            line_out.append(         str(round(sub_grid[2]+N_avg[2],3)).rjust(8))    
            line_out.append((str(occupancy)).rjust(6))
            line_out.append(('0.00').rjust(6))
            # print(''.join(x for x in line_out) + '\n')
            out_file.write(''.join(x for x in line_out) + '\n')      




print('INFO: Generating average water files for PROT and WAT.')
N_avg = []
with open(out_pdb_prot_wat_avg, 'r') as in_file:
    for line in in_file:
        if len(line.split()) > 4:
            if line.split()[1] == str(N_atom):
                N_avg = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])  

grid = []
counter = []
for x in grid_axis:
    for y in grid_axis:
        for z in grid_axis:
            grid.append([x, y, z])
            counter.append(0)
  
plus_minus = water_grid_step/2            
count = -1
with open(out_pdb_wat, 'r') as in_file:
    for line in in_file:
        if 'MODEL' in line:
            count += 1
        if debug_mode == True:
            if count > 10:
                break
        if 'SOL' in line and hbond_info[count] == [True, True]:
            if line.split()[2] in ['O','OW']: 
                temp_sol = np.array([float(line.split()[6]), float(line.split()[7]), float(line.split()[8])])
                temp_vec = temp_sol - N_avg
                temp_dist_list = []
                for i in range(0, len(grid)):
                    if (temp_vec[0] >= grid[i][0] - unit_length) and temp_vec[0] < grid[i][0] + unit_length:
                        if (temp_vec[1] >= grid[i][1] - unit_length) and temp_vec[1] < grid[i][1] + unit_length:
                            if (temp_vec[2] >= grid[i][2] - unit_length) and temp_vec[2] < grid[i][2] + unit_length:
                                counter[i] += 1
                                

count = 0
print(max(counter))
with open(out_pdb_prot_wat_wat, 'w') as out_file:
    out_file.write('MODEL 1 \n')
    for sub_counter, sub_grid in zip(counter,grid):
        if sub_counter > max(counter)*water_cutoff:
            occupancy = round(sub_counter/N_frames_prot_wat,2)
            count += 1
            line_out = []
            line_out.append(('ATOM').ljust(6))
            line_out.append((str(count)).rjust(5))
            line_out.append('  ' + ('O').ljust(4))
            line_out.append(('SOL').ljust(3))
            line_out.append(('A').rjust(2))
            line_out.append((str(count)).rjust(4))
            line_out.append('    ' + str(round(sub_grid[0]+N_avg[0],3)).rjust(8))
            line_out.append(         str(round(sub_grid[1]+N_avg[1],3)).rjust(8))
            line_out.append(         str(round(sub_grid[2]+N_avg[2],3)).rjust(8))    
            line_out.append((str(occupancy)).rjust(6))
            line_out.append(('0.00').rjust(6))
            # print(''.join(x for x in line_out) + '\n')
            out_file.write(''.join(x for x in line_out) + '\n')                     