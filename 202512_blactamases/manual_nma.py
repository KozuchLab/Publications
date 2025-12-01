import numpy as np
from subprocess import run
import os
import shutil

def isint(string):
    try:
        int(string)
    except:
        return False
    return True

def isfloat(string):
    try:
        float(string)
    except:
        return False
    return True

def isatomname(string):
    if string in ['H','C','N','O','F','S','Cl']:
        return True
    else:
        return True

def read_tink_xyz_as_template(file):
    tink_lines = []
    with open(file, 'r') as in_file:
        for line in in_file:
            tink_lines.append(line)
    return tink_lines

def read_tink_xyz_coords(file):
    tink_coords = []
    tink_masses = []
    
    atom2mass = {'H' : 1.00784,
                 'C' : 12.011,
                 'N' : 14.0067,
                 'O' : 15.999,
                 'F' : 18.9984,
                 'S' : 32.065,
                 'Cl' : 35.453}
    
    with open(file, 'r') as in_file:       
        for line in in_file:    
            if len(line.split()) > 2:
                ele1 = line.split()[0]
                ele2 = line.split()[1]
                if isint(ele1) and isatomname(ele2):
                    tink_coords.append([float(x) for x in line.split()[2:5]])
                    tink_masses.append(atom2mass[line.split()[1]])
    return tink_coords, tink_masses

def xyz_linear(xyz):

    np_xyz = np.array(xyz).flatten()

    return np_xyz.tolist()


def xyz_Nby3(xyz):

    length = int(len(xyz)/3)

    new_xyz = []
    for i in range(0, length):
        new_xyz.append([xyz[3*i], xyz[3*i+1], xyz[3*i+2]])
        
    return new_xyz

def write_tink_xyz(filename,template,coords):
    
    num_header_lines = int(len(template) - len(coords))

    with open(filename, 'w') as out_file:
        for line in template[0:num_header_lines]:
            out_file.write(line)
        for line1, line2 in zip(template[num_header_lines:],coords):
            ele1 = line1.split()[0]
            ele2 = line1.split()[1]
            ele3 = str(line2[0])
            ele4 = str(line2[1])
            ele5 = str(line2[2])
            ele6 = ' '.join(x for x in line1.split()[5:])
            string = ' '.join(x for x in [ele1, ele2, ele3, ele4, ele5, ele6])
            out_file.write(string + '\n')
            
def write_tink_xyz_for_nma(filename_E,filename_M,template,in_coords,in_masses):
    
    dr = 0.01   
    num_header_lines = int(len(template) - len(in_coords))
    coords_lin = xyz_linear(in_coords)
    idx_jdx_dr_dr_m_m = []
    masses = []
    for m in in_masses:
        for i in range(3):
            masses.append(m)

    with open(filename_E, 'w') as out_file_E,open(filename_M, 'w') as out_file_M:
        for idx in range(len(coords_lin)):
            for jdx in  range(len(coords_lin)):
                if idx <= jdx:
                    displ_coords_lin = np.array(coords_lin.copy())
                    if idx == jdx:
                        displ_coords_lin[idx] = displ_coords_lin[idx] + dr
                    else:
                        displ_coords_lin[idx] = displ_coords_lin[idx] + dr
                        displ_coords_lin[jdx] = displ_coords_lin[jdx] + dr
                    idx_jdx_dr_dr_m_m.append([idx,jdx,dr,dr,masses[idx],masses[jdx]])
                    use_coords = xyz_Nby3(displ_coords_lin)
                    for line in template[0:num_header_lines]:
                        out_file_E.write(line)
                        if idx == jdx:
                            out_file_M.write(line)
                    for line1, line2 in zip(template[num_header_lines:],use_coords):
                        ele1 = line1.split()[0]
                        ele2 = line1.split()[1]
                        ele3 = str(line2[0])
                        ele4 = str(line2[1])
                        ele5 = str(line2[2])
                        ele6 = ' '.join(x for x in line1.split()[5:])
                        string = ' '.join(x for x in [ele1, ele2, ele3, ele4, ele5, ele6])
                        out_file_E.write(string + '\n')
                        if idx == jdx:
                            out_file_M.write(string + '\n')
        for idx in range(len(coords_lin)):
            for jdx in  range(len(coords_lin)):
                if idx <= jdx:
                    displ_coords_lin = np.array(coords_lin.copy())
                    if idx == jdx:
                        displ_coords_lin[idx] = displ_coords_lin[idx]
                    else:
                        displ_coords_lin[idx] = displ_coords_lin[idx] + dr
                        displ_coords_lin[jdx] = displ_coords_lin[jdx] - dr
                    # idx_jdx_dr_dr.append([idx,jdx,dr,-dr,masses[idx],masses[jdx]])
                    use_coords = xyz_Nby3(displ_coords_lin)
                    for line in template[0:num_header_lines]:
                        out_file_E.write(line)
                    for line1, line2 in zip(template[num_header_lines:],use_coords):
                        ele1 = line1.split()[0]
                        ele2 = line1.split()[1]
                        ele3 = str(line2[0])
                        ele4 = str(line2[1])
                        ele5 = str(line2[2])
                        ele6 = ' '.join(x for x in line1.split()[5:])
                        string = ' '.join(x for x in [ele1, ele2, ele3, ele4, ele5, ele6])
                        out_file_E.write(string + '\n')
        for idx in range(len(coords_lin)):
            for jdx in  range(len(coords_lin)):
                if idx <= jdx:
                    displ_coords_lin = np.array(coords_lin.copy())
                    if idx == jdx:
                        displ_coords_lin[idx] = displ_coords_lin[idx]
                    else:
                        displ_coords_lin[idx] = displ_coords_lin[idx] - dr
                        displ_coords_lin[jdx] = displ_coords_lin[jdx] + dr
                    # idx_jdx_dr_dr.append([idx,jdx,-dr,dr,masses[idx],masses[jdx]])
                    use_coords = xyz_Nby3(displ_coords_lin)
                    for line in template[0:num_header_lines]:
                        out_file_E.write(line)
                    for line1, line2 in zip(template[num_header_lines:],use_coords):
                        ele1 = line1.split()[0]
                        ele2 = line1.split()[1]
                        ele3 = str(line2[0])
                        ele4 = str(line2[1])
                        ele5 = str(line2[2])
                        ele6 = ' '.join(x for x in line1.split()[5:])
                        string = ' '.join(x for x in [ele1, ele2, ele3, ele4, ele5, ele6])
                        out_file_E.write(string + '\n')
        for idx in range(len(coords_lin)):
            for jdx in  range(len(coords_lin)):
                if idx <= jdx:
                    displ_coords_lin = np.array(coords_lin.copy())            
                    if idx == jdx:
                        displ_coords_lin[idx] = displ_coords_lin[idx] - dr
                    else:
                        displ_coords_lin[idx] = displ_coords_lin[idx] - dr
                        displ_coords_lin[jdx] = displ_coords_lin[jdx] - dr
                    # idx_jdx_dr_dr.append([idx,jdx,-dr,-dr,masses[idx],masses[jdx]])
                    use_coords = xyz_Nby3(displ_coords_lin)
                    for line in template[0:num_header_lines]:
                        out_file_E.write(line)
                        if idx == jdx:
                            out_file_M.write(line)
                    for line1, line2 in zip(template[num_header_lines:],use_coords):
                        ele1 = line1.split()[0]
                        ele2 = line1.split()[1]
                        ele3 = str(line2[0])
                        ele4 = str(line2[1])
                        ele5 = str(line2[2])
                        ele6 = ' '.join(x for x in line1.split()[5:])
                        string = ' '.join(x for x in [ele1, ele2, ele3, ele4, ele5, ele6])
                        out_file_E.write(string + '\n')
                        if idx == jdx:
                            out_file_M.write(string + '\n')
    return idx_jdx_dr_dr_m_m

def vT_M_v(vT, M, v):
    max_idx = len(M)
    
    pre_M = np.zeros((max_idx, max_idx, 3))

    
    for i in range(0, max_idx):
        for j in range(0, max_idx):
            for k in range(0, max_idx):
                pre_M[i][j] += M[i][k]*v[k][j]
                
    new_M = np.zeros((max_idx, max_idx, 3))
    
    for i in range(0, max_idx):
        for j in range(0, max_idx):
            for k in range(0, max_idx):
                new_M[i][j] += vT[i][k]*pre_M[k][j]
                
    return new_M
         
def analyze(filename_E,filename_M,keyfile,idx_jdx_dr_dr_m_m,in_masses):
    num_atoms = []
    with open(filename_E, 'r') as in_file:
        for line in in_file:
            num_atoms=int(line)
            break
        
    masses = []
    for m in in_masses:
        for i in range(3):
            masses.append(m)
    
    run('analyze %s -k %s E > analyze_E.out' %(filename_E, keyfile),shell=True)
    nrgs = []
    with open('analyze_E.out', 'r') as in_file:
        for line in in_file:
            if 'Total Potential Energy' in line:
                nrgs.append(float(line.split()[-2]))
                
    nrgs = np.array(nrgs).reshape((4,-1))  
    ddnrgs_dqdq = []
    idx_jdx = []
    
    for nrg1,nrg2,nrg3,nrg4,ijrrmm in zip(nrgs[0],nrgs[1],nrgs[2],nrgs[3],idx_jdx_dr_dr_m_m):
        if ijrrmm[0] == ijrrmm[1]:
            idx_jdx.append([ijrrmm[0], ijrrmm[1]])
            ddEdqdq_num = nrg1-nrg2-nrg3+nrg4
            ddEdqdq_den = np.sqrt(ijrrmm[4])*np.sqrt(ijrrmm[5])*ijrrmm[2]*ijrrmm[3]
            ddEdqdq = ddEdqdq_num/ddEdqdq_den
            ddnrgs_dqdq.append(ddEdqdq)
            # print(str(ijrrmm[0]) + ' ' + str(ijrrmm[1]) + ' ' + str(ddEdqdq))
        else:
            idx_jdx.append([ijrrmm[0], ijrrmm[1]])
            ddEdqdq_num = nrg1-nrg2-nrg3+nrg4
            ddEdqdq_den = np.sqrt(ijrrmm[4])*np.sqrt(ijrrmm[5])*4*ijrrmm[2]*ijrrmm[3]
            ddEdqdq = ddEdqdq_num/ddEdqdq_den
            ddnrgs_dqdq.append(ddEdqdq)
            # print(str(ijrrmm[0]) + ' ' + str(ijrrmm[1]) + ' ' + str(ddEdqdq))
        
    '''Using 3*N+1 x 3*N+1 matrix to be able to calculate dipole derivatives and frequencies'''
    num_coords = 3*num_atoms
    hessian_matrix = np.zeros((num_coords+1,num_coords+1))
        
    for ij,dn in zip(idx_jdx,ddnrgs_dqdq):
        idx = ij[0]+1
        jdx = ij[1]+1
        hessian_matrix[idx][jdx] = dn
        hessian_matrix[jdx][idx] = dn
        
    eigenvalues, eigenvectors = np.linalg.eig(hessian_matrix)
    
    
    run('analyze %s -k %s M > analyze_M.out' %(filename_M, keyfile),shell=True)
    dipoles = []
    with open('analyze_M.out', 'r') as in_file:
        for line in in_file:
            if 'Dipole X,Y,Z-Components' in line:
                dipoles.append([float(x) for x in line.split()[-3:]])      
                
    dipoles_p = dipoles[0:int(len(dipoles)/2)]
    dipoles_m = dipoles[int(len(dipoles)/2):]
    dr = 0.01   
    dipders = np.subtract(dipoles_p,dipoles_m)/2/dr           # in D/Ang
    
    dipder_matrix = np.zeros((len(dipders)+1,len(dipders)+1, 3))
    for idx in range(len(dipders)):
        dipder_matrix[idx+1][0] = dipders[idx]
        dipder_matrix[0][idx+1] = dipders[idx]
        
    dipders = vT_M_v(np.array(eigenvectors).transpose(), np.array(dipder_matrix), np.array(eigenvectors))
       
    '''Cropping the results by removing the initial zero force constant'''
    eigenvalues = eigenvalues[0:-1]
    eigenvectors = eigenvectors[1:]
    displacements = np.array(eigenvectors).transpose()
    displacements = displacements[0:-1]
    dipders = dipders[-1]
    dipders = dipders[0:-1]                       
    
    sqrt_masses = np.sqrt(masses)
    red_masses = []
    cart_displs = []
    
    for ev in displacements:
        carts = np.divide(ev,sqrt_masses)
        norm_ev = np.sqrt(np.sum(ev**2))
        norm_cart = np.sqrt(np.sum(carts**2))
        red_mass = (norm_ev/norm_cart)**2
        red_masses.append(red_mass)
        carts = np.divide(carts,norm_cart)
        cart_displs.append(xyz_Nby3(carts))
        
    dipders = dipders                                          # in D/Ang
    dipders = dipders * 0.52917721092                          # in D/Bohr
    dipders = dipders * 0.39343                                # in e
    new_dipders = []
    for dipder,rm in zip(dipders,red_masses):
        new_dipders.append(dipder/np.sqrt(rm))
    dipders = np.array(new_dipders)                            # in e/sqrt(amu)
    dipders = dipders * 31.223068                              # in sqrt(km/mol)
    dipders.tolist()
        
    freqs = []
    for ev, rm in zip(eigenvalues,red_masses):
        
        # A few comments:
        # ev is omega^2 in units of kcal mol^-1 Ang^-2 amu^-1
        # wn = omega/(2*pi*c)
        # expressing wn in cm^-1 requires: wn = sqrt(ev * 11792.0831650134)
        #                              or: wn = 108.591358611141 * sqrt(ev)
        # first bringing ev to units of s^-2:
        # ev                                     # in kcal mol^-1 Ang^-2 amu^-1
        # ev = ev * 4184000                      # in J Ang^-2 kg^-1 
        # ev = ev / 1e-20                        # in J m^-2 kg^-1 = s^-2
        # now calculating wn^2 from this
        # ev = ev / (2*np.pi*2.99792458)**2      # in m^-2
        # ev = ev / 1e-4                         # in cm^-2
        # wn = np.sqrt(ev)                       # in cm^-1
        # alltogether:

        if ev >= 0:
            freq = np.sqrt(ev * 11792.0831650134)
            # print(str(ev) + ' ' + str(np.sqrt(ev)) + ' ' + str(freq))
        else:
            freq = -np.sqrt(-ev * 11792.0831650134)
            # print(str(ev) + ' ' + str(-np.sqrt(-ev)) + ' ' + str(freq))
        freqs.append(freq)
        
    dipders = np.array(dipders).transpose()
    
    freq_redmass_displs = np.array([freqs, dipders[0], dipders[1], dipders[2], red_masses, cart_displs],dtype=object).transpose()
    freq_redmass_displs.tolist()
    freq_redmass_displs = sorted(freq_redmass_displs, key=lambda x: x[0])
    # for ele in freq_redmass_displs:
    #     f = ele[0]
    #     m1 = ele[1]
    #     m2 = ele[2]
    #     m3 = ele[3]
    #     print(str(round(f,2)) + ' ' + str(round(m1,4)) + ' ' + str(round(m2,4)) + ' ' + str(round(m3,4)))

    return freq_redmass_displs   

def print_output_in_tinkerformat(freq_redmass_displs, template, coords, input_xyz, output_file):
    
    with open(output_file,'w') as out_file:
        out_file.write('    Mode     Freq [cm-1]    DipDer_x DipDer_y DipDer_z [sqrt(km/mol)]      Reduced mass [amu] \n')
        i = 0
        for ele in freq_redmass_displs:
            i = i + 1
            out_file.write( str(i).rjust(4) )
            out_file.write( str(round(ele[0],4)).rjust(12) )
            out_file.write( str(round(ele[1],4)).rjust(12) )
            out_file.write( str(round(ele[2],4)).rjust(12) )
            out_file.write( str(round(ele[3],4)).rjust(12) )
            out_file.write( str(round(ele[4],4)).rjust(12) )
            out_file.write('\n')
        out_file.write('\n')
        
        i = 0
        for ele in freq_redmass_displs:
            i = i + 1
            out_file.write(' Vibrational Normal Mode  '+str(i).rjust(4)+' with Frequency'+str(round(ele[0],4)).rjust(12)+' cm-1 \n')
            out_file.write('\n')
            out_file.write('     Atom     Delta X     Delta Y     Delta Z \n')
            out_file.write('\n')
            j = 0
            for coord in ele[5]:
                j = j + 1
                out_file.write(str(round(j,6)).rjust(9))
                for x in coord:
                    out_file.write(str(round(x,6)).rjust(12))
                out_file.write('\n')
            out_file.write('\n')
            
    
    # write_tink_xyz(input_xyz,template,coords)
            
    
            
def manual_nma(filenames):
    
    shutil.copy(filenames['tink_min_aligned_xyz'],filenames['tink_vib_xyz']) 
    input_xyz = filenames['tink_vib_xyz']
    input_key = filenames['tink_temp_key']
    output_file = filenames['tink_vib_out']
                
    template           = read_tink_xyz_as_template(input_xyz)
    coords, masses     = read_tink_xyz_coords(input_xyz)
    idx_jdx_dr_dr_m_m  = write_tink_xyz_for_nma('for_analyze_E.xyz','for_analyze_M.xyz',template,coords,masses)
    freq_redmass_displs = analyze('for_analyze_E.xyz','for_analyze_M.xyz',input_key,idx_jdx_dr_dr_m_m,masses)
    print_output_in_tinkerformat(freq_redmass_displs, template, coords, input_xyz, output_file)
    os.remove('for_analyze_E.xyz')
    os.remove('for_analyze_M.xyz')
    os.remove('analyze_E.out')

# input_xyz = 'tink_NMA_min_aligned.xyz'
# input_key = 'tink_NMA_temp.key'
# output_file = 'manual_nma.out'

# manual_nma(input_xyz, input_key, output_file)
























